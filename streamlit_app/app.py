import streamlit as st
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from shapely.wkt import loads
import folium
from folium import LayerControl, FeatureGroup
from folium.raster_layers import ImageOverlay
from streamlit_folium import st_folium
import rasterio
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from geopy.distance import geodesic

st.set_page_config(layout="wide")
st.title("🚗 EV Charging Station Estimator - South Delhi")

# === Load Data ===
points_df = pd.read_csv("../output/scored_points.csv")
points_df["geometry"] = points_df["geometry"].apply(loads)
points_gdf = gpd.GeoDataFrame(points_df, geometry="geometry", crs="EPSG:32643")
optimal_sites = gpd.read_file("../output/optimal_sites.geojson").to_crs(epsg=4326)
points_gdf = points_gdf.to_crs(epsg=4326)

# Load vector layers
layer_path = "/home/winvin/Documents/geota_final_layers"
layer_files = {
    "Usable Area": "useable_area.gpkg",
    "Existing Ev Stations": "ev_stations.gpkg",
    "Hospitals": "hospitals.gpkg",
    "Universities": "university.gpkg",
    "Markets": "market_places.gpkg",
    "Amenity": "landuse_amenity.gpkg",
    "Forests": "landuse_forest.gpkg",
    "Historic Monuments": "historic_monuments.gpkg",
    "Protected Areas": "protected_areas.gpkg",
    "Waterways": "waterways.gpkg",
    "Power Lines": "powerlines.gpkg",
    "Roads": "roads.gpkg"
}
loaded_layers = {}
for name, file in layer_files.items():
    try:
        gdf = gpd.read_file(f"{layer_path}/{file}").to_crs(epsg=4326)
        loaded_layers[name] = gdf
    except:
        st.warning(f"Failed to load layer: {name}")

# Load South Delhi boundary
south_delhi_boundary = gpd.read_file(f"{layer_path}/south_delhi_layer.gpkg").to_crs(epsg=4326)

# === Tabs ===
tab1, tab2 = st.tabs(["🗺️ Map & Layers", "📍 Find Nearest EV Station"])

with tab1:
    st.sidebar.title("🗂 Layer Options")
    selected_layers = [name for name in loaded_layers if st.sidebar.checkbox(name, value=False)]
    show_south_delhi = st.sidebar.checkbox("South Delhi Boundary", value=True)
    show_pop_raster = st.sidebar.checkbox("Population Density (Raster)", value=True)
    operation = st.sidebar.selectbox("Spatial Operation on Selected Layers", ["None", "Union", "Intersection", "Difference"])
    process_layers = st.sidebar.button("🔁 Process Layers")

    m = folium.Map(location=[28.6, 77.2], zoom_start=12)

    # --- Add Selected Vector Layers ---
    for name in selected_layers:
        gdf = loaded_layers[name]
        layer = FeatureGroup(name=name)
        for _, row in gdf.iterrows():
            if row.geometry.geom_type in ['Polygon', 'MultiPolygon']:
                folium.GeoJson(row.geometry, name=name).add_to(layer)
            elif row.geometry.geom_type == 'Point':
                folium.CircleMarker(
                    location=[row.geometry.y, row.geometry.x],
                    radius=4, color="blue", fill=True, fill_opacity=0.6
                ).add_to(layer)
        m.add_child(layer)

    # --- Add South Delhi Boundary ---
    if show_south_delhi:
        boundary_layer = FeatureGroup(name="South Delhi Boundary")
        folium.GeoJson(
            south_delhi_boundary.geometry,
            style_function=lambda x: {
                "fillColor": "#00000000",
                "color": "black",
                "weight": 2
            }
        ).add_to(boundary_layer)
        m.add_child(boundary_layer)

    # --- Add Population Raster Layer ---
    if show_pop_raster:
        try:
            with rasterio.open(f"{layer_path}/south_delhi_pop.tif") as src:
                img_array = src.read(1)
                bounds = src.bounds
                img_array = np.ma.masked_equal(img_array, src.nodata)

                norm = Normalize(vmin=np.min(img_array), vmax=np.max(img_array))
                cmap = cm.get_cmap("RdYlGn_r")
                colored = cmap(norm(img_array))
                colored = (colored[:, :, :3] * 255).astype(np.uint8)
                plt.imsave("temp_overlay.png", colored)

                image_overlay = ImageOverlay(
                    name="Population Density (Raster)",
                    image="temp_overlay.png",
                    bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
                    opacity=0.6,
                    interactive=True,
                    cross_origin=False
                )
                image_overlay.add_to(m)
        except Exception as e:
            st.error(f"Error loading raster: {e}")

    # --- Show Top EV Estimation Points ---
    top_sites = points_gdf.sort_values("score", ascending=False).head(10)
    ev_layer = FeatureGroup(name="Top 10 EV Sites")
    for _, row in top_sites.iterrows():
        folium.Marker(
            location=[row.geometry.y, row.geometry.x],
            icon=folium.Icon(color="green", icon="bolt", prefix='fa'),
            popup=f"Score: {row['score']}"
        ).add_to(ev_layer)
    m.add_child(ev_layer)

    # --- Spatial Operations ---
    if process_layers and len(selected_layers) >= 2 and operation != "None":
        gdfs = [loaded_layers[name].to_crs(epsg=32643) for name in selected_layers]
        result = gdfs[0]
        for gdf in gdfs[1:]:
            if operation == "Union":
                result = result.overlay(gdf, how='union')
            elif operation == "Intersection":
                result = result.overlay(gdf, how='intersection')
            elif operation == "Difference":
                result = result.overlay(gdf, how='difference')
        result = result.to_crs(epsg=4326)

        operation_layer = FeatureGroup(name=f"{operation} of Selected Layers")
        for _, row in result.iterrows():
            if row.geometry.geom_type in ['Polygon', 'MultiPolygon']:
                folium.GeoJson(row.geometry, style_function=lambda x: {"fillColor": "orange", "color": "orange", "weight": 1}).add_to(operation_layer)
        m.add_child(operation_layer)

    folium.LayerControl(collapsed=False).add_to(m)
    st_data = st_folium(m, width=1100, height=650)

with tab2:
    st.subheader("📍 Find Nearest EV Station")

    lat = st.number_input("Latitude", value=28.61, format="%.6f")
    lon = st.number_input("Longitude", value=77.21, format="%.6f")
    user_point = Point(lon, lat)

    map2 = folium.Map(location=[lat, lon], zoom_start=13)
    folium.Marker([lat, lon], tooltip="Your Location", icon=folium.Icon(color='blue')).add_to(map2)

    # Compute nearest EV station
    optimal_sites["dist"] = optimal_sites.geometry.apply(lambda geom: geodesic((lat, lon), (geom.y, geom.x)).meters)
    nearest = optimal_sites.sort_values("dist").iloc[0]
    nearest_coords = [nearest.geometry.y, nearest.geometry.x]

    folium.Marker(
        location=nearest_coords,
        icon=folium.Icon(color="red", icon="bolt", prefix='fa'),
        popup=f"Nearest EV Station\nDistance: {nearest['dist']:.2f} meters\nScore: {nearest['score']:.2f}"
    ).add_to(map2)

    folium.PolyLine([ [lat, lon], nearest_coords ], color="blue", weight=2.5, opacity=0.7).add_to(map2)

    st.write(f"🚶 **Distance to nearest EV Station**: `{nearest['dist']:.2f} meters`")
    st_folium(map2, width=1000, height=600)
