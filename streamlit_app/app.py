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
import os
from pathlib import Path

# Set up paths properly
base_dir = Path(os.path.dirname(os.path.abspath(__file__)))
output_dir = base_dir.parent / "output"
layer_path = base_dir.parent / "geota_final_layers"

st.set_page_config(layout="wide")
st.title("🚗 EV Infrastructure Estimator - South Delhi")

# === Load Data with error handling ===
try:
    # Load points data for charging stations
    points_df = pd.read_csv(output_dir / "scored_points.csv")
    # Parse geometry column and ensure it's Point
    points_df["geometry"] = points_df["geometry"].apply(lambda x: loads(x) if isinstance(x, str) else x)
    # Verify geometry types
    points_gdf = gpd.GeoDataFrame(points_df, geometry="geometry", crs="EPSG:32643")
    non_points = points_gdf[~points_gdf.geometry.geom_type.isin(['Point'])]
    if not non_points.empty:
        st.warning(f"Found {len(non_points)} non-Point geometries in points_gdf: {non_points.geometry.geom_type.unique()}")
        points_gdf = points_gdf[points_gdf.geometry.geom_type == 'Point']  # Keep only Points
    
    optimal_sites = gpd.read_file(output_dir / "optimal_sites.geojson").to_crs(epsg=4326)
    points_gdf = points_gdf.to_crs(epsg=4326)

    # Showroom data
    showroom_df = pd.read_csv(output_dir / "scored_showroom_points.csv")
    showroom_df["geometry"] = showroom_df["geometry"].apply(lambda x: loads(x) if isinstance(x, str) else x)
    showroom_gdf = gpd.GeoDataFrame(showroom_df, geometry="geometry", crs="EPSG:32643")
    non_points_showroom = showroom_gdf[~showroom_gdf.geometry.geom_type.isin(['Point'])]
    if not non_points_showroom.empty:
        st.warning(f"Found {len(non_points_showroom)} non-Point geometries in showroom_gdf: {non_points_showroom.geometry.geom_type.unique()}")
        showroom_gdf = showroom_gdf[showroom_gdf.geometry.geom_type == 'Point']  # Keep only Points
    
    optimal_showrooms = gpd.read_file(output_dir / "optimal_showroom_sites.geojson").to_crs(epsg=4326)
    showroom_gdf = showroom_gdf.to_crs(epsg=4326)
except Exception as e:
    st.error(f"Error loading data: {e}")
    st.info("Please make sure all data files exist in the correct directories")
    st.stop()

# Load vector layers
layer_files = {
    "Usable Area": "useable_area.gpkg",
    "Existing EV Stations": "ev_stations.gpkg",
    "Existing Showrooms": "car_showrooms.gpkg",
    "Hospitals": "hospitals.gpkg",
    "Universities": "university.gpkg",
    "Markets": "market_places.gpkg",
    "Amenity": "landuse_amenity.gpkg",
    "Forests": "landuse_forest.gpkg",
    "Historic Monuments": "historic_monuments.gpkg",
    "Protected Areas": "protected_areas.gpkg",
    "Waterways": "waterways.gpkg",
    "Power Lines": "powerlines.gpkg",
    "Roads": "roads.gpkg",
    "Major Roads": "major_roads_only.gpkg",
    "Main Roads": "main_roads.gpkg",
    "Petrol Pumps": "petrol_pumps.gpkg",
    "Buffered Petrol Pumps": "buffered_petrol_pumps.gpkg",
    "Police Stations": "police_stations.gpkg",
    "Buildings": "buildings.gpkg",
    "Bus Stops": "bus_stops.gpkg",
    "Favorable Regions": "favourable_regions.gpkg",
    "Union Regions": "union_regions.gpkg",
    "Residential Areas": "landuse_residential.gpkg",
    "Commercial Areas": "landuse_commercial.gpkg",
    "Industrial Areas": "landuse_industrial.gpkg",
    "Places": "places.gpkg"
}

loaded_layers = {}
for name, file in layer_files.items():
    try:
        gdf = gpd.read_file(layer_path / file).to_crs(epsg=4326)
        loaded_layers[name] = gdf
    except Exception as e:
        st.warning(f"Failed to load layer: {name} - {str(e)}")

# Load South Delhi boundary
try:
    south_delhi_boundary = gpd.read_file(layer_path / "south_delhi_layer.gpkg").to_crs(epsg=4326)
except:
    try:
        south_delhi_boundary = gpd.read_file(layer_path / "south_delhi_layer.shp").to_crs(epsg=4326)
    except Exception as e:
        st.error(f"Failed to load South Delhi boundary: {e}")
        south_delhi_boundary = None

# === Tabs ===
tab1, tab2, tab3, tab4 = st.tabs(["🗺️ Map & Layers", "📍 Optimal EV Stations", "🏪 Optimal EV Showrooms", "📊 Analysis"])

with tab1:
    col1, col2 = st.columns([3, 7])
    
    with col1:
        st.subheader("Layer Options")
        
        layer_categories = {
            "🚗 Transportation": ["Roads", "Major Roads", "Main Roads", "Existing EV Stations", "Existing Showrooms", "Petrol Pumps"],
            "🏙️ Land Use": ["Residential Areas", "Commercial Areas", "Industrial Areas", "Amenity", "Forests"],
            "🏢 Points of Interest": ["Hospitals", "Universities", "Markets", "Police Stations", "Bus Stops"],
            "⚠️ Restricted Areas": ["Protected Areas", "Historic Monuments", "Waterways"],
            "🗺️ Planning": ["Usable Area", "Favorable Regions", "Union Regions"]
        }
        
        selected_layers = []
        for category, layers in layer_categories.items():
            with st.expander(category, expanded=False):
                for layer in layers:
                    if layer in loaded_layers and st.checkbox(layer, value=False, key=f"layer_{layer}"):
                        selected_layers.append(layer)
        
        show_south_delhi = st.checkbox("South Delhi Boundary", value=True)
        show_pop_raster = st.checkbox("Population Density (Raster)", value=False)
        show_top_ev = st.checkbox("Top 10 Optimal EV Charging Sites", value=True)
        show_top_showrooms = st.checkbox("Top 10 Optimal EV Showroom Sites", value=True)
        
        operation = st.selectbox("Spatial Operation on Selected Layers", ["None", "Union", "Intersection", "Difference"])
        process_layers = st.button("🔁 Process Layers")
    
    with col2:
        m = folium.Map(location=[28.54, 77.2], zoom_start=12)

        # --- Add Selected Vector Layers ---
        for name in selected_layers:
            if name in loaded_layers:
                gdf = loaded_layers[name]
                layer = FeatureGroup(name=name)
                
                color_mapping = {
                    "Roads": "gray", "Major Roads": "black", "Main Roads": "darkgray",
                    "Existing EV Stations": "green", "Existing Showrooms": "purple",
                    "Petrol Pumps": "red", "Hospitals": "darkred", "Universities": "blue",
                    "Markets": "orange", "Protected Areas": "red", "Historic Monuments": "brown",
                    "Waterways": "blue", "Forests": "darkgreen", "Residential Areas": "lightblue",
                    "Commercial Areas": "orange", "Industrial Areas": "darkpurple",
                    "Favorable Regions": "green", "Police Stations": "darkblue"
                }
                
                color = color_mapping.get(name, "blue")
                
                for _, row in gdf.iterrows():
                    if row.geometry.geom_type in ['Polygon', 'MultiPolygon']:
                        folium.GeoJson(
                            row.geometry, 
                            name=name,
                            style_function=lambda x, color=color: {
                                "fillColor": color,
                                "color": color,
                                "weight": 1,
                                "fillOpacity": 0.3
                            }
                        ).add_to(layer)
                    elif row.geometry.geom_type in ['Point', 'MultiPoint']:
                        folium.CircleMarker(
                            location=[row.geometry.y, row.geometry.x],
                            radius=4, color=color, fill=True, fill_opacity=0.6
                        ).add_to(layer)
                    elif row.geometry.geom_type in ['LineString', 'MultiLineString']:
                        folium.GeoJson(
                            row.geometry,
                            style_function=lambda x, color=color: {
                                "color": color,
                                "weight": 2
                            }
                        ).add_to(layer)
                m.add_child(layer)

        # --- Add South Delhi Boundary ---
        if show_south_delhi and south_delhi_boundary is not None:
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
                raster_path = layer_path / "south_delhi_pop.tif"
                if os.path.exists(raster_path):
                    with rasterio.open(raster_path) as src:
                        img_array = src.read(1)
                        bounds = src.bounds
                        if src.nodata is not None:
                            img_array = np.ma.masked_equal(img_array, src.nodata)
                        temp_img_path = base_dir / "temp_overlay.png"
                        norm = Normalize(vmin=np.min(img_array), vmax=np.max(img_array))
                        cmap = cm.get_cmap("RdYlGn_r")
                        colored = cmap(norm(img_array))
                        colored = (colored[:, :, :3] * 255).astype(np.uint8)
                        plt.imsave(temp_img_path, colored)
                        image_overlay = ImageOverlay(
                            name="Population Density (Raster)",
                            image=str(temp_img_path),
                            bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
                            opacity=0.6,
                            interactive=True,
                            cross_origin=False
                        )
                        image_overlay.add_to(m)
                else:
                    st.warning("Population density raster file not found")
            except Exception as e:
                st.error(f"Error loading raster: {e}")

        # --- Show Top EV Charging Sites ---
        if show_top_ev:
            top_sites = points_gdf.sort_values("score", ascending=False).head(10)
            ev_layer = FeatureGroup(name="Top 10 Optimal EV Charging Sites")
            for idx, row in top_sites.iterrows():
                if row.geometry.geom_type == 'Point':
                    folium.Marker(
                        location=[row.geometry.y, row.geometry.x],
                        icon=folium.Icon(color="green", icon="bolt", prefix='fa'),
                        popup=f"Rank: {idx+1}<br>Score: {row['score']:.2f}"
                    ).add_to(ev_layer)
                else:
                    st.warning(f"Skipping non-Point geometry for EV site at index {idx}: {row.geometry.geom_type}")
            m.add_child(ev_layer)

        # --- Show Top EV Showroom Sites ---
        if show_top_showrooms:
            top_showrooms = showroom_gdf.sort_values("showroom_score", ascending=False).head(10)
            showroom_layer = FeatureGroup(name="Top 10 Optimal EV Showroom Sites")
            for idx, row in top_showrooms.iterrows():
                if row.geometry.geom_type == 'Point':
                    folium.Marker(
                        location=[row.geometry.y, row.geometry.x],
                        icon=folium.Icon(color="purple", icon="car", prefix='fa'),
                        popup=f"Rank: {idx+1}<br>Score: {row['showroom_score']:.2f}"
                    ).add_to(showroom_layer)
                else:
                    st.warning(f"Skipping non-Point geometry for showroom at index {idx}: {row.geometry.geom_type}")
            m.add_child(showroom_layer)

        # --- Spatial Operations ---
        if process_layers and len(selected_layers) >= 2 and operation != "None":
            try:
                gdfs = [loaded_layers[name].to_crs(epsg=32643) for name in selected_layers]
                result = gdfs[0]
                for gdf in gdfs[1:]:
                    if operation == "Union":
                        result = result.overlay(gdf, how='union')
                    elif operation == "Intersection":
                        result = result.overlay(gdf, how='intersection')
                    elif operation == "Difference":
                        result = result.overlay(gdf, how='difference')
                
                if not result.empty:
                    result = result.to_crs(epsg=4326)
                    operation_layer = FeatureGroup(name=f"{operation} of Selected Layers")
                    for _, row in result.iterrows():
                        if row.geometry.geom_type in ['Polygon', 'MultiPolygon']:
                            folium.GeoJson(
                                row.geometry, 
                                style_function=lambda x: {
                                    "fillColor": "orange", 
                                    "color": "orange", 
                                    "weight": 1,
                                    "fillOpacity": 0.5
                                }
                            ).add_to(operation_layer)
                    m.add_child(operation_layer)
                else:
                    st.warning(f"The {operation} operation resulted in an empty geometry")
            except Exception as e:
                st.error(f"Error performing spatial operation: {e}")

        folium.LayerControl(collapsed=False).add_to(m)
        st_data = st_folium(m, width=800, height=650)

with tab2:
    st.subheader("📍 Optimal EV Charging Station Locations")
    
    col1, col2 = st.columns([4, 6])
    
    with col1:
        st.write("""
        This analysis uses a weighted overlay approach to identify optimal locations 
        for new EV charging stations. The scoring considers:
        
        **Positive Factors:**
        - 🅿️ Proximity to parking areas (+10)
        - 🏥 Hospitals (+8)
        - 🎓 Universities (+8)
        - 🛒 Markets (+7)
        - 🏘️ Residential areas (+10)
        - 🏢 Commercial areas (+9)
        - 🚓 Police stations (+6)
        - 🛣️ Major roads (+8)
        - 🚌 Bus stops (+7)
        - ✅ Favorable regions (+15)
        
        **Negative Factors:**
        - 🌳 Forest areas (-10)
        - 🏛️ Historic monuments (-10)
        - ⚠️ Protected areas (-10)
        - 💧 Water bodies (-8)
        - ⚡ Existing EV stations (-20)
        - ⛽ Petrol pumps (-5)
        - 🏪 Optimal showrooms (-15)
        """)
        
        lat = st.number_input("Latitude", value=28.54, format="%.6f")
        lon = st.number_input("Longitude", value=77.21, format="%.6f")
        
        fig, ax = plt.subplots(figsize=(8, 4))
        points_gdf["score"].hist(bins=20, ax=ax)
        ax.set_title("Distribution of EV Charging Station Scores")
        ax.set_xlabel("Score")
        ax.set_ylabel("Frequency")
        st.pyplot(fig)
        
        st.subheader("Top Recommended EV Charging Sites")
        top_df = points_gdf.sort_values("score", ascending=False).head(10)
        top_df["rank"] = range(1, len(top_df) + 1)
        top_df = top_df[["rank", "score"]].reset_index(drop=True)
        st.dataframe(top_df)
        
    with col2:
        map2 = folium.Map(location=[lat, lon], zoom_start=12)
        
        folium.Marker([lat, lon], tooltip="Your Location", icon=folium.Icon(color='blue')).add_to(map2)
        
        all_scores = points_gdf["score"].values
        min_score = min(all_scores)
        max_score = max(all_scores)
        score_range = max_score - min_score
        
        visible_sites = points_gdf[points_gdf["score"] > 0].sort_values("score", ascending=False)
        
        for idx, row in visible_sites.head(50).iterrows():
            if row.geometry.geom_type == 'Point':
                norm_score = (row["score"] - min_score) / score_range
                icon_color = "green"
                if idx < 10:
                    icon = folium.Icon(color=icon_color, icon="bolt", prefix='fa')
                else:
                    icon = folium.Icon(color=icon_color, icon="dot-circle", prefix='fa')
                
                folium.Marker(
                    location=[row.geometry.y, row.geometry.x],
                    icon=icon,
                    popup=f"Score: {row['score']:.2f}"
                ).add_to(map2)
            else:
                st.warning(f"Skipping non-Point geometry in visible sites at index {idx}: {row.geometry.geom_type}")
        
        if not visible_sites.empty:
            visible_sites["dist"] = visible_sites.geometry.apply(
                lambda geom: geodesic((lat, lon), (geom.y, geom.x)).meters if geom.geom_type == 'Point' else float('inf')
            )
            nearest = visible_sites.sort_values("dist").iloc[0]
            if nearest.geometry.geom_type == 'Point':
                nearest_coords = [nearest.geometry.y, nearest.geometry.x]
                folium.Marker(
                    location=nearest_coords,
                    icon=folium.Icon(color="red", icon="bolt", prefix='fa'),
                    popup=f"Nearest Optimal EV Station<br>Distance: {nearest['dist']:.2f} meters<br>Score: {nearest['score']:.2f}"
                ).add_to(map2)
                folium.PolyLine(
                    [[lat, lon], nearest_coords], 
                    color="blue", 
                    weight=2.5, 
                    opacity=0.7,
                    tooltip=f"Distance: {nearest['dist']:.2f} meters"
                ).add_to(map2)
                st.info(f"🚶 **Distance to nearest optimal EV Charging Station**: {nearest['dist']:.2f} meters")
        
        if "Existing EV Stations" in loaded_layers:
            existing_ev = loaded_layers["Existing EV Stations"]
            for _, row in existing_ev.iterrows():
                if row.geometry.geom_type == 'Point':
                    folium.CircleMarker(
                        location=[row.geometry.y, row.geometry.x],
                        radius=6,
                        color="black",
                        fill=True,
                        fill_color="yellow",
                        fill_opacity=0.7,
                        tooltip="Existing EV Station"
                    ).add_to(map2)
        
        st_folium(map2, width=800, height=600)

with tab3:
    st.subheader("🏪 Optimal EV Showroom Locations")
    
    col1, col2 = st.columns([4, 6])
    
    with col1:
        st.write("""
        This analysis identifies optimal locations for new EV showrooms based on:
        
        **Positive Factors:**
        - 🅿️ Proximity to parking areas (+10)
        - 🏢 Commercial areas (+12)
        - 🏘️ Residential areas (+7)
        - 🛒 Markets (+8)
        - 🛣️ Major roads (+12)
        - 🛣️ Main roads (+9)
        - 🚓 Police stations (+4)
        - 🚌 Bus stops (+6)
        - ✅ Favorable regions (+10)
        - ⚡ Near existing EV stations (+5)
        
        **Negative Factors:**
        - ⚠️ Protected areas (-10)
        - 🌳 Forest areas (-10)
        - 💧 Water bodies (-8)
        - 🏛️ Historic monuments (-5)
        - 📍 Optimal charging stations (-15)
        """)
        
        lat = st.number_input("Latitude (Showroom)", value=28.54, format="%.6f", key="lat3")
        lon = st.number_input("Longitude (Showroom)", value=77.21, format="%.6f", key="lon3")
        
        fig, ax = plt.subplots(figsize=(8, 4))
        showroom_gdf["showroom_score"].hist(bins=20, ax=ax)
        ax.set_title("Distribution of EV Showroom Scores")
        ax.set_xlabel("Score")
        ax.set_ylabel("Frequency")
        st.pyplot(fig)
        
        st.subheader("Top Recommended EV Showroom Sites")
        top_df = showroom_gdf.sort_values("showroom_score", ascending=False).head(10)
        top_df["rank"] = range(1, len(top_df) + 1)
        top_df = top_df[["rank", "showroom_score"]].reset_index(drop=True)
        st.dataframe(top_df)
        
    with col2:
        map3 = folium.Map(location=[lat, lon], zoom_start=12)
        
        folium.Marker([lat, lon], tooltip="Your Location", icon=folium.Icon(color='blue')).add_to(map3)
        
        all_scores = showroom_gdf["showroom_score"].values
        min_score = min(all_scores)
        max_score = max(all_scores)
        score_range = max_score - min_score
        
        visible_sites = showroom_gdf[showroom_gdf["showroom_score"] > 0].sort_values("showroom_score", ascending=False)
        
        for idx, row in visible_sites.head(50).iterrows():
            if row.geometry.geom_type == 'Point':
                norm_score = (row["showroom_score"] - min_score) / score_range
                icon_color = "purple"
                if idx < 10:
                    icon = folium.Icon(color=icon_color, icon="car", prefix='fa')
                else:
                    icon = folium.Icon(color=icon_color, icon="dot-circle", prefix='fa')
                
                folium.Marker(
                    location=[row.geometry.y, row.geometry.x],
                    icon=icon,
                    popup=f"Score: {row['showroom_score']:.2f}"
                ).add_to(map3)
            else:
                st.warning(f"Skipping non-Point geometry in showroom visible sites at index {idx}: {row.geometry.geom_type}")
        
        if not visible_sites.empty:
            visible_sites["dist"] = visible_sites.geometry.apply(
                lambda geom: geodesic((lat, lon), (geom.y, geom.x)).meters if geom.geom_type == 'Point' else float('inf')
            )
            nearest = visible_sites.sort_values("dist").iloc[0]
            if nearest.geometry.geom_type == 'Point':
                nearest_coords = [nearest.geometry.y, nearest.geometry.x]
                folium.Marker(
                    location=nearest_coords,
                    icon=folium.Icon(color="red", icon="car", prefix='fa'),
                    popup=f"Nearest Optimal EV Showroom<br>Distance: {nearest['dist']:.2f} meters<br>Score: {nearest['showroom_score']:.2f}"
                ).add_to(map3)
                folium.PolyLine(
                    [[lat, lon], nearest_coords], 
                    color="purple", 
                    weight=2.5, 
                    opacity=0.7,
                    tooltip=f"Distance: {nearest['dist']:.2f} meters"
                ).add_to(map3)
                st.info(f"🚶 **Distance to nearest optimal EV Showroom**: {nearest['dist']:.2f} meters")
        
        if "Existing Showrooms" in loaded_layers:
            existing_showrooms = loaded_layers["Existing Showrooms"]
            for _, row in existing_showrooms.iterrows():
                if row.geometry.geom_type == 'Point':
                    folium.CircleMarker(
                        location=[row.geometry.y, row.geometry.x],
                        radius=6,
                        color="black",
                        fill=True,
                        fill_color="orange",
                        fill_opacity=0.7,
                        tooltip="Existing Showroom"
                    ).add_to(map3)
        
        st_folium(map3, width=800, height=600)

with tab4:
    st.subheader("📊 Spatial Analysis")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("#### EV Charging Infrastructure Analysis")
        
        total_points = len(points_gdf)
        optimal_points = len(points_gdf[points_gdf["score"] > 20])
        optimal_percentage = (optimal_points / total_points) * 100 if total_points > 0 else 0
        
        st.metric("Total Grid Points", f"{total_points:,}")
        st.metric("Optimal EV Charging Locations", f"{optimal_points:,} ({optimal_percentage:.1f}%)")
        
        st.write("#### Average Scores by Feature Proximity")
        
        fig, ax = plt.subplots(figsize=(10, 6))
        feature_scores = {
            "Near Commercial": points_gdf[points_gdf["score"] > 0]["score"].mean(),
            "Near Residential": points_gdf[points_gdf["score"] > 5]["score"].mean(),
            "Near Markets": points_gdf[points_gdf["score"] > 10]["score"].mean(),
            "Near Roads": points_gdf[points_gdf["score"] > 15]["score"].mean(),
            "All Points": points_gdf["score"].mean()
        }
        
        features = list(feature_scores.keys())
        scores = list(feature_scores.values())
        
        ax.bar(features, scores, color='skyblue')
        ax.set_title("Average EV Station Scores by Context")
        ax.set_ylabel("Average Score")
        ax.tick_params(axis='x', rotation=45)
        st.pyplot(fig)
        
    with col2:
        st.write("#### EV Showroom Location Analysis")
        
        total_points = len(showroom_gdf)
        optimal_points = len(showroom_gdf[showroom_gdf["showroom_score"] > 20])
        optimal_percentage = (optimal_points / total_points) * 100 if total_points > 0 else 0
        
        st.metric("Total Grid Points", f"{total_points:,}")
        st.metric("Optimal EV Showroom Locations", f"{optimal_points:,} ({optimal_percentage:.1f}%)")
        
        st.write("#### Comparison: Charging vs Showroom Locations")
        
        optimal_charging = points_gdf[points_gdf["score"] > 20]
        optimal_showroom = showroom_gdf[showroom_gdf["showroom_score"] > 20]
        
        if not optimal_charging.empty and not optimal_showroom.empty:
            optimal_charging_gdf = gpd.GeoDataFrame(optimal_charging, geometry="geometry", crs="EPSG:4326")
            optimal_showroom_gdf = gpd.GeoDataFrame(optimal_showroom, geometry="geometry", crs="EPSG:4326")
            
            charging_points = set([p.wkt for p in optimal_charging_gdf.geometry if p.geom_type == 'Point'])
            showroom_points = set([p.wkt for p in optimal_showroom_gdf.geometry if p.geom_type == 'Point'])
            both_optimal = charging_points.intersection(showroom_points)
            both_count = len(both_optimal)
            
            labels = ['Charging Only', 'Showroom Only', 'Both']
            sizes = [len(charging_points) - both_count, len(showroom_points) - both_count, both_count]
            colors = ['lightgreen', 'lightblue', 'orange']
            
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90)
            ax.axis('equal')
            ax.set_title("Optimal Location Types")
            st.pyplot(fig)
            
            if both_count > 0:
                st.success(f"**{both_count}** locations are optimal for both EV charging stations and showrooms!")
                st.write("These strategic locations could be developed as comprehensive EV hubs.")
            else:
                st.info("No locations are optimal for both charging stations and showrooms.")
                st.write("This suggests different spatial requirements for these two infrastructure types.")
        else:
            st.warning("Insufficient data to compare optimal locations.")