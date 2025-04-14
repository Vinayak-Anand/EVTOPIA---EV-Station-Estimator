import geopandas as gpd
import pandas as pd
import os
from shapely.geometry import Point
from tqdm import tqdm

# Set paths
layer_path = "geota_final_layers"
output_path = "ev_station_estimator/output"
os.makedirs(output_path, exist_ok=True)

# Load necessary layers
def load_layer(name):
    return gpd.read_file(f"{layer_path}/{name}.gpkg")

layers = {
    "petrol": load_layer("petrol_pumps"),
    "buffered_petrol": load_layer("buffered_petrol_pumps"),
    "ev": load_layer("ev_stations"),
    "parking": load_layer("parking"),
    "hospital": load_layer("hospitals"),
    "university": load_layer("university"),
    "market": load_layer("market_places"),
    "amenities": load_layer("landuse_amenity"),
    "residential": load_layer("landuse_residential"),
    "commercial": load_layer("landuse_commercial"),
    "industrial": load_layer("landuse_industrial"),
    "forest": load_layer("landuse_forest"),
    "historic": load_layer("historic_monuments"),
    "protected": load_layer("protected_areas"),
    "water": load_layer("waterways"),
    "power": load_layer("powerlines"),
    "roads": load_layer("roads"),
    "major_roads": load_layer("major_roads_only"),
    "main_roads": load_layer("main_roads"),
    "showroom": load_layer("car_showrooms"),
    "police": load_layer("police_stations"),
    "buildings": load_layer("buildings"),
    "bus_stops": load_layer("bus_stops"),
    "favourable_regions": load_layer("favourable_regions"),
    "union_regions": load_layer("union_regions"),
    "useable_area": load_layer("useable_area"),
    "places": load_layer("places"),
}

# Load South Delhi boundary
south_delhi_boundary = gpd.read_file(f"{layer_path}/south_delhi_layer.shp")
south_delhi_boundary = south_delhi_boundary.to_crs(epsg=32643)

# Generate grid of points within South Delhi
minx, miny, maxx, maxy = south_delhi_boundary.total_bounds
x_points = int((maxx - minx) / 250)
y_points = int((maxy - miny) / 250)

points = [
    Point(minx + i * 250, miny + j * 250)
    for i in range(x_points)
    for j in range(y_points)
]

grid_gdf = gpd.GeoDataFrame(geometry=points, crs=south_delhi_boundary.crs)
grid_gdf = grid_gdf[grid_gdf.within(south_delhi_boundary.geometry.unary_union)]

# Reproject all layers
for key in layers:
    layers[key] = layers[key].to_crs(epsg=32643)

# Pre-buffer layers once for performance
buffered_layers = {
    "parking": layers["parking"].buffer(100),
    "hospital": layers["hospital"].buffer(200),
    "university": layers["university"].buffer(200),
    "market": layers["market"].buffer(150),
    "residential": layers["residential"].buffer(150),
    "commercial": layers["commercial"].buffer(150),
    "industrial": layers["industrial"].buffer(100),
    "forest": layers["forest"].buffer(100),
    "historic": layers["historic"].buffer(100),
    "protected": layers["protected"].buffer(100),
    "water": layers["water"].buffer(100),
    "ev": layers["ev"].buffer(150),
    "petrol": layers["petrol"].buffer(150),
    "roads": layers["roads"].buffer(100),
    "major_roads": layers["major_roads"].buffer(150),
    "main_roads": layers["main_roads"].buffer(120),
    "police": layers["police"].buffer(200),
    "buildings": layers["buildings"].buffer(50),
    "bus_stops": layers["bus_stops"].buffer(150),
    "favourable_regions": layers["favourable_regions"].buffer(0),
    "useable_area": layers["useable_area"].buffer(0),
    "places": layers["places"].buffer(150)
}

buffered_gdfs = {
    key: gpd.GeoDataFrame(geometry=geom, crs=32643) for key, geom in buffered_layers.items()
}

# Compute EV Charging Station score
def compute_score(point, showroom_buffer=None):
    score = 0
    for key, value in [
        ("parking", 10), ("hospital", 8), ("university", 8),
        ("market", 7), ("residential", 10), ("commercial", 9),
        ("industrial", 5), ("forest", -10), ("historic", -10),
        ("protected", -10), ("water", -8), ("ev", -10),
        ("petrol", -5), ("police", 6), ("major_roads", 8),
        ("main_roads", 6), ("bus_stops", 7), ("favourable_regions", 15),
        ("useable_area", 5)
    ]:
        hits = buffered_gdfs[key].sindex.query(point, predicate="intersects")
        if len(hits) > 0:
            score += value
    # Penalize proximity to optimal showrooms
    if showroom_buffer is not None:
        hits = showroom_buffer.sindex.query(point, predicate="intersects")
        if len(hits) > 0:
            score -= 15  # Negative score to discourage overlap
    return score

# Compute EV Showroom score
def compute_showroom_score(point, charging_buffer=None):
    score = 0
    for key, value in [
        ("parking", 10), ("commercial", 12), ("residential", 7),
        ("market", 8), ("roads", 7), ("major_roads", 12),
        ("main_roads", 9), ("protected", -10), ("forest", -10),
        ("water", -8), ("historic", -5), ("industrial", 6),
        ("police", 4), ("bus_stops", 6), ("favourable_regions", 10),
        ("ev", 10)
    ]:
        hits = buffered_gdfs[key].sindex.query(point, predicate="intersects")
        if len(hits) > 0:
            score += value
    # Penalize proximity to optimal charging stations
    if charging_buffer is not None:
        hits = charging_buffer.sindex.query(point, predicate="intersects")
        if len(hits) > 0:
            score -= 15  # Negative score to discourage overlap
    return score

# Compute initial scores
tqdm.pandas(desc="Initial Scoring")
grid_gdf["score"] = grid_gdf["geometry"].progress_apply(compute_score)
grid_gdf["showroom_score"] = grid_gdf["geometry"].progress_apply(compute_showroom_score)

# Identify high-scoring sites (threshold can be adjusted)
optimal_charging = grid_gdf[grid_gdf["score"] >= 20].copy()
optimal_showrooms = grid_gdf[grid_gdf["showroom_score"] >= 20].copy()

# Create buffers for optimal sites to penalize proximity
charging_buffer = optimal_charging["geometry"].buffer(150)  # 150-meter buffer
showroom_buffer = optimal_showrooms["geometry"].buffer(150)  # 150-meter buffer
charging_buffer_gdf = gpd.GeoDataFrame(geometry=charging_buffer, crs=32643)
showroom_buffer_gdf = gpd.GeoDataFrame(geometry=showroom_buffer, crs=32643)

# Recompute scores with proximity penalties
tqdm.pandas(desc="Adjusted Scoring")
grid_gdf["score"] = grid_gdf["geometry"].progress_apply(lambda x: compute_score(x, showroom_buffer_gdf))
grid_gdf["showroom_score"] = grid_gdf["geometry"].progress_apply(lambda x: compute_showroom_score(x, charging_buffer_gdf))

# Validate geometry types before saving
if not all(grid_gdf.geometry.geom_type == 'Point'):
    non_points = grid_gdf[grid_gdf.geometry.geom_type != 'Point']
    print(f"Warning: Found {len(non_points)} non-Point geometries: {non_points.geometry.geom_type.unique()}")
    grid_gdf = grid_gdf[grid_gdf.geometry.geom_type == 'Point']  # Keep only Points

# Save outputs
grid_gdf.to_csv(f"{output_path}/scored_points.csv", index=False)
grid_gdf[grid_gdf["score"] >= 20].to_file(f"{output_path}/optimal_sites.geojson", driver="GeoJSON")
grid_gdf.to_csv(f"{output_path}/scored_showroom_points.csv", index=False)
grid_gdf[grid_gdf["showroom_score"] >= 20].to_file(f"{output_path}/optimal_showroom_sites.geojson", driver="GeoJSON")

print(f"Total points evaluated: {len(grid_gdf)}")
print(f"Optimal EV charging sites found: {len(grid_gdf[grid_gdf['score'] >= 20])}")
print(f"Optimal EV showroom sites found: {len(grid_gdf[grid_gdf['showroom_score'] >= 20])}")