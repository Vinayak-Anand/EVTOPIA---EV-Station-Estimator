import geopandas as gpd
import pandas as pd
import os
from shapely.geometry import Point
from tqdm import tqdm

# Set paths
# Set paths
layer_path = "../geota_final_layers"
output_path = "../ev_station_estimator/output"
os.makedirs(output_path, exist_ok=True)


# Load necessary layers
def load_layer(name):
    return gpd.read_file(f"{layer_path}/{name}.gpkg")

layers = {
    "petrol": load_layer("petrol_pumps"),
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
    "roads": load_layer("roads")
}

south_delhi_boundary = gpd.read_file(f"{layer_path}/south_delhi_layer.shp")
south_delhi_boundary = south_delhi_boundary.to_crs(epsg=32643)  # UTM Zone 43N

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
grid_gdf = grid_gdf[grid_gdf.within(south_delhi_boundary.unary_union)]

# Reproject all layers
for key in layers:
    layers[key] = layers[key].to_crs(epsg=32643)

def compute_score(point):
    score = 0
    # Positive weights
    if layers["parking"].buffer(100).contains(point).any():
        score += 10
    if layers["hospital"].buffer(200).contains(point).any():
        score += 8
    if layers["university"].buffer(200).contains(point).any():
        score += 8
    if layers["market"].buffer(150).contains(point).any():
        score += 7
    if layers["residential"].buffer(150).contains(point).any():
        score += 10
    if layers["commercial"].buffer(150).contains(point).any():
        score += 7

    # Negative weights
    if layers["forest"].buffer(100).contains(point).any():
        score -= 10
    if layers["historic"].buffer(100).contains(point).any():
        score -= 10
    if layers["protected"].buffer(100).contains(point).any():
        score -= 10
    if layers["water"].buffer(100).contains(point).any():
        score -= 8

    # Already present infra
    if layers["ev"].buffer(150).contains(point).any():
        score -= 20
    if layers["petrol"].buffer(150).contains(point).any():
        score -= 5

    return score

# Compute scores
tqdm.pandas(desc="Scoring")
grid_gdf["score"] = grid_gdf["geometry"].progress_apply(compute_score)

# Save outputs
grid_gdf.to_csv(f"{output_path}/scored_points.csv", index=False)
grid_gdf[grid_gdf["score"] >= 20].to_file(f"{output_path}/optimal_sites.geojson", driver="GeoJSON")
