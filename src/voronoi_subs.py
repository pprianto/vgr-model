import os
import numpy as np
import pandas as pd
import geopandas as gpd
from geovoronoi import voronoi_regions_from_coords, points_to_coords
import folium as fl
from shapely.ops import unary_union
from basics import (prerequisites,
                    savefile,
                    reg_bound)
from IPython.display import display


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram/20678647#20678647

    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def nearest_subs(
        point, 
        subs_gdf
):

    """
    Get index for substation node id according 
    to the voronoi region.

    Maybe not needed after all
    """

    subs_gdf = subs_gdf.to_crs("EPSG:3006")

    distances = subs_gdf.distance(point)
    nearest_idx = distances.idxmin()

    subs_gdf = subs_gdf.to_crs("EPSG:4326")

    return subs_gdf.loc[nearest_idx, "node_id"]


def vormap(
        region_gdf,
        vor_gdf,
        subs_gdf,
        pp_gdf,
        lines_gdf
):

    # Initialize a Folium map centered on the region
    m = fl.Map(
        location=[region_gdf.geometry.centroid.y.mean(), region_gdf.geometry.centroid.x.mean()], 
        zoom_start=7
    )

    # Add the region boundary to the map
    fl.GeoJson(
        region_gdf, 
        name="Region Boundary",
        style_function=lambda feature: {
            "fillColor": "orange",
            "color": "black",
            "weight": 1,
            "dashArray": "5, 5"
            }
    ).add_to(m)

    # Add the Voronoi cells to the map
    for idx, row in vor_gdf.iterrows():
        fl.GeoJson(
            row["geometry"],
            name=f"Voronoi Cell {idx}",
            style_function=lambda feature: {
                "color": "blue", 
                "fillColor": "yellow",
                "weight": 1,
                "fillOpacity": 0.8,
                "dashArray": "5, 5"
                }
        ).add_to(m)

    # Add substations to the map
    fl.GeoJson(
        subs_gdf,
        name="Substation",
        marker=fl.Circle(
            radius=10, 
            fill_color="black", 
            fill_opacity=0.4, 
            color="black", 
            weight=5
        ),
        tooltip=fl.GeoJsonTooltip(
            fields=[
                "node_id", 
                "voltage", 
                "substation"
        ]),
    ).add_to(m)

    # Add power plants to the map
    fl.GeoJson(
        pp_gdf,
        name="Power Plant",
        marker=fl.Circle(
            radius=6, 
            fill_color="green", 
            fill_opacity=0.5, 
            color="green", 
            weight=2
        ),
        tooltip=fl.GeoJsonTooltip(
            fields=[
                "el_MW", 
                "heat_MW", 
                "tech"
        ]),
    ).add_to(m)

    # Add power lines to the map
    fl.GeoJson(
        lines_gdf,
        name="Power Lines",
        style_function=lambda feature: {
            "color": "red",
            "weight": 1
        }
    ).add_to(m)

    # Add layer control to the map
    fl.LayerControl().add_to(m)

    # Save the map to an HTML file
    # m.save("interactive_voronoi_map.html")
    display(m)


def subs_voronoi(
        subs_gdf,
        pp_gdf,
        kommun_gdf=None,
        lan_gdf=None
):

    """
    Process the voronoi diagram
    based on the substations and bounded a certain region
    to assign power plants to substation
    corresponds to the voronoi diagram

    Inputs
    ----------
    subs_gdf : dataframe with geometry
        substation dataframe
    pp_gdf : dataframe with geometry
        power plant dataframe    
    kommun_gdf : dataframe with geometry, optional
        geodataframe for municipalities.
    lan_gdf : dataframe with geometry, optional
        geodataframe for counties.
        

    Returns
    -------
    pp_gdf : updated pp_gdf
        Includes the voronoi regions and assigned substations.
    voronoi_gdf : voronoi dataframe
        If only interesting.    

    """

    # substations as voronoi centres
    vor_nodes = points_to_coords(subs_gdf.geometry)

    # use county boundaries are limits for voronoi regions
    kommun_border = kommun_gdf.unary_union
    lan_border = lan_gdf.unary_union

    # voronoi polygons
    # kommun border is chosen since it gives the most voronoi regions
    polygons, points, unassigned_pts = voronoi_regions_from_coords(
                                            vor_nodes, 
                                            kommun_border,
                                            return_unassigned_points=True,
                                            per_geom=True
                                        )

    # define the voronoi as gdf
    # voronoi_gdf = gpd.GeoDataFrame(geometry=vor_polygons, crs="EPSG:4326")    
    voronoi_gdf = gpd.GeoDataFrame.from_dict(
                                    points, 
                                    orient="index",
                                    columns=["idx"],
                                )
    # voronoi_gdf = gpd.GeoDataFrame({"geometry": polygons}, crs="EPSG:4326")
    voronoi_gdf["geometry"] = polygons
    voronoi_gdf.set_geometry("geometry", inplace=True, crs="EPSG:4326")
    voronoi_gdf["vor_id"] = range(1, len(voronoi_gdf) + 1) 
    voronoi_gdf = voronoi_gdf.to_crs("EPSG:3006")
    voronoi_gdf["area_m2"] = voronoi_gdf.area
    # voronoi_gdf = voronoi_gdf.to_crs("EPSG:4326")

    idx_to_node = subs_gdf["node_id"].to_dict()

    voronoi_gdf["node_id"] = voronoi_gdf["idx"].map(idx_to_node)

    # assign power plant to corresponding voronoi region
    
    # pp_gdf = gpd.sjoin(pp_gdf, voronoi_gdf, how="left", predicate="within")
    pp_gdf = pp_gdf.to_crs("EPSG:3006")
    pp_gdf = gpd.sjoin_nearest(pp_gdf, voronoi_gdf, how="left")
    # pp_gdf["subs_node"] = pp_gdf["geometry"].apply(lambda x: nearest_subs(x, subs_gdf))
    pp_gdf.drop(["index_right", "idx"], axis=1, inplace=True)

    # revert back the CRS
    voronoi_gdf = voronoi_gdf.to_crs("EPSG:4326")
    pp_gdf = pp_gdf.to_crs("EPSG:4326")

    # update subs_gdf with voronoi informations
    subs_gdf = subs_gdf.merge(
                    voronoi_gdf[["vor_id", "area_m2", "node_id"]], 
                    on="node_id", 
                    suffixes=("_subs", "_vor")
            )

    # rearrange columns
    subs_gdf = subs_gdf[[
                    "bus_id",
                    "station_id",
                    "node_id",
                    "voltage",
                    "symbol",
                    "substation",
                    "import_trans",
                    "KnKod",
                    "KnNamn",
                    "municipality",
                    "lon",
                    "lat",
                    "vor_id",
                    "area_m2",
                    "geometry"
                ]]    
        
    pp_gdf = pp_gdf[[
                    "pp_id",
                    "el_MW",
                    "heat_MW",
                    "tech",
                    "method",
                    "operator",
                    "name",
                    "KnKod",
                    "KnNamn",
                    "municipality",
                    "lon",
                    "lat",
                    "node_id",
                    "vor_id",
                    "area_m2",
                    "geometry"
                ]]
    
    return subs_gdf, pp_gdf, voronoi_gdf


if __name__ == "__main__":

    # directory for input files
    input_folder = prerequisites("inputs")
    model_input = prerequisites("modelinput")

    # read files
    kommun_fn = os.path.join(input_folder, "Kommun_Sweref99TM.shp")
    kommun_gdf = gpd.read_file(kommun_fn)
    kommun_gdf = reg_bound(kommun_gdf)

    lan_fn = os.path.join(input_folder, "Lan_Sweref99TM_region.shp")
    lan_gdf = gpd.read_file(lan_fn)
    lan_gdf = reg_bound(lan_gdf)

    # these files are already in EPSG:4326 coordinates
    subs_fn = os.path.join(input_folder, "subs_clean.geojson")
    subs_gdf = gpd.read_file(subs_fn)

    pp_fn = os.path.join(input_folder, "pp_clean.geojson")
    pp_gdf = gpd.read_file(pp_fn)

    lines_fn = os.path.join(input_folder, "lines_clean.geojson")
    lines_gdf = gpd.read_file(lines_fn)

    # process voronoi diagrams
    subs_gdf, pp_gdf, voronoi_gdf = subs_voronoi(
        subs_gdf,
        pp_gdf,
        kommun_gdf,
        lan_gdf
    )

    # save file
    savefile(subs_gdf, "subs_final", model_input)
    savefile(lines_gdf, "lines_final", model_input)
    savefile(pp_gdf, "pp_final", model_input)

    # plot the map
    vormap(
        lan_gdf,
        voronoi_gdf,
        subs_gdf,
        pp_gdf,
        lines_gdf
    )


