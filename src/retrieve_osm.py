import os
import osmnx as ox
import geopandas as gpd
import numpy as np
import pandas as pd
import folium as fl
from IPython.display import display
from basics import (prerequisites,
                    savefile,
                    reg_bound,
                    getinputfiles)

pd.set_option("mode.copy_on_write", True)
pd.options.mode.copy_on_write = True

"""
The following list of functions are adapted from
PyPSA-Earth
https://github.com/pypsa-meets-earth/pypsa-earth
https://arxiv.org/abs/2209.04663
"""

def prep_subs_df(subs_df): # Modified
    """
    Prepare substations dataframe

    Parameters
    ----------
    subs_df : dataframe or geodataframe
        Raw substations geodataframe as downloaded from OpenStreetMap with OSMNX
    """
    
    # Convert the polygon geometry into centroid through projection
    subs_df = subs_df.reset_index(drop=False)
    # subs_df = subs_df.to_crs("EPSG:3035") # https://epsg.io/3035 used in EU
    subs_df = subs_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    subs_df['geometry'] = subs_df['geometry'].centroid
    subs_df = subs_df.to_crs("EPSG:4326") # https://epsg.io/4326 World Geodetic System 1984

    # Modify the naming of the DataFrame columns
    subs_df = subs_df.rename(
        columns={
            "osmid": "bus_id",
            "voltage": "voltage",
            # "dc", # not needed in vgr case
            "power": "symbol",
            # "under_construction", not needed in vgr case
            "substation": "substation",
            "Country": "country", # should be replaced by Kommuns
            "Area": "area",
            "geometry": "geometry",
        }
    )

    # Add longitude/latitude and connection to transmission
    subs_df["lon"] = subs_df["geometry"].x
    subs_df["lat"] = subs_df["geometry"].y
    subs_df["import_trans"] = subs_df["voltage"].str.contains(
                                        "400000|420000|220000", na = False)

    # # Initialize columns to default value
    # subs_df["dc"] = False
    # subs_df["under_construction"] = False

    # Rearrange columns
    clist = [
        "bus_id",
        # "station_id",
        "voltage",
        # "dc",
        "symbol",
        # "under_construction",
        "substation",
        # "area",
        "lon",
        "lat",
        "geometry",
        # "country",
        "import_trans",
    ]

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in subs_df:
            subs_df[c] = np.nan

    # ASSUMPTION that an empty cable info has at least 3 cables
    if subs_df['substation'].isna().any:
        # subs_df['substation'].fillna('transmission', inplace=True)
        subs_df.fillna({'substation' : 'transmission'}, inplace=True)

    subs_df.drop(
        subs_df.columns[~subs_df.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    return subs_df


def prep_lines_df(lines_df): # Modified
    """
    Prepares lines and cables dataframe.

    Parameters
    ----------
    lines_df : dataframe
        Raw lines or cables dataframe as downloaded from OpenStreetMap
    """

    # Convert the polygon geometry into centroid through projection
    lines_df = lines_df.reset_index(drop=False)
    # lines_df = lines_df.to_crs("EPSG:3857") # https://epsg.io/3857 for distance in m
    lines_df = lines_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    lines_df['length_km'] = lines_df.length / 1000
    lines_df = lines_df.to_crs("EPSG:4326")

    # Modification - create final dataframe layout
    lines_df = lines_df.rename(
        columns={
            "osmid": "line_id",
            "voltage": "voltage",
            "circuits": "circuits",
            "cables": "cables",
            "wires": "wires",
            "frequency": "frequency",
            "location": "underground",
            "power": "types",
            "geometry": "geometry",
            "Country": "country",
            "length_km": "length_km",
        }
    )

    # Rearrange columns
    clist = [
        "line_id",
        "bus0",
        "bus1",
        "voltage",
        "circuits",
        "length_km",
        # "under_construction",
        "types",
        # "frequency",
        # "dc",
        "cables",
        "wires",
        # "underground",
        "geometry",
        # "country",
    ]

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in lines_df:
            lines_df[c] = np.nan

    lines_df.drop(
        lines_df.columns[~lines_df.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    return lines_df


def prep_gens_df(gens_df): # Modified
    """
    Prepare generators dataframe.
    """
    # reset index and convert the polygon shapes into centroids
    gens_df = gens_df.reset_index(drop=False)
    # gens_df = gens_df.to_crs("EPSG:3035") # https://epsg.io/3035 used in EU
    gens_df = gens_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    gens_df['geometry'] = gens_df['geometry'].centroid
    gens_df = gens_df.to_crs("EPSG:4326")

    check_fields = ["generator:output:electricity", 
                    "generator:output:hot_water",
                    # "generator:output:biogas" # not included
                   ]

    for field_to_add in check_fields:
        if field_to_add not in gens_df.columns.tolist():
            gens_df[field_to_add] = ""

    gens_df = gens_df.rename(
        columns={
            "osmid": "pp_id",
            "generator:output:electricity": "el_output_MW",
            "generator:output:hot_water": "heat_output_MW",
            "plant:source": "gen_source",
            "geometry": "geometry",
            "operator": "operator",
            "name": "name",
        }
    )

    # Rearrange columns
    clist = [
        "pp_id",
        "el_output_MW",
        "heat_output_MW",
        "gen_source",
        "geometry",
        "operator",
        "name",
    ]

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in gens_df:
            gens_df[c] = np.nan

    gens_df.drop(
        gens_df.columns[~gens_df.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    # drop rows that are NaN and do not have numeric information (e.g: yes, el, heat, etc.)
    gens_df.dropna(subset=["el_output_MW", "heat_output_MW"], how="all", inplace=True)

    # convert electricity column from string to float value
    # TODO: this filtering can be improved
    cols = ["el_output_MW", "heat_output_MW"] # after renaming

    for col in cols:
        gens_df = gens_df[
            (gens_df[col].astype(str).str.contains("MW")) | gens_df[col].isna()
        ]
        
        gens_df[col] = (
            gens_df[col]
            .astype(str)
            .str.extract("(\\d+)")
            .astype(float)
        )

    return gens_df


def prep_plant_df(plant_df): # Modified
    """
    Prepare generators dataframe.
    """
    # reset index and convert the polygon shapes into centroids
    plant_df = plant_df.reset_index(drop=False)
    # plant_df = plant_df.to_crs("EPSG:3035") # https://epsg.io/3035 used in EU
    plant_df = plant_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    plant_df['geometry'] = plant_df['geometry'].centroid
    plant_df = plant_df.to_crs("EPSG:4326")

    check_fields = ["plant:output:electricity", 
                    "plant:output:hot_water",
                    # "plant:output:cold_water" # not included
                   ]

    for field_to_add in check_fields:
        if field_to_add not in plant_df.columns.tolist():
            plant_df[field_to_add] = ""

    # Modification - create final dataframe layout
    plant_df = plant_df.rename(
        columns={
            "osmid": "pp_id",
            "plant:output:electricity": "el_output_MW",
            "plant:output:hot_water": "heat_output_MW",
            # "plant:output:cold_water": "cold_output_MW",
            "plant:source": "gen_source",
            "geometry": "geometry",
            "operator": "operator",
            "name": "name",
        }
    )

    # Rearrange columns
    clist = [
        "pp_id",
        "el_output_MW",
        "heat_output_MW",
        "gen_source",
        "geometry",
        "operator",
        "name",
    ]

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in plant_df:
            plant_df[c] = np.nan

    plant_df.drop(
        plant_df.columns[~plant_df.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    # drop rows that are NaN and do not have numeric information (e.g: yes, el, heat, etc.)
    plant_df.dropna(subset=["el_output_MW", "heat_output_MW"], how="all", inplace=True)

    # convert electricity column from string to float value
    # TODO: this filtering can be improved
    cols = ["el_output_MW", "heat_output_MW"] # after renaming

    for col in cols:
        plant_df = plant_df[
            (plant_df[col].astype(str).str.contains("MW")) | plant_df[col].isna()
        ]
        
        plant_df[col] = (
            plant_df[col]
            .astype(str)
            .str.extract("(\\d+)")
            .astype(float)
        )

    return plant_df


def prep_pp_df(gens_df, plant_df): # Modified, combine both generators and plants
    """
    Prepare power generation dataframe.
    """
    # reset index and convert the polygon shapes into centroids
    gens_df = gens_df.reset_index(drop=False)
    # gens_df = gens_df.to_crs("EPSG:3035") # https://epsg.io/3035 used in EU
    gens_df = gens_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM    
    gens_df['geometry'] = gens_df['geometry'].centroid
    gens_df = gens_df.to_crs("EPSG:4326")
    
    plant_df = plant_df.reset_index(drop=False)
    # plant_df = plant_df.to_crs("EPSG:3035") # https://epsg.io/3035 used in EU
    plant_df = plant_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    plant_df['geometry'] = plant_df['geometry'].centroid
    plant_df = plant_df.to_crs("EPSG:4326")

    check_fields_gens = ["generator:output:electricity", 
                    "generator:output:hot_water",
                    # "generator:output:biogas" # not included
                   ]

    for field_to_add in check_fields_gens:
        if field_to_add not in gens_df.columns.tolist():
            gens_df[field_to_add] = ""

    check_fields_plants = ["plant:output:electricity", 
                    "plant:output:hot_water",
                    # "plant:output:cold_water" # not included
                   ]

    for field_to_add in check_fields_plants:
        if field_to_add not in plant_df.columns.tolist():
            plant_df[field_to_add] = ""

    # Add longitude/latitude
    gens_df["lon"] = gens_df["geometry"].x
    gens_df["lat"] = gens_df["geometry"].y

    plant_df["lon"] = plant_df["geometry"].x
    plant_df["lat"] = plant_df["geometry"].y

    # Modification - create final dataframe layout
    gens_df = gens_df.rename(
        columns={
            "osmid": "pp_id",
            "generator:method": "method",
            "generator:output:electricity": "el_MW",
            "generator:output:hot_water": "heat_MW",
            "generator:source": "gen_source",
            "geometry": "geometry",
            "operator": "operator",
            "name": "name",
        }
    )

    plant_df = plant_df.rename(
        columns={
            "osmid": "pp_id",
            "plant:output:electricity": "el_MW",
            "plant:output:hot_water": "heat_MW",
            # "plant:output:cold_water": "cold_output_MW",
            "plant:source": "plant_source",
            "plant:method": "method",
            "geometry": "geometry",
            "operator": "operator",
            "name": "name",
        }
    )

    # Rearrange columns
    clist = [
        "pp_id",
        "el_MW",
        "heat_MW",
        "gen_source",
        "plant_source",
        "method",
        "lon",
        "lat",
        "geometry",
        "operator",
        "name",
    ]

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in gens_df:
            gens_df[c] = np.nan

    for c in clist:
        if c not in plant_df:
            plant_df[c] = np.nan

    gens_df.drop(
        gens_df.columns[~gens_df.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    plant_df.drop(
        plant_df.columns[~plant_df.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    # drop rows that are NaN and do not have numeric information (e.g: yes, el, heat, etc.)
    gens_df.dropna(subset=["el_MW", "heat_MW"], how="all", inplace=True)
    plant_df.dropna(subset=["el_MW", "heat_MW"], how="all", inplace=True)

    # convert electricity column from string to float value
    # TODO: this filtering can be improved
    cols = ["el_MW", "heat_MW"] # after renaming

    for col in cols:
        gens_df = gens_df[
            (gens_df[col].astype(str).str.contains("MW")) | gens_df[col].isna()
        ]
        
        gens_df[col] = (
            gens_df[col]
            .astype(str)
            .str.extract("(\\d+)")
            .astype(float)
        )

        plant_df = plant_df[
            (plant_df[col].astype(str).str.contains("MW")) | plant_df[col].isna()
        ]
        
        plant_df[col] = (
            plant_df[col]
            .astype(str)
            .str.extract("(\\d+)")
            .astype(float)
        )

    # gens_df = gens_df[clist]
    # plant_df = plant_df[clist]

    pp_df = pd.concat([gens_df, plant_df], axis=0).reset_index(drop=True)
    pp_df = pp_df[clist]

    return pp_df


def add_line_endings_tosubstations(substations, lines):
    if lines.empty:
        return substations

    # extract columns from substation df
    bus_s = gpd.GeoDataFrame(columns=substations.columns, crs=substations.crs)
    bus_e = gpd.GeoDataFrame(columns=substations.columns, crs=substations.crs)

    # Read information from line.csv
    # bus_s[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_s["voltage"] = lines["voltage"].astype(str)
    bus_s["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[0] if len(p.geoms) >= 2 else None
    )
    bus_s["lon"] = bus_s["geometry"].map(lambda p: p.x if p != None else None)
    bus_s["lat"] = bus_s["geometry"].map(lambda p: p.y if p != None else None)
    bus_s["bus_id"] = (
        (substations["bus_id"].max() if "bus_id" in substations else 0)
        + 1
        + bus_s.index
    )
    # bus_s["dc"] = lines["dc"] # not relevant here

    # bus_e[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_e["voltage"] = lines["voltage"].astype(str)
    bus_e["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[1] if len(p.geoms) >= 2 else None
    )
    bus_e["lon"] = bus_e["geometry"].map(lambda p: p.x if p != None else None)
    bus_e["lat"] = bus_e["geometry"].map(lambda p: p.y if p != None else None)
    bus_e["bus_id"] = bus_s["bus_id"].max() + 1 + bus_e.index
    # bus_e["dc"] = lines["dc"] # not relevant here

    bus_all = pd.concat([bus_s, bus_e], ignore_index=True)

    # Initialize default values
    bus_all["station_id"] = np.nan
    # # Assuming substations completed for installed lines
    # bus_all["under_construction"] = False # not relevant here
    # bus_all["tag_area"] = 0.0 # not relevant here
    bus_all["symbol"] = "substation"
    # TODO: this tag may be improved, maybe depending on voltage levels
    # bus_all["tag_substation"] = "transmission"
    bus_all["substation"] = "transmission"

    buses = pd.concat([substations, bus_all], ignore_index=True)

    # Assign index to bus_id
    buses["bus_id"] = buses.index

    return buses


def clean_voltage(df): # Modified
    """
    Clean the raw voltage column: manual fixing and drop nan values
    """
    # # replace raw values

    repl_voltage = {
        # "medium": "33000",
        # "19.1 kV": "19100",
        # "high": "220000",
        # "240 VAC": "240",
        # "2*220000": "220000;220000",
        # "KV30": "30kV",
        "130;40;10": "130000;40000;10000", # only relevant one here
    }

    df.dropna(subset=["voltage"], inplace=True)

    df["voltage"] = (
        df["voltage"]
        .astype(str)
        .str.lower()
        .str.replace(" ", "")
        .str.replace("_", "")
        .str.replace("kv", "000")
        .str.replace("v", "")
        .replace(repl_voltage)
        # .str.replace("/", ";")  # few OSM entries are separated by / instead of ;
        # this line can be a fix for that if relevant
    )

    return df


def clean_lines(df): # Modified
    """
    Clean the raw circuits/lines/cables column: manual fixing and clean nan values
    """
    # replace raw values
    repl_circuits = {
        "1/3": "1",
        "2/3": "2",
        # assumption of two lines and one grounding wire
        "2-1": "2",
        "single": "1",
        "partial": "1",
        "1;1 disused": "1;0",
        # assuming that in case of a typo there is at least one line
        "`": "1",
        "^1": "1",
        "e": "1",
        "d": "1",
        "1.": "1",
    }

    repl_cables = {
        "1 disused": "0",
        "ground": "0",
        "single": "1",
        "triple": "3",
        "3;3 disused": "3;0",
        "1 (Looped - Haul & Return) + 1 power wire": "1",
        "2-1": "3",
        "3+3": "6",
        "6+1": "6",
        "2x3": "6",
        "3x2": "6",
        "2x2": "4",
        # assuming that in case of a typo there is at least one line
        "partial": "1",
        "`": "1",
        "^1": "1",
        "e": "1",
        "d": "1",
        "line": "1",
    }

    # ASSUMPTION that an empty cable info has at least 3 cables
    if df['cables'].isna().any:
        # df['cables'].fillna(3, inplace=True)
        df.fillna({'cables' : 3}, inplace=True)
    
    # remove rows that has NaN in all circuits, wires, cables columns
    df.dropna(subset=["circuits", "wires", "cables"], how="all", inplace=True)

    # note: no string conversion here! it is performed later on
    df["circuits"] = (
        df["circuits"]
        # .replace(repl_circuits)
        .map(lambda x: x.replace(" ", "") if isinstance(x, str) else x)
    )

    df["cables"] = df["cables"].map(
        lambda x: x.replace(" ", "") if isinstance(x, str) else x
    )

    return df


def set_unique_id(df, col):
    """
    Create unique id's, where id is specified by the column "col" The steps
    below create unique bus id's without losing the original OSM bus_id.

    Unique bus_id are created by simply adding -1,-2,-3 to the original bus_id
    Every unique id gets a -1
    If a bus_id exist i.e. three times it it will the counted by cumcount -1,-2,-3 making the id unique

    Parameters
    ----------
    df : dataframe
        Dataframe considered for the analysis
    col : str
        Column name for the analyses; examples: "bus_id" for substations or "line_id" for lines
    """
    # operate only if id is not already unique (nunique counts unique values)
    if df[col].count() != df[col].nunique():
        # create cumcount column. Cumcount counts 0,1,2,3 the number of duplicates
        df["cumcount"] = df.groupby([col]).cumcount()
        # avoid 0 value for better understanding
        df["cumcount"] = df["cumcount"] + 1
        # add cumcount to id to make id unique
        df[col] = df[col].astype(str) + "-" + df["cumcount"].values.astype(str)
        # remove cumcount column
        df.drop(columns="cumcount", inplace=True)

    return df


def split_cells(df, cols=["voltage"]):
    """
    Split semicolon separated cells i.e. [66000;220000] and create new
    identical rows.

    Parameters
    ----------
    df : dataframe
        Dataframe under analysis
    cols : list
        List of target columns over which to perform the analysis

    Example
    -------
    Original data:
    row 1: '66000;220000', '50'

    After applying split_cells():
    row 1, '66000', '50'
    row 2, '220000', '50'
    """
    if df.empty:
        return df

    x = df.assign(**{col: df[col].str.split(";") for col in cols})

    return x.explode(cols, ignore_index=True)


def filter_voltage(df, bot_treshold, up_threshold): # Modified
    """
    Filters df to contain only assets within voltage threshold.
    """
    # Convert to numeric and drop any row with N/A voltage
    df["voltage"] = pd.to_numeric(df["voltage"], errors="coerce").astype(float)
    df.dropna(subset=["voltage"], inplace=True)

    # convert voltage to int
    df["voltage"] = df["voltage"].astype(int)

    # drop assets outside of threshold_voltage
    df.drop(
        df[(df.voltage < bot_treshold) | (df.voltage > up_threshold)].index,
        axis=0,
        inplace=True,
        errors="ignore",
    )

    return df


def filter_by_geometry(df, region): # Modified
    if df.empty:
        return df
    # drop None geometries
    df.dropna(subset=["geometry"], axis=0, inplace=True)

    # projection from lon lat to distance
    # df = df.to_crs("EPSG:3857")
    df = df.to_crs("EPSG:3006")
    # region = region.to_crs("EPSG:3857")
    region = region.to_crs("EPSG:3006")

    # remove lines without endings (Temporary fix for a Tanzanian line TODO: reformulation?)
    df = df[
        df["geometry"].map(lambda g: len(g.boundary.geoms) >= 2)
    ]
    
    # columns to keep in the lines df
    cols_to_keep = len(df.columns) 

    # remove lines not within the region of interest
    df = gpd.sjoin(df, region, predicate='within')
    # df = gpd.sjoin_nearest(df, region)
    df = df.iloc[:, 0:cols_to_keep+1]
    df.drop(['index_right'], axis=1, inplace=True)

    # projecting back to lon lat
    df = df.to_crs("EPSG:4326")

    return df


def assign_region(df, region): # Modified
    if df.empty:
        return df
    # drop None geometries
    df.dropna(subset=["geometry"], axis=0, inplace=True)

    # projection from lon lat to distance
    # df = df.to_crs("EPSG:3857")
    df = df.to_crs("EPSG:3006")
    # region = region.to_crs("EPSG:3857")
    region = region.to_crs("EPSG:3006")
    
    # columns to keep in the lines df
    cols_to_keep = len(df.columns)

    # # remove lines not within the region of interest
    # df = gpd.sjoin(df, region, predicate='within')
    df = gpd.sjoin_nearest(df, region)
    df = df.iloc[:, 0:cols_to_keep+3] # to include KnKod and KnNamn
    df.drop(['index_right'], axis=1, inplace=True)

    # projecting back to lon lat
    df = df.to_crs("EPSG:4326")

    return df


def finalize_subs(subs_df):
    """
    Specify bus_id and voltage columns as integer.
    """
    subs_df["bus_id"] = subs_df["bus_id"].astype(int)
    subs_df["voltage"] = subs_df["voltage"].astype(int)

    return subs_df


def finalize_lines(lines_df):
    """
    Specify line_id and voltage as integer.
    """
    lines_df["line_id"] = lines_df["line_id"].astype(int)
    lines_df["voltage"] = lines_df["voltage"].astype(int)

    return lines_df


def maptofolium(reg, subs, pp, lines):
    """
    Mapping with folium package
    """
    # define base layer, location as per VGR centroid 
    # https://latlong.info/sweden/vastra-gotaland-county
    regmap = fl.Map(location=[58.25279260, 13.05964250], zoom_start=7)
    
    # add markers to the map
    # currently collect kommun borders, substations, power plants, and power lines 
    fl.GeoJson(
        reg['geometry'],
        style_function=lambda feature: {
        "fillColor": "orange",
        "color": "black",
        "weight": 1,
        "dashArray": "5, 5"
        }
    ).add_to(regmap)
    
    fl.GeoJson(
        subs['geometry'],
        name="Substation",
        marker=fl.Circle(radius=10, fill_color="black", fill_opacity=0.4, color="black", weight=5),
    ).add_to(regmap)
    
    fl.GeoJson(
        pp['geometry'],
        name="Generation",
        marker=fl.Circle(radius=6, fill_color="red", fill_opacity=0.4, color="red", weight=2),
    ).add_to(regmap)
    
    fl.GeoJson(
        lines['geometry'],
        style_function=lambda feature: {
            "color": "blue",
            "weight": 1
        }
    ).add_to(regmap)
    
    return regmap


def collect_from_osm(area):
    """
    Retrieve raw infrastructure from OSM
    Download the grid data using OSMNX library
    """
    subs_df = ox.features_from_place(
                    query=area, 
                    tags={"power": "substation"},
                )

    lines_df = ox.features_from_place(
                    query=area, 
                    tags={"power": ["line", "cable"]},
                )

    gens_df = ox.features_from_place(
                    query=area, 
                    tags={"power": ["generator"]},
                )

    plants_df = ox.features_from_place(
                    query=area, 
                    tags={"power": ["plant"]},
                )

    return subs_df, lines_df, gens_df, plants_df


def retrieve_osm(
    rawlines_df,    # lines data
    rawsubs_df,     # subs data
    rawgens_df,     # gens tag data
    rawplants_df,   # plants tag data
    region,         # region area, for lines spatial join
    kommun,         # kommun area, for subs and pp spatial join, give kommun cols
    bot_thres,      # lower voltage threshold for voltage filter 
    up_thres        # upper voltage threshold for voltage filter
):
    """
    Sequence to prepare df and clean the OSM infra data
    """

    lines_df = prep_lines_df(rawlines_df)
    lines_df = clean_voltage(lines_df)
    lines_df = clean_lines(lines_df)
    lines_df = gpd.GeoDataFrame(
                split_cells(pd.DataFrame(lines_df)),
                crs=lines_df.crs,
            )
    lines_df = filter_voltage(lines_df, bot_thres, up_thres)
    lines_df = filter_by_geometry(lines_df, region)
    lines_df = finalize_lines(lines_df)
    lines_df = set_unique_id(lines_df, "line_id")

    subs_df = prep_subs_df(rawsubs_df)
    subs_df = clean_voltage(subs_df)
    subs_df = set_unique_id(subs_df, "bus_id")
    subs_df = gpd.GeoDataFrame(
                split_cells(pd.DataFrame(subs_df)),
                crs=subs_df.crs,
            )
    subs_df = add_line_endings_tosubstations(
                    subs_df, lines_df
                )
    subs_df = filter_voltage(subs_df, bot_thres, up_thres)
    subs_df = assign_region(subs_df, kommun)
    subs_df = finalize_subs(subs_df)
    subs_df = set_unique_id(subs_df, "bus_id")

    pp_df = prep_pp_df(rawgens_df, rawplants_df)
    pp_df = assign_region(pp_df, kommun)

    return lines_df, subs_df, pp_df


"""
The following are the sequence to retrieve and clean the OSM data
"""

if __name__ == "__main__":

    input_folder = prerequisites("inputs")

    rawsubs_df, rawlines_df, rawgens_df, rawplants_df = collect_from_osm('Västra Götaland')

    lans_file = getinputfiles('Lan_Sweref99TM_region.shp', input_folder)
    kommuns_file = getinputfiles('Kommun_Sweref99TM.shp', input_folder)

    lans = gpd.read_file(lans_file)
    kommuns = gpd.read_file(kommuns_file)

    vgr_reg = ox.geocode_to_gdf('Västra Götaland')

    vgr_lan = reg_bound(lans)
    vgr_kommun = reg_bound(kommuns)

    lines_df, subs_df, pp_df = retrieve_osm(rawlines_df, 
                                            rawsubs_df, 
                                            rawgens_df, 
                                            rawplants_df, 
                                            vgr_reg, 
                                            vgr_kommun, 
                                            100000, 
                                            200000
                                )

    savefile(lines_df, "lines_raw", input_folder)
    savefile(subs_df, "subs_raw", input_folder)
    savefile(pp_df, "pp_raw", input_folder)

    # osm_map = maptofolium(vgr_kommun, rawsubs_df, rawgens_df, rawlines_df)
    # display(osm_map)

    vgr_map = maptofolium(vgr_kommun, subs_df, pp_df, lines_df)
    map_file = os.path.join(input_folder, 'vgr_osm.html')
    vgr_map.save(map_file)
    display(vgr_map)

