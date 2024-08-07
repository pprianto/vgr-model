import os
import warnings
import geopandas as gpd
import numpy as np
import pandas as pd
import folium as fl
from shapely.geometry import LineString, Point
from shapely.ops import linemerge, split
from tqdm import tqdm
from IPython.display import display
from basics import (prerequisites,
                    savefile,
                    reg_bound,
                    getinputfiles)
from unidecode import unidecode

pd.set_option("mode.copy_on_write", True)
pd.options.mode.copy_on_write = True
warnings.simplefilter(action="ignore", category=FutureWarning)

"""
The following list of functions are adapted from
PyPSA-Earth
https://github.com/pypsa-meets-earth/pypsa-earth
https://arxiv.org/abs/2209.04663
"""


def line_endings_to_bus_conversion(lines_df):
    # Assign to every line a start and end point

    lines_df["bounds"] = lines_df["geometry"].boundary  # create start and end point

    lines_df["bus_0_coors"] = lines_df["bounds"].map(lambda p: p.geoms[0])
    lines_df["bus_1_coors"] = lines_df["bounds"].map(lambda p: p.geoms[1])

    # splits into coordinates
    lines_df["bus0_lon"] = lines_df["bus_0_coors"].x
    lines_df["bus0_lat"] = lines_df["bus_0_coors"].y
    lines_df["bus1_lon"] = lines_df["bus_1_coors"].x
    lines_df["bus1_lat"] = lines_df["bus_1_coors"].y

    return lines_df


# tol in m
def set_substations_ids(subs_df, tol):
    """
    Function to set substations ids to buses, accounting for location
    tolerance.

    The algorithm is as follows:

    1. initialize all substation ids to -1
    2. if the current substation has been already visited [substation_id < 0], then skip the calculation
    3. otherwise:
        1. identify the substations within the specified tolerance (tol)
        2. when all the substations in tolerance have substation_id < 0, then specify a new substation_id
        3. otherwise, if one of the substation in tolerance has a substation_id >= 0, then set that substation_id to all the others;
           in case of multiple substations with substation_ids >= 0, the first value is picked for all
    """

    subs_df["station_id"] = -1

    # create temporary series to execute distance calculations using m as reference distances
    # temp_bus_geom = subs_df.geometry.to_crs("EPSG:3857") # https://epsg.io/3857 for distance in m
    temp_bus_geom = subs_df.geometry.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM

    # set tqdm options for substation ids
    tqdm_kwargs_substation_ids = dict(
        ascii=False,
        unit=" buses",
        total=subs_df.shape[0],
        desc="Set substation ids ",
    )

    station_id = 0
    for i, row in tqdm(subs_df.iterrows(), **tqdm_kwargs_substation_ids):
        if subs_df.loc[i, "station_id"] >= 0:
            continue

        # get substations within tolerance
        close_nodes = np.flatnonzero(
            temp_bus_geom.distance(temp_bus_geom.loc[i]) <= tol
        )

        if len(close_nodes) == 1:
            # if only one substation is in tolerance, then the substation is the current one iì
            # Note that the node cannot be with substation_id >= 0, given the preliminary check
            # at the beginning of the for loop
            subs_df.loc[subs_df.index[i], "station_id"] = station_id
            # update station id
            station_id += 1
        else:
            # several substations in tolerance
            # get their ids
            subset_substation_ids = subs_df.loc[subs_df.index[close_nodes], "station_id"]
            # check if all substation_ids are negative (<0)
            all_neg = subset_substation_ids.max() < 0
            # check if at least a substation_id is negative (<0)
            some_neg = subset_substation_ids.min() < 0

            if all_neg:
                # when all substation_ids are negative, then this is a new substation id
                # set the current station_id and increment the counter
                subs_df.loc[subs_df.index[close_nodes], "station_id"] = station_id
                station_id += 1
            elif some_neg:
                # otherwise, when at least a substation_id is non-negative, then pick the first value
                # and set it to all the other substations within tolerance
                sub_id = -1
                for substation_id in subset_substation_ids:
                    if substation_id >= 0:
                        sub_id = substation_id
                        break
                subs_df.loc[subs_df.index[close_nodes], "station_id"] = sub_id

    return subs_df


def set_lines_ids(lines_df, subs_df):
    """
    Function to set line buses ids to the closest bus in the list.
    """
    # set tqdm options for set lines ids
    tqdm_kwargs_line_ids = dict(
        ascii=False,
        unit=" lines",
        total=lines_df.shape[0],
        desc="Set line bus ids ",
    )

    # initialization
    lines_df["bus0"] = -1
    lines_df["bus1"] = -1

    # busesepsg = subs_df.to_crs("EPSG:3857") # https://epsg.io/3857 for distance in m
    # linesepsg = lines_df.to_crs("EPSG:3857") # https://epsg.io/3857 for distance in m
    busesepsg = subs_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    linesepsg = lines_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM

    for i, row in tqdm(linesepsg.iterrows(), **tqdm_kwargs_line_ids):
        # select buses having the voltage level of the current line
        buses_sel = busesepsg[
            (subs_df["voltage"] == row["voltage"]) # & (subs_df["dc"] == row["dc"]) # dc not relevant
        ]

        # find the closest node of the bus0 of the line
        bus0_id = buses_sel.geometry.distance(row.geometry.boundary.geoms[0]).idxmin()
        lines_df.loc[i, "bus0"] = subs_df.loc[bus0_id, "station_id"]

        # check if the line starts exactly in the node, otherwise modify the linestring
        distance_bus0 = busesepsg.geometry.loc[bus0_id].distance(
            row.geometry.boundary.geoms[0]
        )
        if distance_bus0 > 0.0:
            # the line does not start in the node, thus modify the linestring
            lines_df.geometry.loc[i] = linemerge(
                [
                    LineString(
                        [
                            subs_df.geometry.loc[bus0_id],
                            lines_df.geometry.loc[i].boundary.geoms[0],
                        ]
                    ),
                    lines_df.geometry.loc[i],
                ]
            )

        # find the closest node of the bus1 of the line
        bus1_id = buses_sel.geometry.distance(row.geometry.boundary.geoms[1]).idxmin()
        lines_df.loc[i, "bus1"] = subs_df.loc[bus1_id, "station_id"]

        # check if the line ends exactly in the node, otherwise modify the linestring
        distance_bus1 = busesepsg.geometry.loc[bus1_id].distance(
            row.geometry.boundary.geoms[1]
        )
        if distance_bus1 > 0.0:
            # the line does not end in the node, thus modify the linestring
            lines_df.geometry.loc[i] = linemerge(
                [
                    lines_df.geometry.loc[i],
                    LineString(
                        [
                            lines_df.geometry.loc[i].boundary.geoms[1],
                            subs_df.geometry.loc[bus1_id],
                        ]
                    ),
                ]
            )

    return lines_df, subs_df


def merge_stations_same_station_id( # Modified
    subs_df, delta_lon=0.001, delta_lat=0.001, precision=4
):
    """
    Function to merge buses with same voltage and station_id This function
    iterates over all substation ids and creates a bus_id for every substation
    and voltage level.

    Therefore, a substation with multiple voltage levels is represented
    with different buses, one per voltage level
    """
    # initialize list of cleaned buses
    buses_clean = []

    # initialize the number of buses
    n_buses = 0

    for g_name, g_value in subs_df.groupby(by="station_id"):
        # average location of the buses having the same station_id
        station_point_x = np.round(g_value.geometry.x.mean(), precision)
        station_point_y = np.round(g_value.geometry.y.mean(), precision)

        # loop for every voltage level in the bus
        # The location of the buses is averaged; in the case of multiple voltage levels for the same station_id,
        # each bus corresponding to a voltage level and each polarity is located at a distance regulated by delta_lon/delta_lat
        v_it = 0
        for v_name, bus_row in g_value.groupby(by=["voltage"]):# , "dc"]): not relevant here
            lon_bus = np.round(station_point_x + v_it * delta_lon, precision)
            lat_bus = np.round(station_point_y + v_it * delta_lat, precision)

            # add the bus
            buses_clean.append(
                [
                    n_buses,  # "bus_id"
                    g_name,  # "station_id"
                    v_name[0],  # "voltage"
                    # bus_row["dc"].all(),  # "dc"
                    "|".join(bus_row["symbol"].unique()),  # "symbol"
                    # bus_row["under_construction"].any(),  # "under_construction"
                    "|".join(bus_row["substation"].unique()),  # "substation"
                    # bus_row["tag_area"].sum(),  # "tag_area"
                    lon_bus,  # "lon"
                    lat_bus,  # "lat"
                    "|".join([str(val) for val in bus_row["import_trans"].unique()]),  # "transmission bus"
                    bus_row["KnKod"].iloc[0],  # "Kommun Kod"
                    bus_row["KnNamn"].iloc[0],  # "Kommun Namn"
                    Point(
                        lon_bus,
                        lat_bus,
                    ),  # "geometry"
                ]
            )

            # increase counters
            v_it += 1
            n_buses += 1

    # names of the columns
    buses_clean_columns = [
        "bus_id",
        "station_id",
        "voltage",
        # "dc",
        "symbol",
        # "under_construction",
        "substation",
        # "area",
        "lon",
        "lat",
        "import_trans",
        "KnKod",
        "KnNamn",
        "geometry",
    ]

    return gpd.GeoDataFrame(buses_clean, columns=buses_clean_columns).set_crs(
        crs=subs_df.crs, inplace=True
    )


# Not used yet, might be relevant later when different voltage levels included
def get_transformers(buses, lines): 
    """
    Function to create fake transformer lines that connect buses of the same
    station_id at different voltage.
    """

    # ac_freq = get_ac_frequency(lines)
    df_transformers = []

    # Transformers should be added between AC buses only
    buses_ac = buses[~buses["dc"]]
    for g_name, g_value in buses_ac.sort_values("voltage", ascending=True).groupby(
        by="station_id"
    ):
        # note: by construction there cannot be more that two buses with the same station_id and same voltage
        n_voltages = len(g_value)

        if n_voltages > 1:
            for id in range(0, n_voltages - 1):
                # when g_value has more than one node, it means that there are multiple voltages for the same bus
                geom_trans = LineString(
                    [g_value.geometry.iloc[id], g_value.geometry.iloc[id + 1]]
                )

                df_transformers.append(
                    [
                        f"transf_{g_name}_{id}",  # "line_id"
                        g_value["bus_id"].iloc[id],  # "bus0"
                        g_value["bus_id"].iloc[id + 1],  # "bus1"
                        g_value.voltage.iloc[id],  # "voltage_bus0"
                        g_value.voltage.iloc[id + 1],  # "voltage_bus0"
                        g_value.country.iloc[id],  # "country"
                        geom_trans,  # "geometry"
                    ]
                )

    # name of the columns
    trasf_columns = [
        "line_id",
        "bus0",
        "bus1",
        "voltage_bus0",
        "voltage_bus1",
        "country",
        "geometry",
    ]

    df_transformers = gpd.GeoDataFrame(df_transformers, columns=trasf_columns)
    if not df_transformers.empty:
        init_index = 0 if lines.empty else lines.index[-1] + 1
        df_transformers.set_index(init_index + df_transformers.index, inplace=True)
    # update line endings
    df_transformers = line_endings_to_bus_conversion(df_transformers)

    return df_transformers


def connect_stations_same_station_id(lines_df, subs_df):
    """
    Function to create fake links between substations with the same
    substation_id.
    """
    # ac_freq = get_ac_frequency(lines)
    station_id_list = subs_df.station_id.unique()

    add_lines = []
    from shapely.geometry import LineString

    for s_id in station_id_list:
        buses_station_id = subs_df[subs_df.station_id == s_id]

        if len(buses_station_id) > 1:
            for b_it in range(1, len(buses_station_id)):
                add_lines.append(
                    [
                        # f"link{buses_station_id}_{b_it}",  # "line_id"
                        f"link{buses_station_id.index[0]}_{buses_station_id.index[b_it]}",  # "line_id"
                        "line", # types
                        3, # cables
                        132000,  # "voltage"
                        None, # wires
                        None, # circuits
                        0.0, # length_km
                        buses_station_id.index[0],  # "bus0"
                        buses_station_id.index[b_it],  # "bus1"
                        LineString(
                            [
                                buses_station_id.geometry.iloc[0],
                                buses_station_id.geometry.iloc[b_it],
                            ]
                        ),  # "geometry"
                        LineString(
                            [
                                buses_station_id.geometry.iloc[0],
                                buses_station_id.geometry.iloc[b_it],
                            ]
                        ).bounds,  # "bounds"
                        buses_station_id.geometry.iloc[0],  # "bus_0_coors"
                        buses_station_id.geometry.iloc[b_it],  # "bus_1_coors"
                        buses_station_id.lon.iloc[0],  # "bus0_lon"
                        buses_station_id.lat.iloc[0],  # "bus0_lat"
                        buses_station_id.lon.iloc[b_it],  # "bus1_lon"
                        buses_station_id.lat.iloc[b_it],  # "bus1_lat"
                    ]
                )

    # name of the columns
    add_lines_columns = [
        "line_id",
        "types",
        "cables",
        "voltage",
        "wires",
        "circuits",
        "length_km",
        "bus0",
        "bus1",
        "geometry",
        "bounds",
        "bus_0_coors",
        "bus_1_coors",
        "bus0_lon",
        "bus0_lat",
        "bus1_lon",
        "bus1_lat",
    ]

    # df_add_lines = gpd.GeoDataFrame(pd.concat(add_lines), columns=add_lines_columns).set_crs(
    #     crs=subs_df.crs, inplace=True
    # )
    df_add_lines = gpd.GeoDataFrame(add_lines, columns=add_lines_columns).set_crs(
        crs=subs_df.crs, inplace=True
    )
    
    lines_df = pd.concat([lines_df, df_add_lines], ignore_index=True)

    return lines_df


def merge_stations_lines_by_station_id_and_voltage(
    lines_df, subs_df, tol
):
    """
    Function to merge close stations and adapt the line datasets to adhere to
    the merged dataset.
    """
    # set substation ids
    # set_substations_ids(subs_df, distance_crs, tol=tol)
    subs_df = set_substations_ids(subs_df, tol)

    # merge buses with same station id and voltage
    # if not subs_df.empty:
    subs_df = merge_stations_same_station_id(subs_df)

    # set the bus ids to the line dataset
    lines_df, subs_df = set_lines_ids(lines_df, subs_df)

    # drop lines starting and ending in the same node
    lines_df.drop(lines_df[lines_df["bus0"] == lines_df["bus1"]].index, inplace=True)
    
    # remove duplicates (if any) based on exactly same bus0 and bus1 and keep the first entry 
    lines_df.drop_duplicates(subset=['bus0', 'bus1'], keep='first', inplace=True)
    
    # drop subs based on station_id existence in bus0 and bus1 in lines_df
    buses = pd.Series(list(lines_df['bus0']) + list(lines_df['bus1']))
    subs_df.drop(subs_df[~subs_df['station_id'].isin(buses.values)].index, inplace=True)

    # update line endings
    lines_df = line_endings_to_bus_conversion(lines_df)

    # set substation_lv
    # set_lv_substations(buses) # not relevant for now

    # append fake converters
    # lines = pd.concat([lines, converters], ignore_index=True)

    # reset index
    lines_df.reset_index(drop=True, inplace=True)
    subs_df.reset_index(drop=True, inplace=True)
    # if len(links) > 0:
    #     links.reset_index(drop=True, inplace=True)

    return lines_df, subs_df


def create_station_at_equal_bus_locations(
    lines_df, subs_df, tol
):
    """
    V1. Create station_id at same bus location
    - We saw that buses are not connected exactly at one point, they are
      usually connected to a substation "area" (analysed on maps)
    - Create station_id at exactly the same location might therefore be not
      always correct
    - Though as you can see below, it might be still sometime the case.
      Examples are **station 4** (2 lines with the same voltage connect at the
      same point) and **station 23** (4 lines with two different voltages connect
      at the same point)
    TODO: Filter out the generator lines - defined as going from generator to
          the next station which is connected to a load. Excluding generator
          lines make probably sense because they are not transmission expansion
          relevant. For now we simplify and include generator lines.

    """

    # If same location/geometry make station
    bus_all = subs_df

    # set substation ids
    subs_df = set_substations_ids(subs_df, tol)

    # set the bus ids to the line dataset
    lines_df, subs_df = set_lines_ids(lines_df, subs_df)

    # update line endings
    lines_df = line_endings_to_bus_conversion(lines_df)

    # For each station number with multiple buses make lowest voltage `substation_lv = TRUE`
    # set_lv_substations(bus_all)

    # TRY: Keep only buses that are not duplicated & lv_substation = True
    # TODO: Check if this is necessary. What effect do duplicates have?
    bus_all = bus_all[bus_all["substation_lv"] == True]

    lines_df = connect_stations_same_station_id(lines_df, subs_df)

    return lines_df, subs_df


def _split_linestring_by_point(linestring, points):
    """
    Function to split a linestring geometry by multiple inner points.

    Parameters
    ----------
    lstring : LineString
        Linestring of the line to be split
    points : list
        List of points to split the linestring

    Return
    ------
    list_lines : list
        List of linestring to split the line
    """

    list_linestrings = [linestring]

    for p in points:
        # execute split to all lines and store results
        temp_list = [split(l, p) for l in list_linestrings]
        # nest all geometries
        list_linestrings = [lstring for tval in temp_list for lstring in tval.geoms]

    return list_linestrings


def fix_overpassing_lines(lines, buses, tol):
    """
    Function to avoid buses overpassing lines with no connection when the bus
    is within a given tolerance from the line.

    Parameters
    ----------
    lines : GeoDataFrame
        Geodataframe of lines
    buses : GeoDataFrame
        Geodataframe of substations
    tol : float
        Tolerance in meters of the distance between the substation and the line
        below which the line will be split
    """

    lines_to_add = []  # list of lines to be added
    lines_to_split = []  # list of lines that have been split

    # lines_epsgmod = lines.to_crs(distance_crs)
    # buses_epsgmod = buses.to_crs(distance_crs)
    lines_epsgmod = lines.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    buses_epsgmod = buses.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM


    # set tqdm options for substation ids
    tqdm_kwargs_substation_ids = dict(
        ascii=False,
        unit=" lines",
        total=lines.shape[0],
        desc="Verify lines overpassing nodes ",
    )

    for l in tqdm(lines.index, **tqdm_kwargs_substation_ids):
        # bus indices being within tolerance from the line
        bus_in_tol_epsg = buses_epsgmod[
            buses_epsgmod.geometry.distance(lines_epsgmod.geometry.loc[l]) <= tol
        ]

        # exclude endings of the lines
        bus_in_tol_epsg = bus_in_tol_epsg[
            (
                (
                    bus_in_tol_epsg.geometry.distance(
                        lines_epsgmod.geometry.loc[l].boundary.geoms[0]
                    )
                    > tol
                )
                | (
                    bus_in_tol_epsg.geometry.distance(
                        lines_epsgmod.geometry.loc[l].boundary.geoms[1]
                    )
                    > tol
                )
            )
        ]

        if not bus_in_tol_epsg.empty:
            # add index of line to split
            lines_to_split.append(l)

            buses_locs = buses.geometry.loc[bus_in_tol_epsg.index]

            # get new line geometries
            new_geometries = _split_linestring_by_point(lines.geometry[l], buses_locs)
            n_geoms = len(new_geometries)

            # create temporary copies of the line
            df_append = gpd.GeoDataFrame([lines.loc[l]] * n_geoms)
            # update geometries
            df_append["geometry"] = new_geometries
            # update name of the line
            df_append["line_id"] = [
                str(df_append["line_id"].iloc[0]) + f"_{id}" for id in range(n_geoms)
            ]

            lines_to_add.append(df_append)

    if not lines_to_add:
        return lines, buses

    df_to_add = gpd.GeoDataFrame(pd.concat(lines_to_add, ignore_index=True))
    df_to_add.set_crs(lines.crs, inplace=True)
    df_to_add.set_index(lines.index[-1] + df_to_add.index, inplace=True)

    # # update length # already done in retrieve_osm
    # df_to_add["length"] = df_to_add.to_crs("EPSG:3006").geometry.length / 1000 # in km

    # update line endings
    df_to_add = line_endings_to_bus_conversion(df_to_add)

    # remove original lines
    lines.drop(lines_to_split, inplace=True)

    lines = gpd.GeoDataFrame(
        pd.concat([lines, df_to_add], ignore_index=True).reset_index(drop=True),
        crs=lines.crs,
    )

    return lines, buses


def cleanbuses(subs_df):
    """
    Final cleaning such as:
     1. substation into transmission only
     2. import transmission node as True/False
     3. set voltage to kV
     4. Create a column of node_id
     5. Set the node_id according to municipality name
     6. 3 capital letters and numbered index
     7. rearrange columns
    """

    # Set substation cols into transmission
    subs_df["voltage"] /= 1e3 # in kV
    subs_df["substation"] = "transmission"
    subs_df['import_trans'] = (subs_df['import_trans']
                               .apply(lambda x: True if pd.notna(x) and '1' in x else False))

    # Change datatype for further purposes
    subs_df["symbol"] = subs_df["symbol"].astype('str')
    subs_df["substation"] = subs_df["substation"].astype('str')
    subs_df["KnKod"] = subs_df["KnKod"].astype('int64')
    subs_df["KnNamn"] = subs_df["KnNamn"].astype('str')
    subs_df["municipality"] = subs_df["KnNamn"].apply(unidecode)

    # create node_id column
    subs_df["helper_col"] = subs_df["municipality"].apply(base_node_id)

    subs_df["node_id"] = subs_df.groupby("helper_col").cumcount() + 1
    subs_df["node_id"] = subs_df["helper_col"] + subs_df["node_id"].astype(str)

    subs_df = subs_df.drop(columns=["helper_col"])

    # column rearrange
    subs_df = subs_df[[
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
                    "geometry"
                ]]

    # subs_df = subs_df.rename(
    #     columns={
    #         "KnKod": "kn_kod",
    #         "KnNamn": "kn_namn",
    #     }
    # )    

    return subs_df

def base_node_id(municipality):
    knlist = {
        "Goteborg" : "GBG"
        # add more if any unique cases
    }

    if municipality in knlist:
        return knlist[municipality]
    else:
        return municipality[:3].upper() # needs 3 letters
    

def cleanlines(lines_df, subs_df):
    """
    Final cleaning such as:
     1. fill wires and circuits rows assuming the cables value
     2. set voltage to 130kV
     3. add lines_id column
     4. add node_from node_id column
     5. node_from and node_id based on the substation node_id
    """

    # Set substation cols into transmission
    lines_df["voltage"] /= 1e3 # in kV
   
    lines_df.loc[lines_df['cables'] == '3', ['wires', 'circuits']] = ['single', 1]
    lines_df.loc[lines_df['cables'] == '6', ['wires', 'circuits']] = ['double', 2]

    # # Change datatype for further purposes
    lines_df["types"] = lines_df["types"].astype('str')
    lines_df["cables"] = lines_df["cables"].astype('int64')
    lines_df["wires"] = lines_df["wires"].astype('str')
    lines_df["circuits"] = lines_df["circuits"].astype('int64')

    # new column for lines id instead of osm_id
    lines_df["lines_id"] = ["L" + str(i + 1) for i in range(len(lines_df))]

    lines_df = getnodesforlines(lines_df, subs_df)

    # column rearrange
    lines_df = lines_df[[
                    "line_id",
                    "lines_id",
                    "types",
                    "voltage",
                    "cables",
                    "wires",
                    "circuits",
                    "node_from",
                    "bus0",
                    "bus0_lon",
                    "bus0_lat",
                    "node_to",
                    "bus1",
                    "bus1_lon",
                    "bus1_lat",
                    "length_km",
                    "geometry"
                ]]
    
    lines_df = lines_df.rename(
        columns={
            "line_id": "osm_id",
            "bus0": "station_from",
            "bus0_lon": "node_from_lon", 
            "bus0_lat": "node_from_lat",
            "bus1": "station_to",
            "bus1_lon": "node_to_lon",
            "bus1_lat": "node_to_lat"
        }
    )    

    return lines_df

def getnodesforlines(lines_df, subs_df):

    station_to_node = subs_df.set_index("station_id")["node_id"].to_dict()

    lines_df["node_from"] = lines_df["bus0"].map(station_to_node)
    lines_df["node_to"] = lines_df["bus1"].map(station_to_node)

    return lines_df



def cleanpp(pp_df):
    """
    Final cleaning such as:
     1. remove pp that has no el or heat output
     2. create tech column that consists of gen and plant source column
     3. column rearrange
    """

    # Set 0 values of output into Nan for removal
    powcols = ['el_MW', 'heat_MW']
    pp_df[powcols] = pp_df[powcols].replace(0, np.nan)

    pp_df = pp_df.dropna(subset = powcols, how='all')
    
    # replace NaN in the source columns with empty strings
    fill_gen_source = pp_df["gen_source"].fillna('')
    fill_plant_source = pp_df["plant_source"].fillna('')

    pp_df["tech"] = fill_gen_source + ";" + fill_plant_source
    pp_df['tech'] = pp_df['tech'].str.strip('; ').replace('; ;', ';')

    pp_df["municipality"] = pp_df["KnNamn"].apply(unidecode)

    # column rearrange
    pp_df = pp_df[[
                    "pp_id",
                    "el_MW",
                    "heat_MW",
                    "tech",
                    "gen_source",
                    "plant_source",
                    "method",
                    "operator",
                    "name",
                    "KnKod",
                    "KnNamn",
                    "municipality",
                    "lon",
                    "lat",
                    "geometry"
                ]]

    # pp_df = pp_df.rename(
    #     columns={
    #         "KnKod": "kn_kod",
    #         "KnNamn": "kn_namn",
    #     }
    # )    

    pp_df.drop(columns=["gen_source", "plant_source"], inplace=True)

    return pp_df


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
        subs,
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
    ).add_to(regmap)
    
    fl.GeoJson(
        pp,
        name="Generation",
        marker=fl.Circle(
                         radius=6, 
                         fill_color="red", 
                         fill_opacity=0.4, 
                         color="red", 
                         weight=2
        ),
        tooltip=fl.GeoJsonTooltip(
            fields=[
                    "el_MW", 
                    "heat_MW", 
                    "tech"
        ]),
    ).add_to(regmap)
    
    fl.GeoJson(
        lines,
        name="Power Lines",
        style_function=lambda feature: {
            "color": "blue",
            "weight": 1
        }
    ).add_to(regmap)
    
    return regmap


def built_network(
    lines,
    buses,
    tol_merge, # 2000
    tol_overpass # 1
):

    # Stage 1: Assign start/end buses for lines
    lines = line_endings_to_bus_conversion(lines)

    # Stage 2: Fix overpassing lines
    lines, buses = fix_overpassing_lines(lines, buses, tol_overpass)

    # # Add bus to countries with no buses
    # buses = add_buses_to_empty_countries(countries_config, inputs.country_shapes, buses)

    # Stage 3: Merging substations by distance tolerance
    lines, buses = merge_stations_lines_by_station_id_and_voltage(
            lines, buses, tol_merge
        )

    # # get transformers: modelled as lines connecting buses with different voltage
    # transformers = get_transformers(buses, lines)

    # # get converters: currently modelled as links connecting buses with different polarity
    # converters = get_converters(buses, lines)

    # return None
    return lines, buses

"""
The following are the sequence to build the el network
"""

if __name__ == "__main__":

    # current_dir = os.getcwd()
    # subfolder = "inputs"
    # input_folder = os.path.join(current_dir, subfolder)
    # mode = 0o666
    # os.makedirs(input_folder, mode, exist_ok=True)

    input_folder = prerequisites("inputs")

    osmlines_fn = getinputfiles("lines_raw.geojson", input_folder)
    osmsubs_fn = getinputfiles("subs_raw.geojson", input_folder)
    osmpp_fn = getinputfiles("pp_raw.geojson", input_folder)

    osmlines_df = gpd.read_file(osmlines_fn)
    osmsubs_df = gpd.read_file(osmsubs_fn)
    osmpp_df = gpd.read_file(osmpp_fn)

    lines_df, subs_df = built_network(osmlines_df, osmsubs_df, 2000, 1)

    lans_file = getinputfiles('Lan_Sweref99TM_region.shp', input_folder)
    kommuns_file = getinputfiles('Kommun_Sweref99TM.shp', input_folder)
    lans = gpd.read_file(lans_file)
    kommuns = gpd.read_file(kommuns_file)
    vgr_lan = reg_bound(lans)
    vgr_kommun = reg_bound(kommuns)

    # vgr_map = maptofolium(vgr_kommun, subs_df, osmpp_df, lines_df)
    # map_file = os.path.join(input_folder, 'vgr_clean.html')
    # vgr_map.save(map_file)


    subs_clean = cleanbuses(subs_df)
    lines_clean = cleanlines(lines_df, subs_clean)
    pp_clean = cleanpp(osmpp_df)


    # so that the index starts from 1 before going into Julia
    subs_clean['bus_id'] = subs_clean['bus_id'] + 1
    subs_clean['station_id'] = subs_clean['station_id'] + 1

    lines_clean['station_from'] = lines_clean['station_from'] + 1
    lines_clean['station_to'] = lines_clean['station_to'] + 1

    savefile(subs_clean, "subs_clean", input_folder)
    savefile(lines_clean, "lines_clean", input_folder)
    savefile(pp_clean, "pp_clean", input_folder)

    vgr_map = maptofolium(vgr_kommun, subs_clean, pp_clean, lines_clean)
    map_file = os.path.join(input_folder, 'vgr_clean.html')
    vgr_map.save(map_file)


