import os
import geopandas as gpd
import numpy as np
import pandas as pd
import folium as fl
import polars as pl
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

"""
Functions to appoint demand on nodes
1. assign nearest nodes for missing municipalities
2. annual demand based on SCB
2. profile based on GE 2019 data with randomisation
3. ...
"""

def get_nodes(subs_df, kommun_df):
    '''
    function to assign missing kommun to the nearest node
    '''

    # distance projection and convert shape of municipalities to point
    subs_df = subs_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    kommun_df = kommun_df.to_crs("EPSG:3006") # https://epsg.io/3006 used in Sweden / SWEREF99 TM
    kommun_df["geometry"] = kommun_df["geometry"].centroid
    
    # dummy dfs for missing and existing municipalities
    missing_kn = kommun_df[~kommun_df["KnNamn"].isin(subs_df["KnNamn"])]
    existing_kn = kommun_df[kommun_df["KnNamn"].isin(subs_df["KnNamn"])]

    # get nearest node
    missing_kn["nearest_idx"] = missing_kn["geometry"].apply(lambda geom: subs_df.distance(geom).idxmin())
    near_nodes = subs_df.loc[missing_kn["nearest_idx"]].reset_index(drop=True)

    # add new columns
    missing_kn["new_KnKod"] = near_nodes["KnKod"].values
    missing_kn["new_KnNamn"] = near_nodes["KnNamn"].values
    missing_kn["new_node_id"] = near_nodes["node_id"].values

    # assign knkod/namn and node to unchanged municipalities
    existing_kn = existing_kn.merge(
        subs_df[["KnKod", "KnNamn", "node_id"]],
        on="KnNamn",
        how="left",
        suffixes=('', '_subs_df')
    )

    existing_kn['new_KnKod'] = existing_kn['KnKod']
    existing_kn['new_KnNamn'] = existing_kn['KnNamn']
    existing_kn['new_node_id'] = existing_kn['node_id']

    # combine changed and unchanged municipalities
    kommun_df = pd.concat(
        [
            missing_kn[["KnKod", "KnNamn", "geometry", "new_KnKod", "new_KnNamn", "new_node_id"]],
            existing_kn[["KnKod", "KnNamn", "geometry", "new_KnKod", "new_KnNamn", "new_node_id"]]
        ],
        ignore_index=True
    )

    kommun_df = kommun_df.to_crs("EPSG:4326") # https://epsg.io/4326 World Geodetic System 1984

    # column rearrange
    kommun_df = kommun_df[[
                    "KnKod",
                    "KnNamn",
                    "new_KnKod",
                    "new_KnNamn",
                    "new_node_id",
                    "geometry"
                ]]

    type_dict = {
                "KnKod" : int,
                "KnNamn": str,
                "new_KnKod": int,
                "new_KnNamn": str,
                "new_node_id": str,

    }

    kommun_df = kommun_df.astype(type_dict)

    kommun_df.sort_values(by=["KnKod"], inplace=True)
    kommun_df.reset_index(drop=True, inplace=True)

    return kommun_df


def profile_kommun(profile_df, sumdemand_df):
    '''
    function to create excel of kommun hourly demand based on
    Goteborg energi profile
    '''

    # empty df to save as excel and further processing
    demand_df = pd.DataFrame()
    demand_df["hour"] = profile_df["hour"]

    # create dataframe for further processing
    for _, row in sumdemand_df.iterrows():
        knkod = row["KnKod"]
        sum_demand = row["all"]

        hr_demand = profile_df["scale"] * sum_demand

        demand_df[knkod] = hr_demand

    return demand_df


def profile_nodes(kn_df, profile_df):
    '''
    function to assign municipal demand to nodal demand
    '''

    # column for number of nodes
    kn_df["nodes_count"] = kn_df.groupby("KnKod")["new_node_id"].transform("count")

    # df with knkod column for demand over time
    profile_melted = profile_df.melt(
                        id_vars=["hour"], 
                        var_name="KnKod", 
                        value_name="demand"
                        )

    # merge the melted with kommun df
    merged_df = profile_melted.merge(
                    kn_df, 
                    on=["KnKod"], 
                    how="left"
                    )

    # adjust the demand according to no of nodes 
    merged_df["demand"] = merged_df["demand"] / merged_df["nodes_count"]

    # pivot the table according to new_node_id as columns
    node_df = merged_df.pivot_table(
                        index="hour", 
                        columns="new_node_id", 
                        values="demand",
                        aggfunc="sum"
                        ).reset_index()
  
    node_df = node_df.rename_axis(None, axis=1)

    return node_df


def group_demand(
    subs_df,
    kommuns,
    profile,
    annual_demand
):

    kommun_nodes = get_nodes(subs_df, kommuns)
    

    kommun_demand = profile_kommun(profile, annual_demand)
    

    nodal_demand = profile_nodes(kommun_nodes, kommun_demand)
    

    return kommun_nodes, kommun_demand, nodal_demand


if __name__ == "__main__":

    input_folder = prerequisites("inputs")
    model_input = prerequisites("modelinput")

    # read substation
    subs_fn = getinputfiles("subs_clean.geojson", input_folder)
    subs_df = gpd.read_file(subs_fn)

    # get municipality
    kommuns_file = getinputfiles('Kommun_Sweref99TM.shp', input_folder)
    kommuns = gpd.read_file(kommuns_file)
    vgr_kommun = reg_bound(kommuns)

    # load profile and demand for municipality
    el_profile_fn = getinputfiles("el_profile.xlsx", input_folder)
    el_demand_fn = getinputfiles("VGR_eldemand.xlsx", input_folder)
    heat_profile_fn = getinputfiles("heat_profile.xlsx", input_folder)
    heat_demand_fn = getinputfiles("VGR_heatdemand.xlsx", input_folder)
    h2_profile_fn = getinputfiles("h2_profile.xlsx", input_folder)
    h2_demand_fn = getinputfiles("VGR_h2demand.xlsx", input_folder)

    el_profile = pd.read_excel(el_profile_fn, sheet_name="all")
    el_demand = pd.read_excel(el_demand_fn, sheet_name="El_demand_MWh")

    heat_profile = pd.read_excel(heat_profile_fn, sheet_name=0)
    heat_demand = pd.read_excel(heat_demand_fn, sheet_name="DH_demand_MWh")

    h2_profile = pd.read_excel(h2_profile_fn, sheet_name=0)
    h2_demand = pd.read_excel(h2_demand_fn, sheet_name="H2_demand_MWh")

    kommun_nodes, el_kommun_demand, el_nodal_demand = group_demand(
                                                    subs_df,
                                                    vgr_kommun,
                                                    el_profile,
                                                    el_demand
                                                    )

    kommun_nodes, heat_kommun_demand, heat_nodal_demand = group_demand(
                                                    subs_df,
                                                    vgr_kommun,
                                                    heat_profile,
                                                    heat_demand
                                                    )

    kommun_nodes, h2_kommun_demand, h2_nodal_demand = group_demand(
                                                    subs_df,
                                                    vgr_kommun,
                                                    h2_profile,
                                                    h2_demand
                                                    )

    savefile(kommun_nodes, "kommun_newnode", model_input)
    savefile(el_kommun_demand, "el_kommun_demand", model_input)
    savefile(el_nodal_demand, "el_nodal_demand", model_input)
    savefile(heat_kommun_demand, "heat_kommun_demand", model_input)
    savefile(heat_nodal_demand, "heat_nodal_demand", model_input)
    savefile(h2_kommun_demand, "h2_kommun_demand", model_input)
    savefile(h2_nodal_demand, "h2_nodal_demand", model_input)
