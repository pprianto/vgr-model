import os
import requests
import pandas as pd
import json
import time
import warnings

# Suppress FutureWarning messages
warnings.simplefilter(action="ignore", category=FutureWarning)

"""
The following is used
to fetch PV and Wind hourly generation profile 
from Renewables Ninja API
https://www.renewables.ninja/
https://doi.org/10.1016/j.energy.2016.08.060
https://doi.org/10.1016/j.energy.2016.08.068
"""


def api_params(
    year,
    tech,
    api_base,
    axes_type=None,
    wind_type=None
):
    '''
    define the API parameters for PV and Wind
    '''

    # arguments depending on tech type
    if tech == "pv":

        url = api_base + "data/pv"

        com_args = {
            # "lat": 34.125,
            # "lon": 39.814,
            "local_time": True,
            "date_from": f"{year}-01-01",
            "date_to": f"{year}-12-31",
            "dataset": "merra2",
            "capacity": 1.0,
            "system_loss": 0.1,
            "tracking": axes_type,            # 0: fixed, 1: single axis, 2: double axis
            "tilt": 30,
            "azim": 180,
            "format": "json",
            "raw": True
        }

    else:

        url = api_base + "data/wind"

        if wind_type == "onshore":

            com_args = {
                # "lat": 34.125,
                # "lon": 39.814,
                "local_time": True,
                "date_from": f"{year}-01-01",
                "date_to": f"{year}-12-31",
                "dataset": "merra2",
                "capacity": 1.0,
                "height": 80,                  # assumed hub height
                "turbine": "Vestas V110 2000",   # turbine database, assumed
                "format": "json",
                "raw": True
            }
        
        else:

            com_args = {
                # "lat": 34.125,
                # "lon": 39.814,
                "local_time": True,
                "date_from": f"{year}-01-01",
                "date_to": f"{year}-12-31",
                "dataset": "merra2",
                "capacity": 1.0,
                "height": 150,                   # assumed hub height
                "turbine": "Vestas V164 9500",   # turbine database, assumed
                "format": "json",
                "raw": True
            }            

    return url, com_args


def fetch_loc(
    url,
    lat, 
    lon, 
    params
):
    '''
    generic API request
    '''
    
    params = params.copy()
    params.update({
                    "lat": lat,
                    "lon": lon
    })
    
    # Send token header with each request
    s = requests.session()
    s.headers = {"Authorization": "Token " + token}
    r = s.get(url, params=params)

    if r.status_code == 200:
        parsed_response = json.loads(r.text)
        data = pd.read_json(json.dumps(parsed_response["data"]), orient="index")

        return data

    else:
        print(f"Error fetching data for lat={lat}, lon={lon}:", r.status_code, r.text)
        return None


def fetch_profile(
    tech,
    year,
    url,
    params,
    subs_df,
    hourly_req = 50,            # 50 requests per hour
    time_interval = 3600 / 50   # fetch data interval so that the limit not exceeded
):
    '''
    fetch the profile in multiple locations of substations
    '''

    # Iterate over the DataFrame and fetch data
    req_count = 0
    start_time = time.time()

    nodal_profile_df = pd.DataFrame()

    for _, row in subs_df.iterrows():
        lat = row["lat"]
        lon = row["lon"]
        node = row["node_id"]
        df = fetch_loc(url, lat, lon, params)

        nodal_profile_df["hour"] = range(1, len(df) + 1)

        if df is not None:

            # electricity generation profile on each node        
            if len(nodal_profile_df) == len(df):
                nodal_profile_df[node] = df["electricity"].values

            else:
                print(f"Data length mismatch for {node}. Expected {len(df)}, got {len(nodal_profile_df)}.")


            if tech == "pv":
                
                if params["tracking"] == 0:
                    tilt = "fixed"
                elif params["tracking"] == 1:
                    tilt = "single_axis"
                else:
                    tilt = "double_axis"

                fn = f"{tech}_{tilt}_{year}_{node}.csv"
                
            else:
                
                if params["height"] == 80:
                    turb_type = "onshore"
                else:
                    turb_type = "offshore"                
                
                
                fn = f"{turb_type}_{tech}_{year}_{node}.csv"
            
            filename = os.path.join(profile_folder, fn)
            df.to_csv(filename, sep=";", index=False)
            print(f"{year} {tech} profile for {node} (lat={lat}, lon={lon}) saved to {filename}")
        
        # increment request
        req_count += 1

        # Check if reached the request limit
        if req_count >= hourly_req:
            elapsed_time = time.time() - start_time
            if elapsed_time < 3600:
                sleep_time = 3600 - elapsed_time
                print(f"Request limit reached. Sleeping for {sleep_time:.2f} seconds.")
                time.sleep(sleep_time)
            # Reset counters
            req_count = 0
            start_time = time.time()
        else:
            # Sleep between requests to respect rate limit
            time.sleep(time_interval)
            # time.sleep(1)             # for trial only


    if tech == "pv":
        profile_fn = f"nodal_profile_{tech}_{tilt}_{year}.csv"
                
    else:
        profile_fn = f"nodal_profile_{turb_type}_{tech}_{year}.csv"

    profile_filename = os.path.join(input_folder, profile_fn)
    nodal_profile_df.to_csv(profile_filename, sep=";", index=False)
    print(f"{tech} nodal profile for saved to {profile_filename}")

    return nodal_profile_df


def fetch_sequence(
    api,
    subs_df,
    techs,
    pv_types,
    wt_types,
    years,
):

    for tech in techs:

        print("=========================")
        print(f"Fetching {tech} profile")
        print("=========================")

        if tech == "pv":

            for pv_type in pv_types:
                
                if pv_type == 0:
                    axes = "fixed"
                else:
                    axes = "double_axis"


                print("=========================")
                print(f"Fetching {axes} profile")
                print("=========================")

                for year in years:

                    print("=========================")
                    print(f"Fetching {year} profile")
                    print("=========================")

                    url, com_args = api_params(year, tech, api, axes_type=pv_type)
                    nodal_profile = fetch_profile(tech, year, url, com_args, subs_df)

        else:

            for wt_type in wt_types:
                
                print("=========================")
                print(f"Fetching {wt_type} profile")
                print("=========================")

                for year in years:

                    print("=========================")
                    print(f"Fetching {year} profile")
                    print("=========================")

                    url, com_args = api_params(year, tech, api, wind_type=wt_type)
                    nodal_profile = fetch_profile(tech, year, url, com_args, subs_df)

        print("=========================")
        print(f"Fetching {tech} profile finished")
        print("=========================")
        
    print("=========================")
    print("Fetching completed")
    print("=========================")


if __name__ == "__main__":

    start_ = time.time()

    # df of substation for profiles
    current_dir = os.getcwd()
    input_folder = os.path.join(current_dir, "modelinput")
    mode = 0o666
    os.makedirs("ninjaprofile", mode, exist_ok=True)
    profile_folder = os.path.join(current_dir, "ninjaprofile")

    subs_fn = os.path.join(input_folder, "subs_clean.csv")
    subs = pd.read_csv(subs_fn, sep=";")

    # subs = subs[:2] # dummy for trial


    # API token from profile
    token = "insert API token here"
    api_base = "https://www.renewables.ninja/api/"


    # fetching data in loop for pv and wind
    re_techs = ["pv", "wind"]
    pv_axes = [0, 2]                       # fixed and double axis for opt tilt
    wt_techs = ["onshore", "offshore"]
    profile_years = [2019]#, 2020]#, 2021, 2022, 2023]  # fetched profile year


    fetch_sequence(
        api_base,
        subs,
        re_techs,
        pv_axes,
        wt_techs,
        profile_years
    )


    end_ = time.time()

    print("=========================")
    print(f"Entire fetching process completed in {(end_ - start_)/3600: .2f} hours")
    print("=========================")
