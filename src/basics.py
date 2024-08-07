import os
import geopandas as gpd
import pandas as pd

def prerequisites(folder_name):

    current_dir = os.getcwd()
    # subfolder = f"{folder_name}" # "inputs"
    path_folder = os.path.join(current_dir, folder_name)
    mode = 0o666
    os.makedirs(path_folder, mode, exist_ok=True)

    return path_folder


def count_geometry_cols(df):
    geometry_cols = [col for col, dtype in df.dtypes.items() if dtype == 'geometry']
    no_geo_cols = len(geometry_cols)
    
    return no_geo_cols


def reg_bound(nat_df):
    """
    Define areas and municipality as boundary
    """

    # Filter the region of interest, 14 in the case of VGR
    # iloc[:, 0] because the 1st column gives the kommun or län kod
    reg_df = nat_df[nat_df.iloc[:, 0].astype(str).str.contains('^14')]
    
    # Convert crs format into lon and lat
    reg_df = reg_df.to_crs("EPSG:4326")

    return reg_df


def getinputfiles(fn, folder_name):

    current_dir = os.getcwd()
    # subfolder = "inputs"
    file_path = os.path.join(current_dir, folder_name, fn)

    return file_path


def savefile(df, filename, folder_name):
    """
    Saves as geojson for geopandas or xls for pandas
    """

    file_path = getinputfiles(filename, folder_name)

    no_geo_cols = count_geometry_cols(df)

    if no_geo_cols > 1:
        col_to_keep = 'geometry'
        geometry_cols = [col for col in df.columns if df[col].dtype == 'geometry' and col != col_to_keep]
        df.drop(geometry_cols, axis=1, inplace=True)


    if isinstance(df, gpd.GeoDataFrame):

        # df.to_file(path, f"{filename}.geojson", driver="GeoJSON")
        df.to_file(f"{file_path}.geojson", driver="GeoJSON")
        df_to_save = df.drop('geometry', axis=1)
        
        df_to_save.to_csv(f"{file_path}.csv", sep=';', index=False, encoding='utf-8-sig')
        df_to_save.to_excel(f"{file_path}.xlsx", sheet_name=f"{filename}", index=False)

    else:
        df.to_csv(f"{file_path}.csv", sep=';', index=False, encoding='utf-8-sig')
        df.to_excel(f"{file_path}.xlsx", sheet_name=f"{filename}", index=False)

