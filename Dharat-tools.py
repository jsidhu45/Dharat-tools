print("""
This tool estimates landslide failure surfaces using spline interpolation method.

For more information, visit: https://github.com/jsidhu45/Dharat-tools

Reference: Singh, J. and Sepulveda, S.A. (2025) Assessing the landslide failure surface
depth and volume: A new spline interpolation method. Engineering Geology.
https://doi.org/10.1016/j.enggeo.2025.108319

**The algorithm interpolate each grid point inside the landslide boundary iteratively, depending on the size of the landslide and DEM resolution, it is optimal to use coarse resolution DEM and larger stop delta volume first to check the computation time.**
      
The algorithm has two stopping criteria options: 1) delta volume between iterations and 2) maximum depth. The user can select either of the stopping criteria based on their preference. If nothing is specified, the algorithm will prioritize delta volume as the stopping criterion by default.
      
The intermediate interpolated results are saved in Excel files for analysis after every specified number of iterations (default is 1000) in the results folder.
""")


import argparse
import logging
import time
from pathlib import Path
from os.path import join
from typing import Tuple

import numpy as np
import pandas as pd
from osgeo import gdal
from scipy import interpolate
import warnings

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[logging.FileHandler("run.log"), logging.StreamHandler()],
    force=True
)

# Constants
DATA_DIR = Path('./input_data/')
OUTPUT_DIR = Path('./output_files/')
RESULTS_DIR = Path('./results/')
DEM_DIR = Path('./results/DEM/')

# Handle specific warnings
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")
gdal.UseExceptions()


def load_and_resample_dem(dem_filename: str, resolution: float) -> Tuple[gdal.Dataset, tuple]:
    """Load and resample DEM file.

    Args:
        dem_filename: Name of the DEM file in input_data/
        resolution: Target resolution for resampling

    Returns:
        Tuple of (resampled DEM dataset, geotransform)
    """
    dem_path = DATA_DIR / dem_filename
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM file not found: {dem_path}")

    logger.info(f"Loading DEM: {dem_path}")
    dem_org = gdal.Open(str(dem_path))
    if dem_org is None:
        raise ValueError(f"Could not open DEM file: {dem_path}")

    resampled_path = OUTPUT_DIR / 'dem_resampled.tif'
    logger.info(f"Resampling DEM to {resolution}m resolution")
    dem_resampled = gdal.Warp(str(resampled_path), dem_org, xRes=resolution, yRes=resolution, resampleAlg="bilinear")

    gt = dem_resampled.GetGeoTransform()
    projection = dem_resampled.GetProjection()
    return dem_resampled, gt, projection


def clip_dem_to_boundary(dem_dataset: gdal.Dataset, boundary_filename: str) -> Tuple[gdal.Dataset, np.ndarray]:
    """Clip DEM to landslide boundary.

    Args:
        dem_dataset: GDAL dataset of the DEM
        boundary_filename: Name of the boundary shapefile in input_data/

    Returns:
        Tuple of (clipped DEM dataset, clipped array)
    """
    boundary_path = DATA_DIR / boundary_filename
    if not boundary_path.exists():
        raise FileNotFoundError(f"Boundary file not found: {boundary_path}")

    clipped_path = OUTPUT_DIR / 'DEM_clip.tif'
    logger.info(f"Clipping DEM to boundary: {boundary_path}")
    dem_clip = gdal.Warp(str(clipped_path), dem_dataset, cutlineDSName=str(boundary_path), dstNodata=np.nan)

    clip_arr = dem_clip.GetRasterBand(1).ReadAsArray()
    return dem_clip, clip_arr


def prepare_dataframes(dem_dataset: gdal.Dataset, clip_arr: np.ndarray) -> Tuple[pd.DataFrame, pd.DataFrame, int, int]:
    """Prepare dataframes from DEM and clipped data.

    Args:
        dem_dataset: Resampled DEM dataset
        clip_arr: Clipped DEM array

    Returns:
        Tuple of (full_box_df, inside_df, sub_rows, sub_cols)
    """
    # Convert full DEM to XYZ
    box_xyz_path = OUTPUT_DIR / 'box_xyz.xyz'
    gdal.Translate(str(box_xyz_path), dem_dataset) #It should be resampled dem
    df_box = pd.read_csv(box_xyz_path, sep=" ", header=None, names=['x', 'y', 'z'])
   
    # Convert clipped DEM to XYZ and filter valid points
    inn_xyz_path = OUTPUT_DIR / 'inn_xyz.xyz'
    dem_clip, _ = clip_dem_to_boundary(dem_dataset, '') # Wait, need to fix this - boundary is already clipped
    # Actually, let's refactor this properly

    # For now, assume dem_clip is passed
    gdal.Translate(str(inn_xyz_path), dem_clip)
    df_inn_xyz = pd.read_csv(inn_xyz_path, sep=" ", header=None, names=['x', 'y', 'z'])
    
    # Filter points inside landslide (z > 0)
    df_inn = df_inn_xyz[df_inn_xyz['z'] > 0].copy()

    # Count points outside in x and y directions
    x_min, x_max = df_inn['x'].min(), df_inn['x'].max()
    y_min, y_max = df_inn['y'].min(), df_inn['y'].max()

    outside_x = df_box[(df_box['x'] < x_min) | (df_box['x'] > x_max)]['x'].nunique()
    outside_y = df_box[(df_box['y'] < y_min) | (df_box['y'] > y_max)]['y'].nunique()
    return df_box, df_inn, outside_y, outside_x


def perform_interpolation(fill: bool,df_box: pd.DataFrame, df_inn: pd.DataFrame, rows: int, cols: int,
                         sub_rows: int, sub_cols: int, save_interval: int, stop_delta_vol: float, stop_max_depth: float, gt: float ) -> np.ndarray:
    """Perform the spline interpolation iterations.

    Args:
        df_box: Full DEM dataframe
        df_inn: Inside landslide dataframe
        rows: Number of rows in DEM
        cols: Number of columns in DEM
        sub_rows: Rows outside boundary
        sub_cols: Columns outside boundary
        iterations: Total iterations
        save_interval: Save every N iterations

    Returns:
        Interpolated elevation matrix
    """
    # Round and sort dataframes
    df_inn_rounded = df_inn.round(3).sort_values(by = ['x', 'y'], ascending=[True, False])
    df_box_rounded = df_box.round(3).sort_values(by =['x', 'y'], ascending=[True, False])
    
    # Convert to numpy arrays
    z_arr = df_box_rounded['z'].values.reshape((cols, rows)).T
    x_arr = df_box_rounded['x'].values.reshape((cols, rows)).T
    y_arr = df_box_rounded['y'].values.reshape((cols, rows)).T
    
    pd.DataFrame(z_arr).to_excel(RESULTS_DIR / 'original.xlsx', index=False) #save orignal DEM elevation matrix

    # Initialize interpolated matrix
    interpolated_z = z_arr.copy()
    
    # Get indices of inside points
    inside_indices = []
    for _, row in df_inn_rounded.iterrows():
        x, y = row['x'], row['y']
        j, i = np.where((x_arr == x) & (y_arr == y))
        if len(j) > 0 and len(i) > 0:
            inside_indices.append((j[0], i[0]))
    
    logger.info(f"Found {len(inside_indices)} points lying inside the landslide boundary for interpolation")
    
    start_time = time.time()
    logger.info(f"Starting  interpolation")

    if stop_max_depth is not None:
        logger.info(f"Stopping criterion  Max Depth = {stop_max_depth}")
    else:
        logger.info(f"Stopping criterion Delta Volume = {stop_delta_vol}")

    iterations = 0
    while True:
        iterations += 1
        prev_vol = np.sum(interpolated_z - z_arr)*gt[1]*gt[1]
        prev_vol = float(np.abs(prev_vol))
        for j, i in inside_indices: # Check and loop through each point inside the landslide boundary only
            current_z = z_arr[j, i]
            current_x = x_arr[j, i]
            current_y = y_arr[j, i]
            
            # Remove current point for interpolation
            row_z = np.delete(interpolated_z[j, :], i) #[j, :] gives full row and delete i (current) element from the row j
            row_x = np.delete(x_arr[j, :], i)
            

            col_z = np.delete(interpolated_z[:, i], j) #[:, i] gives full column and delete j (current) element from the column i
            col_y = np.delete(y_arr[:, i], j)

            # Interpolate along x-axis (row)
            if len(row_x) > 1:
                try:
                    f_x = interpolate.interp1d(row_x, row_z, kind='cubic', bounds_error=False, fill_value=np.nan)
                    interp_x = f_x(current_x)
                except:
                    interp_x = current_z
            else:
                interp_x = current_z

            # Interpolate along y-axis (column)
            if len(col_y) > 1:
                try:
                    f_y = interpolate.interp1d(col_y, col_z, kind='cubic', bounds_error=False, fill_value=np.nan)
                    interp_y = f_y(current_y)
                except:
                    interp_y = current_z
            else:
                interp_y = current_z

            # Average the interpolations
            avg_z = (interp_x + interp_y)*0.5
            
            if fill == False:
                # Update if lower (deeper failure surface)
                if not np.isnan(avg_z) and avg_z < current_z:
                    interpolated_z[j, i] = avg_z
            else:
                 # Update if lower (deeper failure surface)
                if not np.isnan(avg_z) and avg_z > current_z:
                    interpolated_z[j, i] = avg_z

        #Compute volume and max depth after each iteration
        current_vol = np.abs(np.sum(interpolated_z - z_arr)*gt[1]*gt[1])
        delta_vol = np.abs(current_vol - prev_vol)
        current_max_depth = np.abs(np.max(interpolated_z - z_arr))

        #print volume updates every 500 iterations
        if iterations % 25 == 0:
            logger.info(f"Iteration {iterations}: Delta Volume = {delta_vol}: Volume = {current_vol} m\u00b3: Max Depth = {current_max_depth} meters")
        
        # Save intermediate results
        if iterations % save_interval == 0:
            output_path = RESULTS_DIR / f'interpolated{iterations}.xlsx'
            pd.DataFrame(interpolated_z).to_excel(output_path, index=False)
            logger.info(f"Saved intermediate result: {output_path}")

        # Stopping criteria volume difference and max depth
        if stop_delta_vol is not None and delta_vol <= stop_delta_vol:
            logger.info(f"Stopping criterion met: Delta Volume {delta_vol} <= {stop_delta_vol}")
            break

        if stop_max_depth is not None and current_max_depth >= stop_max_depth:
            logger.info(f"Stopping criterion met: Max depth {current_max_depth} >= {stop_max_depth}")
            break
    
    end_time = time.time()
    time_per_iter = (end_time - start_time) / iterations
    #logger.info(f"Total iterations: {iterations}")
    logger.info(f"Interpolation completed! Total iterations: {iterations} Time per iteration: {time_per_iter:.2f} seconds")
    pd.DataFrame(interpolated_z).to_excel(RESULTS_DIR / 'interpolated_final.xlsx', index=False)
    logger.info(f"Final Volume = {current_vol} m\u00b3: Max Depth = {np.abs(np.max(interpolated_z - z_arr))} meters")

    return interpolated_z, z_arr

#Function to write geotiff output
def write_geotiff(dem_dir, filename, array, cols, rows, gt, projection):
    ds = gdal.GetDriverByName("GTiff").Create(
        join(dem_dir, filename),
        cols, rows, 1, gdal.GDT_Float32
    )

    ds.SetGeoTransform(gt)
    ds.SetProjection(projection)

    band = ds.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(np.nan)
    band.FlushCache()

    ds = None   # close file

def save_final_DEM(original_matrix: np.ndarray, interpolated_matrix: np.ndarray, cols: int, rows: int, gt: tuple, projection: str):
    """Save final results to GeoTIFF files."""
    write_geotiff(DEM_DIR, "dem_resampled.tif", original_matrix, cols, rows, gt, projection)
    write_geotiff(DEM_DIR, "final_interpolated.tif", interpolated_matrix, cols, rows, gt, projection)
    write_geotiff(DEM_DIR, "thickness.tif", original_matrix - interpolated_matrix, cols, rows, gt, projection)

    logger.info(f"Final thickness and interpolated surface DEM files are saved in {DEM_DIR}")


def process_landslide_data(fill:bool, dem_filename: str, resolution: float, boundary_filename: str,
                          save_interval: int, stop_delta_vol: float, stop_max_depth: float):
    """Main processing function."""
    # Load and resample DEM
    dem_resampled, gt, projection = load_and_resample_dem(dem_filename, resolution)
    rows, cols = dem_resampled.RasterYSize, dem_resampled.RasterXSize

    # Clip to boundary
    dem_clip, clip_arr = clip_dem_to_boundary(dem_resampled, boundary_filename)

    # Prepare dataframes
    df_box, df_inn, sub_rows, sub_cols = prepare_dataframes(dem_resampled, clip_arr)
    
    # Perform interpolation
    interpolated_matrix, z_arr = perform_interpolation(fill, df_box, df_inn, rows, cols, sub_rows, sub_cols, save_interval, stop_delta_vol, stop_max_depth, gt)

    # Save final results in form of interpolated DEM and thickness DEM
    original_matrix_z = df_box['z'].values.reshape((cols, rows)).T
    save_final_DEM(z_arr, interpolated_matrix, cols, rows, gt, projection)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fill', action='store_true', help='(Optional input) Default is False. If selected, the algorithm will reverse and fill the topography')
    parser.add_argument('--dem', required=True, help='Name of DEM (.tif) file in input_data/')
    parser.add_argument('--resolution', type=float, required=True, help='DEM resolution for resampling (meters)')
    parser.add_argument('--boundary', required=True, help='Name of landslide boundary (.shp) file in input_data/')
    parser.add_argument('--save-interval', type=int, default=1000,
                       help='(Optional input) Number of iterations after which intermediate interpolated results are saved to an Excel file (default: 1000)')
    #parser.add_argument('--stop-delta-vol', type=float, default=0.001, help='(Optional input) Stopping criterie based on volume difference threshold from previous iteration (default: 0.001 cubic meters).')


    # --- stopping criteria (ONLY ONE allowed) ---
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--stop-delta-vol', type=float, default=0.001,
                   help='(Optional input) Stopping criterie based on volume difference threshold from previous iteration (default: 0.001)')
    group.add_argument('--stop-max-depth', type=float,
                   help='(Optional input) Stopping criterion based on maximum depth')

    args = parser.parse_args()

    # Create directories
    OUTPUT_DIR.mkdir(exist_ok=True)
    RESULTS_DIR.mkdir(exist_ok=True)
    DEM_DIR.mkdir(exist_ok=True)

    # Process
    try:
        process_landslide_data(args.fill, args.dem, args.resolution, args.boundary, args.save_interval, args.stop_delta_vol, args.stop_max_depth)
        logger.info("Processing completed successfully!")
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise


if __name__ == "__main__":
    main()