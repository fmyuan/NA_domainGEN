import netCDF4 as nc
import numpy as np
from pyproj import Transformer
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import pandas as pd
import sys

# An input csv file that contains the locations of land gridcells in Daymet domain;
# 1:  land gridcell ID
# 2:  lat, lon of land gridcell center 
# 3:  Y,X location in LLC projection

# output 1D AOIdomain.nc that can be used to generate 1D AOIsurfdata.nc and 1D AOIforcing.nc 

def read_gridcells(file_name):
    try:
        df = pd.read_csv(file_name, sep=",")
    except pd.errors.ParserError:
        df = pd.read_csv(file_name, sep=" ")
    return df

def find_nearest_points(listA, listB):
    tree = cKDTree(np.array(listB))
    _, indices = tree.query(listA, k=1)
    return indices

def main():
    args = sys.argv[1:]
    AOI_gridcell_file = args[0]  # user provided gridcell csv file
    user_option = args[1]
    AOI = args[2]

    # save to the 1D domain file
    AOIdomain = str(AOI)+'domain.nc'
    dst = nc.Dataset(AOIdomain, 'w', format='NETCDF4')

    # open the 1D domain data
    r_domain = nc.Dataset('Daymet4.1km.1d.domain.nc', 'r', format='NETCDF4')

    # 3) Open a csv file to read a list of points (y, x)
    #df = read_gridcells(AOI_gridcell_file)  

    if user_option == '1': # gridID is used directly

        df = pd.read_csv(AOI_gridcell_file, sep=" ", skiprows=1, names = ['gridIDs'])
        #read gridIds
        AOI_points = list(df['gridIDs'])

        # read gridIDs
        NA_gridIds = r_domain.variables['gridIDs'][:]
        NA_gridcell_list = list(NA_gridIDs)
        print(NA_gridcell_list[0:5])

    if user_option == '2': # use lat lon coordinates
        df = pd.read_csv(AOI_gridcell_file, sep=" ", skiprows=1, names = ['yc', 'xc', 'yc_LLC', 'xc_LLC'])
        # read yc, xc (lat, lon)
        NA_yc = r_domain.variables['yc'][:]
        NA_xc = r_domain.variables['xc'][:]

        #read in y, x coordinate (lat, lon)
        AOI_points = list(zip(df['yc'], df['xc']))
        # create list for all land gridcell (lat, lon)
        NA_gridcell_list = list(zip(NA_yc, NA_xc))
        print(NA_gridcell_list[0:5])

    if user_option == '3': # gridID is used directly

        df = pd.read_csv(AOI_gridcell_file, sep=" ", skiprows=1, names = ['yc', 'xc', 'yc_LLC', 'xc_LLC'])
        # read yc_llc, xc_llc (y, x in LLC)
        NA_yc_LLC = r_domain.variables['yc_LLC'][:]
        NA_xc_LLC = r_domain.variables['xc_LLC'][:]

        #read in y, x coordinate (in LLC projection)
        AOI_points = list(zip(df['yc_LLC'], df['xc_LLC']))
        # create list for all land gridcell (lat, lon)
        NA_gridcell_list = list(zip(NA_yc_LLC, NA_xc_LLC))
        print(NA_gridcell_list[0:5])
    
    #Find the nearest points in the 1D list of points and return the index of these nearest points
    NA_domain_idx = find_nearest_points(AOI_points, NA_gridcell_list)
    print(NA_domain_idx[0:10])

    # Copy the global attributes from the source to the target
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    # Copy the dimensions from the source to the target
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]

        # Copy the variable attributes
        for attr_name in variable.ncattrs():
            dst[name].setncattr(attr_name, variable.getncattr(attr_name))

    # Update the 'ni' dimension with the length of the list
    dst.dimensions['ni'].set_length(len(AOI_points))

    # Sort the points
    AOI_points.sort()

    # Use the points as index to copy all the elements in each variable
    for name in src.variables:
        dst[name][:] = src[name][AOI_points]
    
    # Close the source netCDF file
    r_domain.close()

    # Save the target netCDF file
    dst.close()

if __name__ == '__main__':
    main()