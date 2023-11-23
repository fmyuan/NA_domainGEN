import netCDF4 as nc
import numpy as np
from pyproj import Transformer
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import pandas as pd
import sys, os

from datetime import datetime

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m%d%Y')

# An input csv file that contains the locations of land gridcells in Daymet domain;
# 1:  land gridcell ID
# 2:  xc, yc (lon, lat) of land gridcell center 
# 3:  xc_LCC,yc_LCC (x,y) location in LLC projection

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

    # Check the number of arguments
    if len(sys.argv) != 4  or sys.argv[1] == '--help':  # sys.argv includes the script name as the first argument
        print("Example use: python AOIcreation.py <points_file_name> <AOI_domain_name>")
        print(" <points_file_name>:  <domain>_gridID.csv, <domain>_xcyc.csv, <domain>_xcyc_LCC.csv")
        print(" <user_option>:  1-> gridID, 2-> xcyc, 3-> xcyc_LCC")
        print(" <AOI_domain>:  Uniqe string for AOI domain name")
        exit(0)

    args = sys.argv[1:]
    AOI_gridcell_file = args[0]  # user provided gridcell csv file
    user_option = args[1]
    AOI = args[2]

    # save to the 1D domain file
    AOIdomain = str(AOI)+'domain.nc'

    # check if file exists then delete it
    if os.path.exists(AOIdomain):
        os.remove(AOIdomain)

    source_file = 'Daymet4.1km.1d.domain.nc'
    dst = nc.Dataset(AOIdomain, 'w', format='NETCDF4')

    # open the 1D domain data
    src = nc.Dataset(source_file, 'r', format='NETCDF4')

    # 3) Open a csv file to read a list of points (y, x)
    #df = read_gridcells(AOI_gridcell_file)  

    if user_option == '1': # gridID is used directly
        #AOI_gridcell_file = AOI+'_gridID.csv'  # user provided gridcell IDs
        df = pd.read_csv(AOI_gridcell_file, sep=",", skiprows=1, names = ['gridID'])
        #read gridIds
        AOI_points = list(df['gridID'])
        #AOI_points = list(mygridIDs)

        # read gridIDs
        NA_gridIDs = src.variables['gridID'][:]
        NA_gridcell_list = list(NA_gridIDs)
        print(NA_gridcell_list[0:5])

        domain_idx = np.where(np.in1d(NA_gridcell_list, AOI_points))[0]

    if user_option == '2': # use lat lon coordinates
        #AOI_gridcell_file = AOI+'_xcyc.csv'  # user provided gridcell csv file  (xc, yc) (lon, lat)
        df = pd.read_csv(AOI_gridcell_file, sep=",", skiprows=1, names = ['xc', 'yc'], engine='python')

        #read in x, y coordinate (lon, lat)
        AOI_points = list(zip(df['xc'], df['yc']))
        #AOI_points = list(zip(myxc, myyc))
        
        # read yc, xc (y, x)
        NA_yc = src.variables['yc'][:]
        NA_xc = src.variables['xc'][:]
        # create list for all land gridcell (lat, lon)
        NA_gridcell_list = list(zip(NA_xc, NA_yc))

        AOI_points_arr = np.array(AOI_points)
        print("AOI_points_arr", AOI_points_arr[0:5], "shape", AOI_points_arr.shape)
        NA_gridcell_arr = np.squeeze(np.array(NA_gridcell_list)).transpose()    
        print("NA_gridcell_arr", NA_gridcell_arr[0:5], "shape", NA_gridcell_arr.shape)

        tree = cKDTree(NA_gridcell_arr)
        _, domain_idx = tree.query(AOI_points_arr, k=1)

    if user_option == '3': # xc_LCC and yc_LCC is used directly
        #AOI_gridcell_file = AOI+'_XYLCC.csv'  # user provided gridcell csv file  (xc_LLC, yc_LLC) 
        df = pd.read_csv(AOI_gridcell_file, sep=",", skiprows=1, names = ['xc_LCC', 'yc_LCC'], engine='python')

        #read in x, y coordinate (in LLC projection)
        AOI_points = list(zip(df['xc_LCC'], df['yc_LCC']))
        #AOI_points = list(zip(myxc_lcc, myyc_lcc))

        # read yc_LCC, xc_LCC (y, x in LCC)
        NA_yc_LCC = src.variables['yc_LCC'][:]
        NA_xc_LCC = src.variables['xc_LCC'][:]

        # create list for all land gridcell (lat, lon)
        NA_gridcell_list = list(zip(NA_xc_LCC, NA_yc_LCC))

        AOI_points_arr = np.array(AOI_points)
        print("AOI_points_arr", AOI_points_arr[0:5], "shape", AOI_points_arr.shape)
        NA_gridcell_arr = np.squeeze(np.array(NA_gridcell_list)).transpose()    
        print("NA_gridcell_arr", NA_gridcell_arr[0:5], "shape", NA_gridcell_arr.shape)

        tree = cKDTree(NA_gridcell_arr)
        _, domain_idx = tree.query(AOI_points_arr, k=1)

    domain_idx = np.sort(domain_idx).squeeze()
    print("gridID_idx", domain_idx[0:10])
    #np.savetxt("AOIgridId_idx.csv", domain_idx[0:100], delimiter=",", fmt='%d')

    # Copy the global attributes from the source to the target
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    # Copy the dimensions from the source to the target
    for name, dimension in src.dimensions.items():
        if name != 'ni':
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        else:
            # Update the 'ni' dimension with the length of the list
            #dst.dimensions['ni'].set_length(len(AOI_points))
            ni = dst.createDimension("ni", len(AOI_points))

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)   
        print(name, variable.dimensions)
        
        if (name != 'lambert_conformal_conic'):
            if (variable.dimensions[-1] != 'ni'):
                dst[name][:] = src[name][:]
            else:
                dst[name][:] = src[name][0][domain_idx]
           
        # Copy the variable attributes
        for attr_name in variable.ncattrs():
            dst[name].setncattr(attr_name, variable.getncattr(attr_name))

    dst.title = '1D domain for '+ AOI +', generated on ' +formatted_date + ' with ' + source_file
       
    # Close the source netCDF file
    src.close()

    # Save the target netCDF file
    dst.close()

if __name__ == '__main__':
    main()