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
        print("Example use: python AOI_surfdataGEN.py <input_path> <output_path> <AOI_points_file>")
        print(" <input_path>: path to the 1D source data directory")
        print(" <output_path>:  path for the 1D AOI surface data directory")
        print(" <AOI_points_file>:  <AOI>_gridID.csv or <AOI>_domain.nc")
        print(" The code uses NA forcing to generation 1D AOI forcing")      
        exit(0)

    input_path = args[0]
    output_path = args[1]
    AOI_gridID_file = args[2]
    AOI=AOI_gridID_file.split("_")[0]

    if (AOI_gridID_file.endswith(gridID.csv)):
        #AOI_gridcell_file = AOI+'_gridID.csv'  # user provided gridcell IDs
        df = pd.read_csv(AOI_gridID_file, sep=",", skiprows=1, names = ['gridID'])
        #read gridIds
        AOI_points = np.array(df['gridID'])
    elif filename.endswith('domain.nc'):
        src = nc.Dataset(AOI_gridID_file, 'r')
        AOI_points = src['gridID'][:]
    else:
        print("Error: Invalid AOI_points_file, see help.")

    # save to the 1D domain file
    AOIsurfdata = str(AOI)+'surfdata.nc'

    # check if file exists then delete it
    if os.path.exists(AOIsurfdata):
        os.remove(AOIsurfdata)

    source_file = 'Daymet4.1km.1d.surfdata.nc'
    dst = nc.Dataset(AOIsurfdata, 'w', format='NETCDF4')

    # open the 1D domain data
    src = nc.Dataset(source_file, 'r', format='NETCDF4')

    # read gridIDs from src file
    NA_gridIDs = src.variables['gridID'][:]
    NA_gridcell_list = list(NA_gridIDs)
    print(NA_gridcell_list[0:5])

    # get the index of AOI_points in the NA_gridcell_list
    domain_idx = np.where(np.in1d(NA_gridcell_list, AOI_points))[0]

    # domain_idx = np.sort(domain_idx).squeeze()
    print("gridID_idx", domain_idx[0:10])

    # Copy the global attributes from the source to the target
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    # Copy the dimensions from the source to the target
    for name, dimension in src.dimensions.items():
        if name != 'gridcell':
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        else:
            # Update the 'ni' dimension with the length of the list
            #dst.dimensions['ni'].set_length(len(AOI_points))
            ni = dst.createDimension("gridcell", len(AOI_points))

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)   
        print(name, variable.dimensions)
        
        if (variable.dimensions[-1] != 'gridcell'):
            dst[name][:] = src[name][:]
        else:
            dst[name][:] = src[name][0][domain_idx]
           
        # Copy the variable attributes
        for attr_name in variable.ncattrs():
            dst[name].setncattr(attr_name, variable.getncattr(attr_name))

    dst.title = '1D surfdata for '+ AOI +', generated on ' +formatted_date + ' with ' + source_file
       
    # Close the source netCDF file
    src.close()

    # Save the target netCDF file
    dst.close()

if __name__ == '__main__':
    main()