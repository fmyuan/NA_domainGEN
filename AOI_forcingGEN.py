
# data_partition module for batch processing
# based on array_split and function definition

import os,sys
import glob
import math
import netCDF4 as nc
import numpy as np
import pandas as pd
from time import process_time
from datetime import datetime

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m-%d-%Y')

def AOI_forcing_save_1d(input_path, file, AOI, AOI_points, output_path):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    src = nc.Dataset(source_file, 'r', format='NETCDF4')

    #read gridIDs
    grid_ids = src['gridID'][...]    # gridID for all NA
 
    #  
    AOI_idx = np.where(np.in1d(grid_ids, AOI_points))[0]
    
    #
    dst_name = output_path + '/'+ AOI + '/'+file
    # check if file exists then delete it
    if not os.path.exists(output_path+'/'+AOI): 
        os.makedirs(output_path+'/'+AOI)
    elif os.path.exists(dst_name):
        os.remove(dst_name)

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF4')
    dst.title = './'+ AOI + '/'+file+' creted from '+ source_file +' on ' +formatted_date

    # Copy the global attributes from the source to the target
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    # Copy the dimensions from the source to the target
    for name, dimension in src.dimensions.items():
        if name != 'ni' and name != 'gridcell':
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        else:
            # Update the 'ni' dimension with the length of the list
            #dst.dimensions['ni'].set_length(len(AOI_points))
            if name == 'ni' or name == 'gridcell': 
                ni = dst.createDimension(name, AOI_points.size)

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)   
        print(name, variable.dimensions)
        
        if (name != 'lambert_conformal_conic'):
            if variable.dimensions[-1] != 'ni' and variable.dimensions[-1] != 'gridcell':
                dst[name][...] = src[name][...]
            else:
                dst[name][...] = src[name][..., AOI_idx]
           
        # Copy the variable attributes
        for attr_name in variable.ncattrs():
            dst[name].setncattr(attr_name, variable.getncattr(attr_name))
        
    src.close()  # close the source file 
    dst.close()  # close the new file        
    
def get_files(input_path, ncheader='clmforc'):
    print(input_path)
    #files = os.listdir(input_path)
    files = glob.glob("%s*.%s" % (input_path+ncheader,'nc'))
    
    files.sort() 

    #files_nc = [f for f in files if (f[-2:] == 'nc')] 
    print("total " + str(len(files)) + " files need to be processed")

    return files

def main():
    args = sys.argv[1:]

    if len(sys.argv) != 4  or sys.argv[1] == '--help':  # sys.argv includes the script name as the first argument
        print("Example use: python AOI_forcingGEN.py <input_path> <output_path> <AOI_points_file>")
        print(" <input_path>: path to the 1D source data directory")
        print(" <output_path>:  path for the 1D AOI forcing data directory")
        print(" <AOI_points_file>:  <AOI>_gridID.csv or <AOI>_domain.nc")
        print(" The code uses NA forcing to generation 1D AOI forcing")              
        exit(0)

    input_path = args[0]
    if not input_path.endswith("/"): input_path=input_path+'/'
    output_path = args[1]
    if not output_path.endswith("/"): output_path=output_path+'/'
    AOI_gridID_file = args[2]
    AOI=AOI_gridID_file.split("_")[0]

    if (AOI_gridID_file.endswith('gridID.csv')):
        #AOI_gridcell_file = AOI+'_gridID.csv'  # user provided gridcell IDs
        df = pd.read_csv(AOI_gridID_file, sep=",", skiprows=1, names = ['gridID'])
        #read gridIds
        AOI_points = np.array(df['gridID'])
    elif AOI_gridID_file.endswith('domain.nc'):
        src = nc.Dataset(AOI_gridID_file, 'r')
        AOI_points = src['gridID'][:]
    else:
        print("Error: Invalid AOI_points_file, see help.")
        
    files_nc = get_files(input_path)
    n_files = len(files_nc)
    #
#------------------------------------------------------------------------------------------
# mpi implementation - simply round-robin 'n_files' over cpu_cores
    if HAS_MPI4PY:
        mycomm = MPI.COMM_WORLD
        myrank = mycomm.Get_rank()
        mysize = mycomm.Get_size()
    else:
        mycomm = 0
        myrank = 0
        mysize = 1

    ng = math.floor(n_files/mysize)
    ng_rank = np.full([mysize], int(1))
    ng_rank = np.cumsum(ng_rank)*ng
    xg = int(math.fmod(n_files, mysize))
    xg_rank = np.full([mysize], int(0))
    if xg>0: xg_rank[:xg]=1
    ng_rank = ng_rank + np.cumsum(xg_rank) - 1        # ending file index, starting 0, for each rank
    ng0_rank = np.hstack((0, ng_rank[0:mysize-1]+1))  # starting file index, starting 0, for each rank

    for i in range(ng0_rank[myrank], ng_rank[myrank]+1):
        f = files_nc[i]

        if (not f.split('/')[-1].startswith('clmforc')): continue
        print('processing '+ f )
        start = process_time() 
        AOI_forcing_save_1d(input_path, f.split('/')[-1], AOI, AOI_points, output_path)
        end = process_time()
        print("Generating 1D forcing data for "+AOI+ " domain takes {}".format(end-start))

if __name__ == '__main__':
    main()
    