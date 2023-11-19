import sys 
sys.path.append('/Users/7xw/Documents/work/ELM_ECP/lib/python3.8/site-packages/')
import numpy as np
import netCDF4 as nc

vis = 1 
save_file =0

# Need to change to binary file for large areas
AOI_gridcell_file = 'daymet_AK_all2.csv'
AOIdata = np.loadtxt(AOI_gridcell_file, delimiter = ",")

# Load the Daymet landmask for NA
mask_file = 'Daymet_NA_landmask.csv'
mask_data = np.loadtxt(mask_file, delimiter = ",")

x_file = 'x_space.csv'
x_data = np.loadtxt(x_file, delimiter = ",")
y_file = 'y_space.csv'
y_data = np.loadtxt(y_file, delimiter = ",")

total_rows = len(y_data)
total_cols = len(x_data)
mask_data = mask_data.reshape(total_rows, total_cols)

totoal_pts = len(AOIdata)
i=0
for mydata in AOIdata:
    x_idx = np.where(x_data == mydata[2])
    y_idx = np.where(y_data == mydata[3])
    #coordinates[i][0]= x_idx[0][0]
    #coordinates[i][1]= y_idx[0][0] 
    #print(coordinates[i][0], coordinates[i][1], i)
    mask_data[y_idx[0][0]][x_idx[0][0]] = 30
    i+=1

if vis==1: 
    fig, ax = plt.subplots()
    X,Y = np.meshgrid(x_data, y_data)
    CS = ax.contour(X[1000:2000,1000:2000], Y[1000:2000,1000:2000], mask_data[1000:2000,1000:2000])
    plt.show()
    
if save_file == 1:
    AOI_name = 'AK'
    AOI_mask_file = 'AOI_mask_'+ AOI_name +'.nc'

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(AOI_mask_file, 'w', format='NETCDF4')
    w_nc_fid.title = 'The Area of Interest for ELM simulation: '+ AOI_name 

    # create the gridIDs variable
    x_dim = w_nc_fid.createDimension('x_dim', total_cols)
    y_dim = w_nc_fid.createDimension('y_dim', total_rows)
    w_nc_var = w_nc_fid.createVariable('Mask', np.int32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'AOI_landmask: 0_Ocean, 1_Land, other_AOI'
    w_nc_fid.variables['Mask'][:] = mask_data
    w_nc_fid.close()  # close the new file
    
