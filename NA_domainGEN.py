
# data_partition module for batch processing
# based on array_split and function definition

import os 
import netCDF4 as nc
import numpy as np
from itertools import cycle
from time import process_time

def data_read(file_name, var_name, timesteps):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    r_nc_fid = nc.Dataset(file_name, 'r', format='NETCDF4')

    total_rows = r_nc_fid.dimensions['x'].size
    total_cols = r_nc_fid.dimensions['y'].size
    total_timesteps = r_nc_fid.dimensions['time'].size
    x_dim = r_nc_fid['x']
    y_dim = r_nc_fid['y']
    XC, YC = np.meshgrid(x_dim, y_dim)  # the array is (y,x) to match the mask
    #YX = np.array([Y, X])
    #print('total timesteps is :' + str(total_timesteps))
    if timesteps == -1:
        timesteps = total_timesteps
    data = r_nc_fid[var_name][0:timesteps, :, :] # read (timestep, y, x) format
    lon = r_nc_fid['lon'][:,:]  # read lon(y, x) format
    lat = r_nc_fid['lat'][:,:]  # read lat(y, x) format

    print(lon.shape, lat.shape, XC.shape, YC.shape) 
    return total_rows, total_cols, timesteps, data, lon, lat, XC, YC

def domain_info_1D(total_rows, total_cols, timesteps, data, lon, lat, XC, YC):
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(total_cols,total_rows)

    """
    # create a mask for land grid_ids (1)
    mask = data[0]    # FSDS is in (time, Y, X) format
    mask = np.where(~np.isnan(mask), 1, 0)
    # Create a boolean mask where True indicates non-NaN values
    bool_mask = ~np.isnan(data[0])

    # create an flattened list of land gridID and reduce the size of gridIDs array
    grid_ids = grid_ids.reshape(total_cols,total_rows)
    grid_ids = np.multiply(mask,grid_ids)
    grid_ids = grid_ids[grid_ids != 0]

    # use the size of land gridcells to resize the FSDS matrix

    landcells = len(grid_ids)    
    
    lon = np.multiply(mask,lon)
    lon = lon[lon != 0]
    
    lat = np.multiply(mask,lat)
    lat = lat[lat != 0]

    XC = np.multiply(mask,XC)
    XC = XC[XC != 0]
    """
    # Create a boolean mask where True indicates non-NaN values
    bool_mask = ~np.isnan(data[0])
    #land_mask = np.where(~np.isnan(data[0]), 1, 0)
    # Apply the mask to the data data array
    grid_ids = grid_ids[bool_mask]
    lon = lon[bool_mask]
    lat = lat[bool_mask]
    XC = XC[bool_mask]
    YC = YC[bool_mask]

    print("domain info" + str(lon.shape) + str(lat.shape) + str(XC.shape) \
            + str(YC.shape) + str(bool_mask.shape) + str(grid_ids.shape))
    return grid_ids, data, lon, lat, XC, YC, bool_mask


def data_partition_RR(number_of_subdomains, grid_ids, data):
    # cyclic (round-robin) partition
    domains = [[] for _ in range(number_of_subdomains)]
    for element, domain in zip(grid_ids, cycle(domains)):
        domain.append(element)

    grid_id_domains = domains.copy()
    
    landcells = len(grid_ids) 
    # partition the data over landcells
    # landcell_idx is alse the column_idx of FSDS
    landcell_idx = np.linspace(0, landcells-1, landcells, dtype=int)

    domains = [[] for _ in range(number_of_subdomains)]
    for element, domain in zip(landcell_idx, cycle(domains)):
        domain.append(element)
    
    # save the boundaries of each subdomain (for array_split)
    size_of_subdomains = [ len(domain) for domain in domains]

    # partitioned landcells_idx in subdomains 
    arranged_grid_idx = np.concatenate(domains).ravel()

    # find the original index of landcells for column swap
    np.sort(arranged_grid_idx)
    grid_swap_idx = (np.argsort(arranged_grid_idx))

    # create swap index and arrange data
    idx = np.empty_like(grid_swap_idx)
    idx[grid_swap_idx] = np.arange(len(grid_swap_idx))
    data = data[:,idx]

    # split the FSDS into subdomains using the boundary index
    subdomain_idx = size_of_subdomains
    for i in range(0,len(subdomain_idx)-1):
        subdomain_idx[i+1] +=subdomain_idx[i]
    
    #FSDS_list = np.hsplit(FSDS,subdomain_idx[:-1])
    data = np.hsplit(data,subdomain_idx[:-1])
    
    return grid_id_domains, data


def data_save(number_of_subdomains, grid_id_domains, subdomain_path, i_timesteps, \
              var_name, period, data):
    for i in range(number_of_subdomains):
        # convert local grid_id_lists into an array
        grid_id_arr = np.array(grid_id_domains[i])

        #data_arr = np.array(FSDS_list[i])
        data_arr = np.array(data[i])
        file_name = subdomain_path + 'subdomain'+ str(i) + '.' + var_name + '.' + period +'.nc'

        # Open a new NetCDF file to write the data to. For format, you can choose from
        # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
        w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
        w_nc_fid.title = 'The ELM domain files on individudal process: '+str(i)

        # create the gridIDs variable
        x_dim = w_nc_fid.createDimension('x_dim', grid_id_arr.size)
        time_dim = w_nc_fid.createDimension('time_dim', i_timesteps)
        w_nc_var = w_nc_fid.createVariable('gridIDs', np.int32, ('x_dim',))
        w_nc_var.long_name = 'gridIds in the subdomain'    
        w_nc_fid.variables['gridIDs'][:] = grid_id_arr.reshape(grid_id_arr.size)

        w_nc_var = w_nc_fid.createVariable(var_name, np.float32, ('time_dim', 'x_dim'))
        w_nc_var.long_name = 'FSDS in the subdomain'    
        w_nc_fid.variables[var_name][:] =data_arr.reshape(i_timesteps,grid_id_arr.size)
        w_nc_fid.close()  # close the new file

        
def data_save_1dNA(output_path, grid_ids, i_timesteps, var_name, period, data, lon, lat):
    # convert local grid_id_lists into an array
    grid_id_arr = np.array(grid_ids)

    #data_arr = np.array(FSDS_list[i])
    data_arr = np.array(data)
    file_name = output_path + 'clmforc.Daymet4.1km.1d.' + var_name + '.' + period +'.nc'

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
    w_nc_fid.title = 'forcing : '+ var_name + 'in the NA domain'

    # create the gridIDs, lon, and lat variable
    x_dim = w_nc_fid.createDimension('x_dim', grid_id_arr.size)
    time_dim = w_nc_fid.createDimension('time_dim', i_timesteps)
    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('x_dim',))
    w_nc_var.long_name = 'gridIds in the NA domain'    
    w_nc_fid.variables['gridIDs'][:] = grid_id_arr.reshape(grid_id_arr.size)

    w_nc_var = w_nc_fid.createVariable('lon', np.int32, ('x_dim',))
    w_nc_var.long_name = 'longitude of land gridcells in the NA domain'    
    w_nc_fid.variables['lon'][:] = lon
        
    w_nc_var = w_nc_fid.createVariable('lat', np.int32, ('x_dim',))
    w_nc_var.long_name = 'latitude of land gridcells in the NA domain'    
    w_nc_fid.variables['lat'][:] = lat
        
    # create the var_name variable
    w_nc_var = w_nc_fid.createVariable(var_name, np.float32, ('time_dim', 'x_dim'))
    w_nc_var.long_name = var_name + 'in the NA domain'    
    w_nc_fid.variables[var_name][:] =data_arr.reshape(i_timesteps,grid_id_arr.size)
        
    w_nc_fid.close()  # close the new file        


def domain_save_1dNA(output_path, grid_ids, data, lon, lat, XC, YC):
    # convert local grid_id_lists into an array
    grid_id_arr = np.array(grid_ids)

    #create land gridcell mask, area, and landfrac (all are 1s)
    mask = np.ones(grid_id_arr.size, dtype=np.int32)
    area = np.ones(grid_id_arr.size, dtype=np.float32)
    landfrac = np.ones(grid_id_arr.size, dtype=np.float32)
    file_name = output_path + 'Daymet4.1km.1d.domain.nc'
    print("The domain file is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
    w_nc_fid.title = '1D domain file for the Daymet NA region'

    # create the gridIDs, lon, and lat variable
    ni_dim = w_nc_fid.createDimension('ni', grid_id_arr.size)
    nj_dim = w_nc_fid.createDimension('nj', 1)
    nv_dim = w_nc_fid.createDimension('nv', 4)

    x_dim = w_nc_fid.createDimension('x_dim', 7814)
    y_dim = w_nc_fid.createDimension('y_dim', 8075)   

    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId in the NA domain'
    w_nc_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
    w_nc_fid.variables['gridID'][:] = grid_id_arr.reshape(grid_id_arr.size)

    w_nc_var = w_nc_fid.createVariable('xc', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    #w_nc_var.bounds = "xv" ;    
    w_nc_fid.variables['xc'][:] = lon
        
    w_nc_var = w_nc_fid.createVariable('yc', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"        
    w_nc_fid.variables['yc'][:] = lat
        
    # create the XC, YC variable
    w_nc_var = w_nc_fid.createVariable('xc_LCC', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'X of land gridcell center (Lambert Conformal Conic), increasing from west to east'
    w_nc_var.units = "m"
    #w_nc_var.bounds = "xv" ;    
    w_nc_fid.variables['xc_LCC'][:] = XC
        
    w_nc_var = w_nc_fid.createVariable('yc_LLC', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'Y of land gridcell center (Lambert Conformal Conic), decreasing from north to south'
    w_nc_var.units = "m"        
    w_nc_fid.variables['yc_LCC'][:] = YC

    w_nc_var = w_nc_fid.createVariable('area', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells (Lambert Conformal Conic)'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "km2"        
    w_nc_fid.variables['area'][:] = area

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'mask of land gridcells (1 means land)'
    w_nc_var.units = "unitless"        
    w_nc_fid.variables['mask'][:] = mask  

    w_nc_var = w_nc_fid.createVariable('landfrac', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'fraction of land gridcell that is active'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "unitless" ;        
    w_nc_fid.variables['landfrac'][:] = landfrac  

    w_nc_var = w_nc_fid.createVariable('lambert_conformal_conic', np.short)
    w_nc_var.grid_mapping_name = "lambert_conformal_conic"
    w_nc_var.longitude_of_central_meridian = -100. 
    w_nc_var.latitude_of_projection_origin = 42.5         
    w_nc_var.false_easting = 0. 
    w_nc_var.false_northing = 0.
    w_nc_var.standard_parallel = 25., 60.          
    w_nc_var.semi_major_axis = 6378137.
    w_nc_var.inverse_flattening = 298.257223563

    w_nc_fid.close()  # close the new file  


def domain_save_2dNA(output_path, total_rows, total_cols, data, lon, lat, XC, YC):

    # add the gridcell IDs.
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(total_cols,total_rows)

    #create land gridcell mask, area, and landfrac (otherwise 0)
    mask = np.where(~np.isnan(data[0]), 1, 0)
    area = mask.astype(float)
    landfrac = area 

    file_name = output_path + 'Daymet4.1km.2d.domain.nc'
    print("The 2D domain file is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
    w_nc_fid.title = '2D domain file for the Daymet NA region'

    # create the gridIDs, lon, and lat variable
    x_dim = w_nc_fid.createDimension('x_dim', 7814)
    y_dim = w_nc_fid.createDimension('y_dim', 8075)   

    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'gridId in the NA domain'
    w_nc_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
    w_nc_fid.variables['gridID'][:] = grid_ids 

    w_nc_var = w_nc_fid.createVariable('xc', np.float32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    #w_nc_var.bounds = "xv" ;    
    w_nc_fid.variables['xc'][:] = lon
        
    w_nc_var = w_nc_fid.createVariable('yc', np.float32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"        
    w_nc_fid.variables['yc'][:] = lat
        
    # create the XC, YC variable
    w_nc_var = w_nc_fid.createVariable('xc_LCC', np.float32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'X of land gridcell center (Lambert Conformal Conic), increasing from west to east'
    w_nc_var.units = "m"
    #w_nc_var.bounds = "xv" ;    
    w_nc_fid.variables['xc_LLC'][:] = XC
        
    w_nc_var = w_nc_fid.createVariable('yc_LCC', np.float32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'Y of land gridcell center (Lambert Conformal Conic), decreasing from north to south'
    w_nc_var.units = "m"        
    w_nc_fid.variables['yc_LLC'][:] = YC

    w_nc_var = w_nc_fid.createVariable('area', np.float32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'Area of land gridcells (Lambert Conformal Conic)'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "km2"        
    w_nc_fid.variables['area'][:] = area

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'mask of land gridcells (1 means land)'
    w_nc_var.units = "unitless"        
    w_nc_fid.variables['mask'][:] = mask  

    w_nc_var = w_nc_fid.createVariable('landfrac', np.float32, ('y_dim','x_dim'))
    w_nc_var.long_name = 'fraction of land gridcell that is active'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "unitless" ;        
    w_nc_fid.variables['landfrac'][:] = landfrac  


    w_nc_var = w_nc_fid.createVariable('lambert_conformal_conic', np.short)
    w_nc_var.grid_mapping_name = "lambert_conformal_conic"
    w_nc_var.longitude_of_central_meridian = -100. 
    w_nc_var.latitude_of_projection_origin = 42.5         
    w_nc_var.false_easting = 0. 
    w_nc_var.false_northing = 0.
    w_nc_var.standard_parallel = 25., 60.          
    w_nc_var.semi_major_axis = 6378137.
    w_nc_var.inverse_flattening = 298.257223563

    w_nc_fid.close()  # close the new file  


def launch_job(input_path, file_name, output_path, var_name, period, \
               number_of_subdomains, i_timesteps):
   
    #print('number of subdomains: ('+ str(number_of_subdomains) + \
    #      '), timeseries: ('+ str(i_timesteps) + '), file_name: (' + file_name + \
    #      '), input_path : (' + input_path + ')')

    start = process_time()
    file_name = input_path + file_name
    #print(input_path, file_name)
    [total_rows, total_cols, i_timesteps, data, lon, lat, XC, YC] = data_read(file_name, var_name, i_timesteps)
    end = process_time()
    print("Reading " + file_name + " takes  {}".format(end-start))

    start = process_time()
    domain_save_2dNA(output_path, total_rows, total_cols, data, lon, lat, XC, YC)
    end = process_time()
    print("Saving 2D domain data takes {}".format(end-start))
    
    start = process_time()  
    [grid_ids, data, lon, lat, XC, YC, land_mask] = domain_info_1D(total_rows, total_cols, i_timesteps, data, lon, lat, XC, YC)
    end = process_time()
    print("Creating 1D domain data takes  {}".format(end-start))

    """
    start = process_time()
    [grid_id_domains, data] = data_partition_RR(number_of_subdomains, grid_ids, data)
    end = process_time()
    print("Partitioning data/GridID in " + file_name + "takes  {}".format(end-start))
    """
    
    start = process_time() 
    #data_save(number_of_subdomains, grid_id_domains, output_path, i_timesteps, \
    #          var_name, period, data)
    #data_save_1dNA(output_path, grid_ids, i_timesteps, var_name, period, data, lon, lat)
    domain_save_1dNA(output_path, grid_ids, data, lon, lat, XC, YC)
    end = process_time()
    print("Saving 1D domain data takes {}".format(end-start))
    
def get_files(input_path):
    print(input_path)
    files = os.listdir(input_path) 

    files.sort() 

    file_no =0

    files_nc = [f for f in files if (f[-2:] == 'nc')] 
    print("total " + str(len(files_nc)) + " files need to be processed")
    return files_nc

    
def main():
    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    number_of_subdomains = int(arg[2])
    i_timesteps = int(arg[3])
    """
    
    input_path= './'
    file_name = 'TBOT.201401.nc'
    output_path = input_path
    number_of_subdomains = 1
    i_timesteps = 1
    var_name = 'TBOT'
    period= "2014"
    
    launch_job(input_path, file_name, output_path, var_name, period, number_of_subdomains, i_timesteps)

if __name__ == '__main__':
    main()
    