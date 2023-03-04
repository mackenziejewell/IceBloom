# DEPENDENCIES:
import os
from datetime import datetime
import xarray as xr
import pandas as pd
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs


def grab_monthly_sicCDR(year = 2000, month = 1, 
                main_path = '/Volumes/Jewell_EasyStore/IceBloomData/sicCDRNOAANSIDC/', 
                return_vars = ['xx', 'yy', 'sic', 'proj_info', 'ds', 'proj'], suppress_prints = True):
    
    """Grab daily sea ice concentrations from monthly V4 Climate Data Record NSIDC-G02202 sea ice concentration data (doi: 10.7265/efmz-2t65).

INPUT: 
- year: year of data to open (int)
- month: month of data to open (int)
- main_path: path to directory where SIC data are stored (string)
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['xx', 'yy', 'sic', 'proj_info', 'ds', 'proj']
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT: some or all of the following, as specified in return_vars
- xx: x-coordinates of projected SIC data (M x N array)
- yy: y-coordinates of projected SIC data (M x N array)
- sic: monthly CDR SIC data (M x N array)
- proj_info: projection info
- proj: cartopy projection from sea ice concentration data projection info
- ds: xarray data frame containing opened SIC data

DEPENDENCIES:
import os
from datetime import datetime
import xarray as xr
import pandas as pd
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs

Latest recorded update:
03-03-2022
    """
    
    # assert input variable types and options
    assert type(return_vars) == list, f'return_vars must be list, not {type(return_vars)}'
    assert len(return_vars) > 0, 'return_vars list is empty. Must have length >=1'
    assert type(main_path) == str, f'main_path must be string, not {type(main_path)}'
    assert os.path.isdir(main_path), f'"{main_path}" not an existing directory'
    assert type(suppress_prints) == bool, f'suppress_prints must be bool, not {type(suppress_prints)}'
    
    if suppress_prints == False:
        print(f' >>> opening {main_path}')
    
    # name_format = 'seaice_conc_daily_nh_{}_v04r00.nc'
    name_prefix = 'seaice_conc_monthly_nh_{}{}_'.format(str(year), str(month).zfill(2))

    if suppress_prints == False:
            print(f' >>> searching for file starting with {name_prefix}')

    # list files with name_prefix adjusted for provided date
    possible_files = [file for file in os.listdir(main_path) if file.startswith(name_prefix)]

    if suppress_prints == False:
        print('\n found following files:')
        for file in possible_files:
            print(f'  * {file}')
                
    # if only one file found, select as filename
    assert len(possible_files) == 1, f'should be single filename matching {name_prefix} found {len(possible_files)} matching files: {possible_files} in directory: {sic_datapath}'
    filename = possible_files[0]

    # open SIC data and select date
    data_path = main_path+filename
    ds = xr.open_dataset(data_path)
    ds.close()

    # find index of desired date and select
    ice_conc = ds.cdr_seaice_conc_monthly.values[0]

    if suppress_prints == False:
        print(f'\n time: {ds.time.values}')

    # store all output in dictionary
    vars_dict = {}
    
    # grab projection attributes, geo grid, and variables
    vars_dict['proj_info'] = ds.projection.attrs['spatial_ref']
    vars_dict['xx'], vars_dict['yy'] = np.meshgrid(ds.xgrid, ds.ygrid)
    vars_dict['sic'] = ma.masked_where(ice_conc>1, ice_conc)
    vars_dict['ds'] = ds
    
    # create ice projection from proj info
    #-------------------------------------
    # grab parameters from projection parameters
    semimajor = ds.projection.semimajor_radius
    semiminor = ds.projection.semiminor_radius
    standard_parallel = ds.projection.standard_parallel
    longitude_of_projection_origin = ds.projection.longitude_of_projection_origin
    # projection
    projection_uncropped = vars_dict['proj_info'][vars_dict['proj_info'].find('PROJECTION'):]
    projection = projection_uncropped[projection_uncropped.find("[")+2:projection_uncropped.find("]")-1]
    
    # create ice projection from info
    ice_projection = ccrs.NorthPolarStereo(central_longitude=int(longitude_of_projection_origin), true_scale_latitude=int(standard_parallel), globe=ccrs.Globe(semimajor_axis = semimajor, semiminor_axis = semiminor))
    vars_dict['proj'] = ice_projection
    
    if suppress_prints == False:
        print(f'\nProjection: {projection}')
        print(f' - standard_parallel: {standard_parallel}')
        print(f' - longitude_of_projection_origin: {longitude_of_projection_origin}') 
        
    # grab projection extent info if wanting to use with imshow
    #-----------------------------------------------------------
    xul = ds.projection.attrs['grid_boundary_left_projected_x']
    xlr = ds.projection.attrs['grid_boundary_right_projected_x']
    yul = ds.projection.attrs['grid_boundary_top_projected_y']
    ylr = ds.projection.attrs['grid_boundary_bottom_projected_y']
    EXTENT = [xul, xlr, yul, ylr]
    if suppress_prints == False:
        print(f'EXTENT: {EXTENT}')
    
    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]

    return return_data



# function for averaging SIC data
#--------------------------------
def calc_mean_sic(years_to_average = [2000], months_to_average = [1], main_path = '/Volumes/Jewell_EasyStore/IceBloomData/sicCDRNOAANSIDC/', suppress_prints = True):
    
    """Calculate mean SIC across given months and years. Account for different month lengths with time-weighted mean.

INPUT: 
- years_to_average: list/array of years across which to average (integers)
- months_to_average: list/array of months across which to average (integers)
- main_path: path to directory where SIC data are stored (string)
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- xx: x-coordinates of projected SIC data (M x N array)
- yy: y-coordinates of projected SIC data (M x N array)
- mean_SIC: time-average SIC data (M x N array)
- proj: cartopy projection from sea ice concentration data projection info

DEPENDENCIES:
import os
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs

homemade function: 
grab_monthly_sicCDR

Latest recorded update:
03-03-2022
    """
    

    # sum up SIC data over given years and months
    counter = 0

    # run through given months and years
    # add weighted sic to running sum
    for year in years_to_average:
        for month in months_to_average:

            # find weight of monthly average (by month length)
            if month in [4,6,9,11]:
                weight = 30
            elif month == 2:
                if year%4 == 0:
                    weight = 29
                else:
                    weight = 28
            else:
                weight = 31

            # weighted sum of SIC
            if counter == 0:
                
                # open and add sic data
                xx, yy, sic, proj = grab_monthly_sicCDR(year = year, month = month, main_path = main_path, 
                                                        return_vars = ['xx', 'yy', 'sic', 'proj'],
                                                        suppress_prints=suppress_prints)
                ice_conc_sum =  weight*sic
                sum_weights  = weight
                
            else:
                
                # open and add sic data
                sic = grab_monthly_sicCDR(year = year, month = month, main_path = main_path, 
                                          return_vars = ['sic'],
                                          suppress_prints=suppress_prints)
                ice_conc_sum += weight*sic
                sum_weights  += weight

            counter+=1

    # calculate weighted mean SIC
    mean_SIC = ice_conc_sum/sum_weights
    
    return xx, yy, mean_SIC, proj



# function for grabbing chlorophyll data
#-----------------------------------------
def grab_chlor(filename):
    
    """Grab NASA ocean color (chlorophyll) data.

INPUT: 
- filename: path to file

OUTPUT:
- lons: longitude (M x N array)
- lats: latitude (M x N array)
- Chlor_a: chlorophyll a concentrations (M x N array)
- ds: xarray data frame for chl-a data

DEPENDENCIES:
import xarray as xr
import numpy as np

Latest recorded update:
03-03-2022
    """
    
    # open dataset with xarray
    ds = xr.open_dataset(filename)
    ds.close()
    
    # shift from (-180,180) to (0,360) becuase of plotting issues around dateline
    lons = ds.lon.values
    lons[lons<0]+=360
    ds = ds.assign_coords({"lon": lons}).sortby("lon")

    # drop lons to range around Bering Sea and extract data
    ds = ds.sel(lon=slice(ds.lon.values[3500],ds.lon.values[6000]))
    lons, lats  = np.meshgrid(ds.lon.values, ds.lat.values)
    Chlor_a = ds.chlor_a.values
    
    return lons, lats, Chlor_a, ds

