"""
Created on Fri Jun 19 18:50:08 2020

@author: Alberto S. Rabaneda
"""

from file_explorer import FileManager
from preprocessing import SNAP
from winddir import WIND
from objects import Object
from waves import Wave

#----------------------------------------------------------------------------
#COMPULSORY VARIABLES

'''Folder where S1 GRD products with VV polarisation are located'''
input_folder_full_path = ''

'''Folder to save the Netcdf products of the processing chain'''
output_folder_full_path = ''

'''Set True for including object detction'''
include_objects = True

'''Set True for including wind layers'''
include_wind = True

'''Set True for including wave layers'''
include_waves = True 

#----------------------------------------------------------------------------
#OPTIONAL VARIABLES

'''Folder where to keep temprorarily the preprocessing outputs'''
tempfolder_full_path = ''

'''Set True if object intermidiate and quality layers are meant to be included
in the product NetCDF file'''
object_extra_layers = True

'''Set True if wind intermidiate and quality layers are meant to be included
in the product NetCDF file'''
wind_extra_layers = True

'''Set True if waves intermidiate and quality layers are meant to be included
in the product NetCDF file'''
wave_extra_layers = True

'''Metadata to add to the final netcdf file'''
glob_attr = {'Title':'',
             'Conventions':'CF-1.7',
             'Institution':'University of Hull',
             'Product_type':'Sentinel1 Product',
             'Origin':'Scihub',
             'Project':'DCS4COP',
             'Contact':'A.Rabaneda@hull.ac.uk or R.Forster@hull.ac.uk',
             'History':'Rabaneda Sentinel1 processing chain',
             'time_coverage_start':'',
             'time_coverage_end':'',
             'time_coverage_duration':'0',
             'time_coverage_resolution':'second',
             'geographic_coodinate_system':'EPSG:32631',
             'geospatial_lat_min':'',
             'geospatial_lat_max':'',
             'geospatial_lon_min':'',
             'geospatial_lon_min':'',
             'raster_width':'100',
             'raster_height':'None',
             'raster_resolution':'100',
             'raster_resolution_units':'metres',
             'version':'v1.0'}

#----------------------------------------------------------------------------
#RUN PROCESSING CHAIN

if __name__ == "__main__":
    
    print ('Starting computation of S1 processing chain')
    f = FileManager(input_folder_full_path, output_folder_full_path, tempfolder=tempfolder_full_path)
    s1_images = f.get_source_products(f.Ifolder)

    try:
        a=1
        for img in s1_images:
            
            print ('Preprocessing image '+str(a)+' out of '+str(len(s1_images)))
            pre = SNAP(img, f.Tfolder)
            if include_objects == True and include_wind == True and include_waves == True:
                pre.cmd_order(band='amplitude')
                pre.cmd_order(band='intensity')
                pre.cmd_order(band='waves')
                pre.cmd_order(band='slant')
                pre.cmd_order(band='objects')
            elif include_objects == False and include_wind == True and include_waves == True:
                pre.cmd_order(band='amplitude')
                pre.cmd_order(band='intensity')
                pre.cmd_order(band='waves')
                pre.cmd_order(band='slant')
            elif include_objects == True and include_wind == False and include_waves == True:
                pre.cmd_order(band='amplitude')
                pre.cmd_order(band='waves')
                pre.cmd_order(band='slant')
                pre.cmd_order(band='objects')
            elif include_objects == True and include_wind == True and include_waves == False:
                pre.cmd_order(band='amplitude')
                pre.cmd_order(band='intensity')
                pre.cmd_order(band='objects')
            elif include_objects == False and include_wind == False and include_waves == True:
                pre.cmd_order(band='amplitude')
                pre.cmd_order(band='waves')
                pre.cmd_order(band='slant')
            elif include_objects == False and include_wind == True and include_waves == False:
                pre.cmd_order(band='amplitude')
                pre.cmd_order(band='intensity')
            elif include_objects == True and include_wind == False and include_waves == False:
                pre.cmd_order(band='amplitude')
                pre.cmd_order(band='objects')
                
            preprocessed = f.get_source_products(f.Tfolder)
            for x in preprocessed:
                if x.endswith('_objects.nc'):
                    obj = x
                elif x.endswith('_amplitude.nc'):
                    amp = x
                elif x.endswith('_intensity.nc'):
                    sig = x
                elif x.endswith('_waves.nc'):
                    wav = x
                elif x.endswith('_slant.nc'):
                    sla = x
            
            print ('Processing image '+str(a)+' out of '+str(len(s1_images)))
            if sig and amp:
                #merging amplitude and sigma 100x100 products
                gen_amp = WIND(dataset=amp, output=None)
                gen_sig = WIND(dataset=sig, output=None)
                print(gen_amp.ds)
                print(gen_sig.ds)
                wind_ds = gen_amp.merge_datasets([gen_amp.ds, gen_sig.ds])
                
                #calculation of wind direction and speed
                print ('Calculating wind layers')
                wi = WIND(dataset = wind_ds, output = wind_ds, drop_var_output = True, inter_layer = wind_extra_layers)
                wi.get_direction_matrix()
                wi.get_speed_matrix()
                result = wi.return_ds(wi.out)
            else:
                #just prepare canvas
                gen_amp = WIND(dataset=amp, output=None)
                canvas = gen_amp.drop_var(gen_amp.ds, list(gen_amp.ds.data_vars))
                result = gen_amp.return_ds(canvas)
                
            if wav and sla:   
                #merging sigma and slant 10x10 products
                print ('Calculating wave layers')
                gen_wav = Wave(dataset=wav, output=None)
                gen_sla = Wave(dataset=sla, output=None)
                wave_ds = gen_wav.merge_datasets([gen_wav.ds, gen_sla.ds])
                
                #calculating wave parameters
                wa = Wave(dataset = wave_ds, output = result, inter_layer = wave_extra_layers)
                wa.get_wave_arrays()
                result = wa.return_ds(wa.out)
            
            if obj:
                #calculating objects
                print ('Calculating object layers')
                ob = Object(dataset = obj, output = result, inter_layer = object_extra_layers)
                ob.get_objects_array()
                result = ob.return_ds(ob.out)
            
            #saving netcdf with results
            fi = Object()
            final = fi.add_attributes(result, glob_attr, override=True)
            path = f.get_file_path(f.Ofolder, f.get_basename(img)+'_completed.nc')
            fi.save_ds(path, final)
            #f.remove_folder_contents(f.Tfolder)
            a += 1
        
    finally:
        print ('done')
        if tempfolder_full_path == None:
            f.remove_dir(f.Tfolder)
