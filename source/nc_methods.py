"""
Created on Mon Jun  8 18:39:56 2020

@author: Alberto S. Rabaneda
"""

import xarray as xr
import numpy as np
import sys
import os
import math
from scipy.signal import fftconvolve, convolve2d

class NetCDFManager:
    
    def __init__(self, dataset = None, output = None, drop_var_output = False, inter_layer = False):
        '''Dataset can be mem object (xr.Dataset), nc or zarr format.
        In case of nc or zarr formats are passed, these needs to be the 
        full absolute path to the file'''
        '''Set inter_layer = True to retrieve intermidate layers, such as
        flags, or quality layers'''
        
        self.inter = inter_layer
        if dataset != None:
            self.ds = self.assert_dataset(dataset)
            if output == None:
                #self.out = self.drop_var(self.ds, list(self.ds.data_vars))
                pass
            elif output != None and drop_var_output == False:
                self.out = self.assert_dataset(output)
            elif output != None and drop_var_output == True:
                self.out = self.drop_var(self.ds, list(self.assert_dataset(output).data_vars))
        else:
            pass
        
    def assert_dataset(self, dataset):
        
        if isinstance(dataset, xr.Dataset) == True:
            return dataset
        elif isinstance(dataset, str) == True:
            if os.path.isabs(dataset) == True:
                return self.open_ds(dataset)
            else:
                sys.exit('Please, give absolute, valid path to dataset instead of relative path')
        else:
            sys.exit('Invalid Dataset format')
        
        
    def open_ds(self, filepath):

        if isinstance(filepath, str) == True:
            if filepath.endswith('.zarr') == True:
                ds = xr.open_zarr(filepath)
            elif  filepath.endswith('.nc') == True:
                try:
                    ds = self.beam_to_cf(xr.open_dataset(filepath))
                except KeyError:
                    ds = xr.open_dataset(filepath)
            else:
                sys.exit('Wrong dataset file format')
        else:
            sys.exit('Invalid filepath to dataset')
        return ds
    
    def save_ds(self, filepath, ds):
        
        if isinstance(filepath, str) == True:
            if filepath.endswith('.zarr') == True:
                ds.to_zarr(filepath)
            elif  filepath.endswith('.nc') == True:
                ds.to_netcdf(filepath)
            else:
                sys.exit('Wrong dataset file format')
        else:
            sys.exit('Invalid filepath to dataset')
            
    def return_ds(self, ds):
        '''Just return the created  dataset'''
        
        return ds
            
    def drop_var(self, ds, var_list):
        '''Method to remove unwanted varariables
            
            var_list: List with variable names to drop
            '''

        ds2 = ds.drop_vars(var_list)
        return ds2
    
    def get_var_array(self, ds, var_name):
        '''Returns the amplitude variable as numpy array'''
        
        ampli_array = ds[var_name].values
        return ampli_array
    
    def concat_datasets(self, files, dim):
        '''Method to concat netcdffiles along a dimension.
        
            files: string list of absolute paths to datasets
            dim: name of dimension,  string'''
        
        filepaths = [self.assert_dataset(x) for x in files]
        print ('>>>>>>>>>>>>>>>>>> Appending file '+filepaths[0]+'>>>>>>>>>>>>>>>>>>>>>>>>')
        first = self.open_ds(filepaths[0])
        for i in range(1, len(filepaths), 1):
            print ('>>>>>>>>>>>>>>>>>> Appending file '+filepaths[i]+'>>>>>>>>>>>>>>>>>>>>>>>>')
            first = xr.concat([first, self.open_ds(filepaths[i])], dim = dim)
        ds2 = first.sortby(dim)
        return ds2
    
    def merge_datasets(self, filepaths, overriding=False):
        '''Merge different datasets of same dimensions into one dataset
            that will contain all variables from datasets without repeating
            variables.
            
            filepaths: list of strings, absolute paths, or list of datasets in the mem'''
            
        dss = [self.assert_dataset(x) for x in filepaths]
        if overriding == False:
            try:
                ds2 = xr.merge(dss)
            except xr.MergeError:
                sys.exit('Conflict between datasets: Same variable name with different values')
        elif overriding == True:
            for i in range(1, len(dss), 1):
                dss[0] = ds2
                ds2.update(dss[i])
        return ds2
        
    def beam_to_cf(self, ds):
        '''This method will convert a xarray.dataset from SNAP-beam to cf convention'''
        
        for i in dict(ds['metadata'].attrs).keys():
            if 'platformHeading' in i:
                ds.attrs['azimuth_direction'] = float(ds['metadata'].attrs[i])
                break
            else:
                continue
        ds.attrs['crs'] = ds['crs'].attrs['wkt']
        ds.drop(['crs', 'metadata'])
        return ds
    
    def subset(self, ds, bbox):
        '''Extract a subset of a dataset
            
            ds: xarray dataset
            bbox: list of coordinates as [lat_min, lat_max, lon_min, lon_max]
            lat:[-90, 90] lon[-180, 180] degrees.'''
        
        try:
            polygon = ds.sel(latitude = slice(bbox[0], bbox[1]), longitude = slice(bbox[2], bbox[3]))
        except KeyError:
            polygon = ds.sel(lat = slice(bbox[0], bbox[1]), lon = slice(bbox[2], bbox[3]))
        return polygon
    
    def add_attributes(self, ds, dictio, override = False):
        '''Add attributes to dataset
        
            dictio: dictionary with attributes,
                    key = attribute, value = descrition
            override: if True, drops current attributes and set dictio as
                    new attributes. Default is False.
            '''
            
        if override == False:
            atts = dict(self.ds.attrs).update(dictio)
        elif override == True:
            atts = dictio
        del ds.attrs
        ds.attrs = atts
        return ds
        
    def add_var_attributes(self, ds, dictio, varname, override = False):
        '''Add attributes to variable/coordinate of the dataset
        
            ds: xarray.Dataset, dataset to add the variable
            dictio: dictionary with attributes,
                    key = attribute, value = description
            varname: string, name of the var to modify its attributes
            override: if True, drops current attributes and set dictio as
                    new attributes. Default is False.
            '''
            
        if override == False:
            atts = dict(self.ds[varname].attrs).update(dictio)
        elif override == True:
            atts = dictio
        del ds[varname].attrs
        ds[varname].attrs = atts
        
    def add_var(self, ds, var_name, var_array, var_dictio):
        '''Add variable to dataset.
        
        ds: xarray.Dataset, dataset to add the variable
        var_name: str, name of variable
        var_array: array, array of values
        var_dictio: dictionary with variable attributes'''
        
        prevs = list(ds.data_vars)
        if prevs[0].shape == var_array.shape:
            pass
        else:
            sys.exit('Variable array has different shape')
            
        ds[var_name] = var_array
        self.add_var_attributes(ds, var_dictio, var_name, override=True)
        return ds
        
    def gen_imagettes(self,ds, multilook=10, progressive_multilook=False):
        '''Generates subimages of a given dataset
            
            ds: array, array of values from a dataset variable
            multilook: int, number of pixels to cluster together.
                        Default is 10, then cluster will be 10x10 pixels.
            progessive multilook: Boolean, default is False, each imagette is
                        independent of the others, pixels belong to only one imagette.
                        If True, pixels will belong to multiple imagettes at the same time
                        since imagattes will overlap because an imagette is created for each 
                        pixel independently of multilook value; N of pixels = N of imagettes.
                        Each pixel will belong to multiple imagettes, but there will be
                        one imagette where this pixel will be the centre of the imagette'''
        
        if (multilook % 2) == 0:
            pass
        else:
            sys.exit('multilook must be an even number')
        
        if progressive_multilook == True:
            for i in range(ds.shape[0]):
                for j in range(ds.shape[1]):
                    
                    if i-(multilook//2) < 0:
                        s_i = 0
                    else:
                        s_i = i-(multilook//2)
                    if j-(multilook//2) < 0:
                        s_j = 0
                    else:
                        s_j = j-(multilook//2)
                    roi = ds[s_i:i+multilook//2, s_j:j+multilook//2]
                    coords = (i,j)
                    yield (roi, coords)
                    
        elif progressive_multilook == False:
            for i in range(multilook//2, ds.shape[0], multilook):
                for j in range(multilook//2, ds.shape[1], multilook):
                    
                    roi = ds[i-multilook//2:i+multilook//2, j-multilook//2:j+multilook//2]
                    coords = (i-multilook//2, i+multilook//2, j-multilook//2, j+multilook//2)
                    yield (roi, coords)
                    
    def convolution_2D (self, matrix, kernel):
        '''Returns the convolution of matrix*kernel, using determined from sums.
        Slow for big pics or many pics'''
        
        product = convolve2d(matrix, kernel, mode='same')
        return product
    
    def convolution_fourier (self, matrix, kernel):
        '''Returns the convolution of matrix*kernel, use fast fourier transformation'''
    
        product = fftconvolve(matrix, kernel, mode='same')
        return product
    
    def convolution_torch (self, matrix, kernel):
        '''Returns the convolution of matrix*kernel, useful when a GPU is used'''
    
        '''Under construction'''
        
    def downsampling_2D (self, array, multilook=10, mode='mean'):
        '''To downsample a 2D array
            
            array: array-like, preferible a numpy array.
            multilook: factor to reduce resolution, default to 10.
            mode: 'mean', 'median' or 'angle', default to 'mean'. '''
       
        new_size = array[::multilook, ::multilook].shape
        new_array = np.zeros(shape = new_size)
        if multilook == 1:
            return array
        else:
            for i in range(multilook//2, array.shape[0], multilook):
                for j in range(multilook//2, array.shape[1], multilook):
                    roi = array[i-multilook//2:i+multilook//2, j-multilook//2:j+multilook//2]
                    if np.count_nonzero(np.isnan(roi)) >= (multilook**2)//2:
                        if mode == 'median':
                            pix = np.median(roi[np.logical_not(np.isnan(roi))])
                        elif mode =='angle':
                            pix = self.avg_angles(roi[np.logical_not(np.isnan(roi))].flatten())
                        else:
                            pix = np.mean(roi[np.logical_not(np.isnan(roi))])
                    else:
                        pix = np.nan
                    new_array[i-multilook//2, j-multilook//2] = pix
            
            return new_array
        
    def avg_angles(self, *args):
        '''average of list of angles in degrees'''
        
        radians=[math.radians(x) for x in args]
        si=(sum([math.sin(x) for x in radians]))/len(args)
        co=(sum([math.cos(x) for x in radians]))/len(args)
        pri=math.degrees(math.atan2(si,co))
        if pri<0:
            ang=pri+360
        else:
            ang=pri
        return ang
