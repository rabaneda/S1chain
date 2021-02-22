"""
Created on Wed Jun 17 17:24:15 2020

@author: Alberto S. Rabaneda
"""

import numpy as np
from nc_methods import NetCDFManager
import scipy.stats as st


class Object(NetCDFManager):
    
    kernels = {'50m':np.array([[1,1,1,1,1], [1,10,10,10,1], [1,10,0,10,1], [1,10,10,10,1], [1,1,1,1,1]]),
                '70m':np.array([[1,1,1,1,1,1,1], [1,5,5,5,5,5,1], [1,5,10,10,10,5,1], [1,5,10,0,10,5,1], [1,5,10,10,10,5,1], [1,5,5,5,5,5,1], [1,1,1,1,1,1,1]]),
                '90m':np.array([[1,1,1,1,1,1,1,1,1], [1,3,3,3,3,3,3,3,1], [1,3,5,5,5,5,5,3,1], [1,3,5,10,10,10,5,3,1], [1,3,5,10,0,10,5,3,1], [1,3,5,10,10,10,5,3,1], [1,3,5,5,5,5,5,3,1], [1,3,3,3,3,3,3,3,1], [1,1,1,1,1,1,1,1,1]])}
     
    magnitude_thresholds = {'50m':40,
                             '70m':90,
                             '90m':150}
    
    significance_threshold = 80 #as %
    intensity = 'Intensity_VV'
     
    obj_name = 'Object detected'
    obj_attr = {'Long_name': 'Human-built-Object detected', 
                'Standard_name':'Object detected',
                'units':'none',
                'resolution':'100m',
                'scale_factor':1}
    
    sig_name = 'Object significance'
    sig_attr = {'Long_name': 'Human-built-Object detected significance', 
                'Standard_name':'Object significance',
                'units':'none',
                'resolution':'100 m',
                'scale_factor':1}
    
    mag_name = 'Object magnitude'
    mag_attr = {'Long_name': 'Human-built-Object detected magnitude', 
                'Standard_name':'Object magnitude',
                'units':'none',
                'resolution':'100 m',
                'scale_factor':1}
    
    
    def get_significance(self):
        '''Return the significance of each pixel according to
        mean and stnadard deviation of the whole image'''
        
        x = self.get_var_array(self.ds, Object.intensity)
        mean = np.mean(x)
        std = np.std(x)
        arr = x[np.logical_not(np.isnan(x))]
        N = len(arr.flatten())
        S = (x-mean)*(N**(1/2.))/std
        return S
     
    def apply_threshold(self, array, th):
        '''Applies threshold to array
        
        array: numpy.array
        th: int, threshold as maximum value'''
        
        arr = array[array>th]
        return arr
    
    def get_S_threshold(self):
        '''Calculates and return the threshold value for significance'''
        
        arr = self.get_var_array(self.ds, Object.intensity).flatten()
        th = st.scoreatpercentile(arr, Object.significance_threshold)
        return th
    
    def get_objects_array(self, magnitude_precision='70m'):
        '''Return the object detection array calculated through significance
        and magnitude thresholds and kernel.
        
        magnitude_precision: str. Default is 70m, but could be 50m or 90m.
                            This parameter sets the precision to identify clusters
                            of pixels corresponding to an object'''
        
        Sig = self.get_significance()
        S = Sig
        S_th = self.get_S_threshold()
        
        kernel = Object.kernels[magnitude_precision]
        M_th = Object.magnitude_thresholds[magnitude_precision]

        S[S>S_th] = 1
        S[S<=S_th] = 0
        Mag = self.convolution_fourier(S, kernel)
        M=Mag
        M[M>=M_th] = 1
        M[M<M_th] = 0
        obj = S*M
        obj[obj==0] = np.nan
        
        self.add_var(self.out, Object.obj_name, self.downsampling_2D(obj, mode='median'), Object.obj_attr)
        if self.inter == True:
            self.add_var(self.out, Object.sig_name, self.downsampling_2D(Sig), Object.sig_attr)
            self.add_var(self.out, Object.mag_name, self.downsampling_2D(Mag), Object.mag_attr)
