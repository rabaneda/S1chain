"""
Created on Mon Jun 15 18:13:22 2020

@author: Alberto S. Rabaneda
"""

from nc_methods import NetCDFManager
import numpy as np
import scipy.signal as signal

class Wave(NetCDFManager):
    
    sigma = 'Sigma0_VV'
    incidence = 'incidenceAngleFromEllipsoid' # in degrees
    slant = 'Slant_range'
    S1_V = 7502.4 # in m/s
    
    direc_name = 'Wave direction'
    direc_attr = {'Long_name': 'Wave direction with 180 degrees of ambiguity', 
                'Standard_name':'Wave direction',
                'units':'degrees',
                'resolution':'10 km',
                'scale_factor':100}
    
    phi_calc_name = 'Phi angle'
    phi_calc_attr = {'Long_name': 'Angle between wave direction and azimuth', 
                'Standard_name':'Phi angle',
                'units':'degrees',
                'resolution':'1 km',
                'scale_factor':10}
    
    cutoff_name = 'Cutoff'
    cutoff_attr = {'Long_name': 'Cutoff wavelength', 
                'Standard_name':'Cutoff wavelength',
                'units':'m',
                'resolution':'1 km',
                'scale_factor':10}
    
    height_name = 'Significant wave height'
    height_attr = {'Long_name': 'Significant wave height', 
                'Standard_name':'Significant wave height',
                'units':'m',
                'resolution':'1 km',
                'scale_factor':10}
    
    period_name = 'Wave period'
    period_attr = {'Long_name': 'Mean wave period', 
                'Standard_name':'Wave period',
                'units':'s',
                'resolution':'1 km',
                'scale_factor':10}
    
    def get_wave_arrays(self, pgr=False):
        '''Calculates wave direction (180 degrees of amiguity), 
        significant wave height and wave period.
        
            prg: Boolean, default is False, each imagette is
                        independent of the others, pixels belong to only one imagette.
                        If True, pixels will belong to multiple imagettes at the same time
                        since imagattes will overlap because an imagette is created for each 
                        pixel independently of multilook value; N of pixels = N of imagettes.
                        Each pixel will belong to multiple imagettes, but there will be
                        one imagette where this pixel will be the centre of the imagette'''
        
        
        sigma = self.get_var_array(self.ds, Wave.sigma)
        inci_matrix = self.get_var_array(self.ds, Wave.incidence)
        beta_matrix = self.get_var_array(self.ds, Wave.slant)/Wave.S1_V
        direc_matrix = np.zeros(shape=sigma.shape)
        phi_matrix = np.zeros(shape=sigma.shape)
        cutoff_matrix = np.zeros(shape=sigma.shape)
        height_matrix = np.zeros(shape=sigma.shape)
        period_matrix = np.zeros(shape=sigma.shape)
        
        subimages = list(self.gen_imagettes(sigma, multilook=Wave.direc_attr['scale_factor'], progressive_multilook=pgr))
        for roi in subimages:
            direc = np.angle(np.max(np.fft.fft2(roi[0])), deg=True)
            
            phi = direc - float(self.ds.attrs['azimuth_direction'])
            
            psd = signal.csd(roi[0], np.transpose(roi[0]),scaling='density')
            phases = (np.random.randn(psd.shape[0]*psd.shape[1]).reshape(psd.shape))*2*np.pi
            acf = np.fft.ifft2(np.sqrt(psd*2)*np.exp(1j*phases))
            std = np.std(signal.medfilt2d(np.abs(acf), kernel_size=11))
            cutoff = np.sqrt(2*np.pi*std)
            
            if len(roi[1]) == 2:
                sub_inci = np.deg2rad(np.mean(inci_matrix[roi[1][0],roi[1][1]]))
                sub_beta = np.mean(beta_matrix[roi[1][0],roi[1][1]])
                
                hei = ((cutoff/sub_beta)*(0.48 + 0.26*np.sin(sub_inci) + 0.27*np.cos(2*np.deg2rad(phi)))) + 0.22
                period = (hei*(sub_beta/cutoff)*1.65) + 5.60
                
                direc_matrix[roi[1][0],roi[1][1]] = direc
                phi_matrix[roi[1][0],roi[1][1]] = phi
                cutoff_matrix[roi[1][0], roi[1][1]] = cutoff
                height_matrix[roi[1][0], roi[1][1]] = hei
                period_matrix[roi[1][0], roi[1][1]] = period
                
            elif len(roi[1]) == 4:
                sub_inci = np.deg2rad(np.mean(inci_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]]))
                sub_beta = np.mean(beta_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]])
                
                hei = ((cutoff/sub_beta)*(0.48 + 0.26*np.sin(sub_inci) + 0.27*np.cos(2*np.deg2rad(phi)))) + 0.22
                period = (hei*(sub_beta/cutoff)*1.65) + 5.60
                
                direc_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = direc
                phi_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = phi
                cutoff_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = cutoff
                height_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = hei
                period_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = period
        
        self.add_var(self.out, Wave.direc_name, self.downsampling_2D(direc_matrix, multilook=Wave.direc_attr['scale_factor'], mode='angle'), Wave.direc_attr)
        self.add_var(self.out, Wave.height_name, self.downsampling_2D(height_matrix, multilook=Wave.height_attr['scale_factor']), Wave.height_attr)
        self.add_var(self.out, Wave.period_name, self.downsampling_2D(period_matrix, multilook=Wave.period_attr['scale_factor']), Wave.period_attr)
        if self.inter == True:
            self.add_var(self.out, Wave.phi_name, self.downsampling_2D(phi_matrix, multilook=Wave.phi_attr['scale_factor'], mode='angle'), Wave.phi_attr)
            self.add_var(self.out, Wave.cutoff_name, self.downsampling_2D(cutoff_matrix, multilook=Wave.cutoff_attr['scale_factor']), Wave.cutoff_attr)
