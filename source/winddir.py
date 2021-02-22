"""
Created on Wed May 20 15:52:47 2020

@author: Alberto
"""

import numpy as np
import scipy.stats as st
from nc_methods import NetCDFManager
import warnings

#------------------------------------------------------------------------------

class WIND(NetCDFManager):
    
    kernels = {'op_sobel':np.array([[3,0,3], [10,0,-10], [3,0,-3]])*(1/32.)}
    ampli = 'Amplitude_VV'
    sigma = 'Sigma0_VV'
    incidence = 'incidenceAngleFromEllipsoid' # in degrees
    
    wspd_name = 'Wind speed'
    wspd_attr = {'Long_name': 'Neutral wind speed at 10 metres over the still water level', 
                'Standard_name':'Neutral wind speed at 10 m',
                'units':'m/s',
                'resolution':'1 km',
                'scale_factor':10}
    
    sigma_calc_name = 'NRCS calc'
    sigma_calc_attr = {'Long_name': 'Calculated Normalised Radar Cross Section', 
                'Standard_name':'Sigma nought calculated',
                'units':'m/s',
                'resolution':'100 m',
                'scale_factor':1}
    
    sigma_obs_name = 'NRCS obs'
    sigma_obs_attr = {'Long_name': 'Observed Normalised Radar Cross Section', 
                'Standard_name':'Sigma nought observed',
                'units':'m/s',
                'resolution':'100 m',
                'scale_factor':1}
    
    wdir_name = 'Wind direction'
    wdir_attr = {'Long_name': 'Wind direction with 180 degrees of ambiguity', 
                'Standard_name':'Wind direction',
                'units':'degrees',
                'resolution':'4 km',
                'scale_factor':40}
    
    R_name = 'Alignment'
    R_attr = {'Long_name': 'Mean Resultant Length', 
                'Standard_name':'Alignment',
                'units':'none',
                'resolution':'4 km',
                'scale_factor':40}
    
    ME_name = 'Marginal error'
    ME_attr = {'Long_name': 'Marginal Error of the Mean Resultant Vector', 
                'Standard_name':'Marginal Error',
                'units':'degrees',
                'resolution':'4 km',
                'scale_factor':40}
                
        
    def get_phase_matrix(self, amplitude_matrix):
        
        real = self.convolution_fourier(amplitude_matrix, WIND.kernels['op_sobel'])
        img = self.convolution_fourier(amplitude_matrix, np.transpose(WIND.kernels['op_sobel']))
        lg = real + (img*1j)
        phases = np.angle(lg)
        return phases
    
    def get_direction(self, arr, confidence):
        '''Returns mean direction, mean resultant vector and marginal error'''
        '''array must be contain axial directional data'''
        
        array = arr.flatten()
        angle = np.arctan2(np.mean(np.sin(array)), np.mean(np.cos(array)))*0.5
        R = np.power((np.mean(np.cos(array))**2)+(np.mean(np.sin(array))**2), 0.5)
        print (array.shape)
        print (np.mean(array)*confidence,np.mean(array)*(1-confidence))
        med = st.scoreatpercentile(array, 50, limit=(np.mean(array)*confidence,np.mean(array)*(1-confidence)))
        alpha = np.mean(np.cos(4*(array-angle)))
        ME = 0.5*(np.arcsin(med*np.power((1-alpha)/(2*len(array)*(R**2)), 0.5)))
        return (np.degrees(angle), R, ME)
    
    def get_direction_matrix(self, confidence=0.05, threshold=15, progressive_multilook=False):
        '''Returns mean direction array, mean resultant vector array and marginal error array'''
        ''' confidence: int, 0 to 1. Percintile to remove from its freq. distribution
                        Default 0.05, which means it will remove values within distribution
                        from 0 to 0.05 anf from 0.95 to 1.
            thershold: int, in degrees. Maximum marginal error in degrees to accept.
                        Default is 15 degrees.
            progessive multilook: Boolean, default is False, each imagette is
                        independent of the others, pixels belong to only one imagette.
                        If True, pixels will belong to multiple imagettes at the same time
                        since imagattes will overlap because an imagette is created for each 
                        pixel independently of multilook value; N of pixels = N of imagettes.
                        Each pixel will belong to multiple imagettes, but there will be
                        one imagette where this pixel will be the centre of the imagette'''
         
        axial = 2*self.get_phase_matrix(self.get_var_array(self.ds, WIND.ampli))
        angle_matrix = np.zeros(shape=axial.shape)
        R_matrix = np.zeros(shape=axial.shape)
        ME_matrix = np.zeros(shape=axial.shape)
                    
        subimages = list(self.gen_imagettes(axial, multilook=WIND.wdir_attr['scale_factor'], progressive_multilook=progressive_multilook))
        for roi in subimages:
            angle, R, ME = self.get_direction(roi[0], confidence=confidence)
            if ME > threshold:
                angle = np.nan
        
            if len(roi[1]) == 2:
                angle_matrix[roi[1][0],roi[1][1]] = angle
                R_matrix[roi[1][0], roi[1][1]] = R
                ME_matrix[roi[1][0], roi[1][1]] = ME
            elif len(roi[1]) == 4:
                angle_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = angle
                R_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = R
                ME_matrix[roi[1][0]:roi[1][1], roi[1][2]:roi[1][3]] = ME
        
        
        self.add_var(self.out, WIND.wdir_name, angle_matrix, WIND.wdir_attr)
        if self.inter == True:
            self.add_var(self.out, WIND.R_name, R_matrix, WIND.R_attr)
            self.add_var(self.out, WIND.ME_name, ME_matrix, WIND.ME_attr)

    
    def cmod5n_forward(self,v,phi,theta):
        '''!     ---------
        !     cmod5n_forward(v, phi, theta)
        !         inputs:
        !              v     in [m/s] wind velocity (always >= 0)
        !              phi   in [deg] angle between azimuth and wind direction
        !                    (= D - AZM)
        !              theta in [deg] incidence angle
        !         output:
        !              CMOD5_N NORMALIZED BACKSCATTER (LINEAR)
        !
        !        All inputs must be Numpy arrays of equal sizes
        !---------------------------------------------------------------------
           '''
        # Ignore overflow errors for wind calculations over land
        warnings.simplefilter("ignore", RuntimeWarning) 
        
        DTOR   = 57.29577951
        THETM  = 40.
        THETHR = 25.
        ZPOW   = 1.6
        
        # NB: 0 added as first element below, to avoid switching from 1-indexing to 0-indexing
        C = [0, -0.6878, -0.7957,  0.3380, -0.1728, 0.0000,  0.0040, 0.1103, 0.0159, 
              6.7329,  2.7713, -2.2885, 0.4971, -0.7250, 0.0450, 
              0.0066,  0.3222,  0.0120, 22.7000, 2.0813,  3.0000, 8.3659,
              -3.3428,  1.3236,  6.2437,  2.3893, 0.3249,  4.1590, 1.6930]
        Y0 = C[19]
        PN = C[20]
        A  = C[19]-(C[19]-1)/C[20]
    
        B  = 1./(C[20]*(C[19]-1.)**(3-1))
    
    #  !  ANGLES
        FI=phi/DTOR
        CSFI = np.cos(FI)
        CS2FI= 2.00 * CSFI * CSFI - 1.00
    
        X  = (theta - THETM) / THETHR
        XX = X*X
    
        #  ! B0: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
        A0 =C[ 1]+C[ 2]*X+C[ 3]*XX+C[ 4]*X*XX
        A1 =C[ 5]+C[ 6]*X
        A2 =C[ 7]+C[ 8]*X
    
        GAM=C[ 9]+C[10]*X+C[11]*XX
        S0 =C[12]+C[13]*X
        
        # V is missing! Using V=v as substitute, this is apparently correct
        V=v
        S = A2*V
        S_vec = S.copy() 
        SlS0 = [S_vec<S0]
        S_vec[SlS0]=S0[SlS0]
        A3=1./(1.+np.exp(-S_vec))
        SlS0 = (S<S0)
        A3[SlS0]=A3[SlS0]*(S[SlS0]/S0[SlS0])**( S0[SlS0]*(1.- A3[SlS0]))
        #A3=A3*(S/S0)**( S0*(1.- A3))
        B0=(A3**GAM)*10.**(A0+A1*V)
            
        #  !  B1: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
        B1 = C[15]*V*(0.5+X-np.tanh(4.*(X+C[16]+C[17]*V)))
        B1 = C[14]*(1.+X)- B1
        B1 = B1/(np.exp( 0.34*(V-C[18]) )+1.)
    
        #  !  B2: FUNCTION OF WIND SPEED AND INCIDENCE ANGLE
        V0 = C[21] + C[22]*X + C[23]*XX
        D1 = C[24] + C[25]*X + C[26]*XX
        D2 = C[27] + C[28]*X
    
        V2 = (V/V0+1.)
        V2ltY0 = V2<Y0
        V2[V2ltY0] = A+B*(V2[V2ltY0]-1.)**PN
        B2 = (-D1+D2*V2)*np.exp(-V2)
    
        #  !  CMOD5_N: COMBINE THE THREE FOURIER TERMS
        CMOD5_N = B0*(1.0+B1*CSFI+B2*CS2FI)**ZPOW
        return CMOD5_N
        
    
    def cmod5n_inverse(self, sigma0_obs, phi, incidence, iterations=10):
        '''!     ---------
        !     cmod5n_inverse(sigma0_obs, phi, incidence, iterations)
        !         inputs:
        !              sigma0_obs     Normalized Radar Cross Section [linear units]
        !              phi   in [deg] angle between azimuth and wind direction
        !                    (= D - AZM)
        !              incidence in [deg] incidence angle
        !              iterations: number of iterations to run
        !         output:
        !              Wind speed, 10 m, neutral stratification 
        !
        !        All inputs must be Numpy arrays of equal sizes
        !
        !    This function iterates the forward CMOD5N function
        !    until agreement with input (observed) sigma0 values   
        !---------------------------------------------------------------------
           '''
        # Ignore overflow errors for wind calculations over land
        warnings.simplefilter("ignore", RuntimeWarning) 
        # First guess wind speed
        V = np.array([10.])*np.ones(sigma0_obs.shape);
        step=10.
        
        # Iterating until error is smaller than threshold
        for iterno in range(1, iterations):
            #print iterno
            sigma0_calc = self.cmod5n_forward(V, phi, incidence)
            ind = sigma0_calc-sigma0_obs>0
            V = V + step
            V[ind] = V[ind] - 2*step 
            step = step/2
            
        #mdict={'s0obs':sigma0_obs,'s0calc':sigma0_calc}
        #from scipy.io import savemat
        #savemat('s0test',mdict)
        
        if self.inter == False:
            return (V)
        elif self.inter == True:
            return (V, sigma0_obs, sigma0_calc)
    
    def get_speed_matrix(self):
        
        azimuth = float(self.ds.attrs['azimuth_direction'])
        if self.inter == False:
            speed = self.cmod5n_inverse(self.downsampling_2D(self.get_var_array(self.ds, WIND.sigma), multilook=WIND.wspd_attr['scale_factor']), 
                                        self.downsampling_2D(self.get_var_array(self.out, WIND.wdir_name), multilook=WIND.wspd_attr['scale_factor'])-azimuth, 
                                        self.downsampling_2D(self.get_var_array(self.ds, WIND.incidence)), multilook=WIND.wspd_attr['scale_factor'])
            self.add_var(self.out, WIND.wspd_name, speed, WIND.wspd_attr)
        elif self.inter == True:
            speed, sigma0_obs, sigma0_calc = self.cmod5n_inverse(self.downsampling_2D(self.get_var_array(self.ds, WIND.sigma), multilook=WIND.wspd_attr['scale_factor']), 
                                        self.downsampling_2D(self.get_var_array(self.out, WIND.wdir_name), multilook=WIND.wspd_attr['scale_factor'])-azimuth, 
                                        self.downsampling_2D(self.get_var_array(self.ds, WIND.incidence), multilook=WIND.wspd_attr['scale_factor']))
            self.add_var(self.out, WIND.wspd_name, speed, WIND.wspd_attr)
            self.add_var(self.out, WIND.sigma_obs_name, sigma0_obs, WIND.sigma_obs_attr)
            self.add_var(self.out, WIND.sigma_calc_name, sigma0_calc, WIND.sigma_calc_attr)
