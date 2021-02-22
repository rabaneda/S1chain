#-----------------------------------------------------------------------------

class SNAP:
    '''Class to use snap for command line via gpt, directly from Python'''
    
    graph_amplitude = 'Preprocessing_amplitude.xml'
        #variable outputs: Amplitude_VV, incidence angle
        #pixel resolution = 100x100 m 
    
    graph_intensity = 'Preprocessing_sigmaintensity.xml'
         #variable outputs: sigma_vv from intensity
         #pixel resolution = 100x100 m
         
    graph_waves = 'Preprocessing_sigmaintensity_no_multilook.xml'
         #variable outputs: sigma_vv from intensity, incidence angle
         #pixel resolution = 10x10 m
         
    graph_objects = 'Preprocessing_intensity.xml'
         #variable outputs: intensity with speckle filter
         #pixel resolution = 10x10 m
         
    graph_slant_range = 'Preprocessing_slantrange.xml'
         #variable outputs:slant range distance
         #pixel resolution = 10x10 m
         
         
    def __init__(self, inputpath, ouputfolder, snap_dir=None):
        
        if os.path.isfile(inputpath) == True:
            self.Ifile = inputpath
        else:
            sys.exit('Input file does not exists')
        if os.path.isdir(ouputfolder) == True:
            self.Ofolder=ouputfolder
        else:
            sys.exit('Output folder does not exists')
        if snap_dir != None:
            if os.path.exists(snap_dir) == True:
                os.environ['PATH'] += snap_dir+':$PATH' 
            else:
                sys.exit('Invalid snap/bin path. Add .snap/bin to environment path')
        
    def cmd_order(self, band='amplitude'):
        '''Defines the order to use in the commnad line'''
        
        if band == 'amplitude':
            operations = SNAP.graph_amplitude
            end = '_amplitude.nc'
        elif band =='intensity':
            operations = SNAP.graph_intensity
            end = '_intensity.nc'
        elif band =='waves':
            operations = SNAP.graph_waves
            end = '_waves.nc'
        elif band =='objects':
            operations = SNAP.graph_objects
            end = '_objects.nc'
        elif band =='slant':
            operations = SNAP.graph_slant_range
            end = '_slant.nc'
        else:
            sys.exit('Invalid band name. Only amplitude and intensity can be passed as band')
            
        outputpath = os.path.join(self.Ofolder, os.path.basename(self.Ifile)[:-4]+end)
        line = ''.join(['gpt ', str(operations), ' -Pinput="', self.Ifile, '" -Poutput="', outputpath, '"'])
        #line = ['gpt', self.operations, '-Pinput="'+inputpath+'"', '-Poutput="'+outputpath+'"']
        subprocess.run(line, shell=True)

#--------------------------------------------------------------------------
