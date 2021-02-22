"""
Created on Mon Jun  8 17:40:49 2020

@author: Alberto
"""

import os
import sys

class FileManager:
    
    def __init__(self, inputfolder, ouputfolder, tempfolder=None):
        
        if os.path.isdir(str(inputfolder)) == True:
            self.Ifolder = str(inputfolder)
        else:
            sys.exit('Input folder does not exist')
        if os.path.isdir(str(ouputfolder)) == True:
            self.Ofolder = str(ouputfolder)
        else:
            sys.exit('Output folder does not exist')
        if tempfolder == None:
            self.Tfolder=os.path.join(os.getcwd(), 'temp')
            os.mkdir(self.Tfolder)
        elif os.path.isdir(str(tempfolder)) == True:
            self.Tfolder = str(tempfolder)
        else:
            self.Tfolder=os.path.join(os.getcwd(), 'temp')
            os.mkdir(self.Tfolder)
            
    def get_source_products(self, folder, fileformat=None):
        '''returns the generator'''

        g = list(self.gen_file_from_folder(folder, fileformat=fileformat))
        return g
        
    def gen_file_from_folder(self, folder, fileformat=None):
        '''Generator retrieveing files in a folder'''

        for subdir, dirs, files in os.walk(folder):
            for item in files:
                if fileformat == None:
                    yield self.get_file_path(folder, item)
                else:
                    if fileformat in item:
                        print (item)
                        yield self.get_file_path(folder, item)
                    else:
                        continue
                
    def get_file_path(self, folder, filename):
        '''Creates a path to file within a folder'''
        
        path = os.path.join(str(folder), str(filename))
        return path
    
    def remove_folder_contents(self, folder):
        
        filelist = list(self.gen_file_from_folder(folder))
        for f in filelist:
            os.remove(f)
            
    def remove_dir(self, folder):
        
        os.rmdir(folder)
        
    def get_basename(self, path):
        
        return os.path.basename(path)
