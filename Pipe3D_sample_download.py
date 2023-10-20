#For downloading a set of MaNGA datacubes from the Pipe3D archive storage

#Necessary imports

import numpy as np
import urllib.request
import os
import warnings
import gzip
import shutil

#MaNGA_IDs whose Pipe3D datacubes are to be downloaded. Supply the address of the text file containing MaNGA_IDs in 'Address1'
#Each MaNGA_ID should be a string of the form "aaaa-aaaa"

MaNGA_IDs = np.genfromtxt('Address1', dtype = 'str')

def get_file(Manga_ID): #code to get fits file from the Pipe3D archive storage and save it in a local location. Input the MaNGA_ID as a string.
    
    First = Manga_ID.split('-')[0] #the part of the MaNGA ID before the hyphen
    Last = Manga_ID.split('-')[1] #the part of the MaNGA ID after the hyphen

    #For getting the standard name of the file on the Pipe3D archive storage -- this is a zip file inside which the fits file corresponding to the datacube is present
    ID = 'manga' + '-' + First + '-' + Last
    Ext = '.Pipe3D.cube.fits'
    filename = ID+Ext
    gzip_ext = '.Pipe3D.cube.fits.gz'
    zip_file = ID + gzip_ext
    std_url = 'https://data.sdss.org/sas/dr16/manga/spectro/pipe3d/v2_4_3/2.4.3/' #SDSS DR16
    complete_url = std_url + First + '/' + zip_file #the final URL which will be launched for an input MaNGA_ID
    
    destination = 'Address2' #Supply the address of the directory where the zipped datacubes are to be stored
    file_add = 'Address3' #Supply the address of the directory where the zipped datacubes are to extracted

    if(os.path.isfile(file_add) == False): #checks for the fits file not being present already
        if(os.path.isfile(destination) == False): #checks for the gzip file not being present already
            urllib.request.urlretrieve(complete_url, destination) #downloads the gzip file
        
        #extract the fits file from the gzip file
        with gzip.open(destination,"rb") as f_in, open(file_add,"wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        
    return file_add #returns the address of the fits file

#Looping over all MaNGA_IDs to get all the datacubes
for Manga_ID in MaNGA_IDs:
    get_file(Manga_ID)
    
#Done !

##############################################################################################################################################################
