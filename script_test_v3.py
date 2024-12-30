import pci
from pci.pcimod import pcimod
from pci.sieve import *
from pci.speclass import speclass
import numpy as np
from pci.api import gobs, datasource as ds
import os
import time

import tkinter as tk
from tkinter import filedialog


################################### PROCESSING FUNCTIONS ############################################################################

#make_ndwi: Takes a multispectral image, applies NDWI index calculation, and creates a new .pix file containing index
### in_img: the filepath of the multispectral image (input)
### out_img: the filepath of the NDWI image to be created(output)
### bandnum1: the band (either coastal or green)
### bandnum2: the second band (either NIR, NIR-1, or NIR-2)
def make_ndwi(in_img, out_img, bandnum1, bandnum2):
    print("--Starting Make NDWI Function--")

    if os.path.isfile(out_img):
        os.remove(out_img)

    

    ### Read in original raster 
    with ds.open_dataset(in_img) as dataset:

        #Read in the channels you need (1 = coastal, 8 = NIR2)
        reader = ds.BasicReader(dataset, [bandnum1, bandnum2])
        in_aux = reader.aux_data
        in_coord_sys = reader.crs
        in_geocode = reader.geocoding
        
        in_raster = reader.read_raster(0, 0, reader.width, reader.height)
        print ("Finished Reading in raster")

    ### Give the bands an appropriate name 
    coastal_band = in_raster.data[:, :, 0]
    coastal_band = np.float32(coastal_band)
    coastal_band[coastal_band <= 0] = np.nan

    nir_band = in_raster.data[:, :, 1]
    nir_band = np.float32(nir_band)
    nir_band[nir_band <= 0] = np.nan
    
    print("Finished converting from 16S to 32R")

    ### Calculate NDWI
    
    with np.errstate(divide='ignore', invalid='ignore'):
        ndwi_array = np.nan_to_num((coastal_band - nir_band)/(coastal_band + nir_band), nan = -32768)

    ndwi_raster = gobs.array_to_raster(ndwi_array)
    print("Finished calculating NDWI")


    ### Output
    with ds.new_dataset(out_img, 'PCIDSK', '') as write_dataset:
        writer = ds.BasicWriter(write_dataset)

        writer.create(ndwi_raster)
        writer.write_raster(ndwi_raster)

    
        in_aux.set_file_metadata_value('NO_DATA_VALUE', '-32768')
        mdmap = [['BandDescription', 'NDWI'], ['NO_DATA_VALUE', '-32768']]
        in_aux.set_chan_metadata(mdmap, 1)
        in_aux.set_chan_description('NDWI', 1)
        writer.aux_data = in_aux
        
        writer.crs = in_coord_sys
        writer.geocoding = in_geocode

    print("--NDWI Function Finished--")
    return




#make_mask: Taking an image with NDWI, applies a threshold value to specify only values bigger than the threshold is given a "1" value; all others NO_DATA
### in_img: the filepath of the NDWI image (input)
### out_img: the filepath of the mask image to be created (output)
### thresh_val: the threshold value. For purposes of separating land from water, either use otsu or 0.2
def make_mask(in_img, out_img, thresh_val):
    print("--Starting Creating Mask Function--")

    if os.path.isfile(out_img):
        os.remove(out_img)

    
    ### Read in original raster
    with ds.open_dataset(in_img) as dataset:
        reader = ds.BasicReader(dataset, [1])
        in_aux = reader.aux_data
        in_coord_sys = reader.crs
        in_geocode = reader.geocoding
        
        in_raster = reader.read_raster(0, 0, reader.width, reader.height)
        print ("Finished Reading in raster")

    ### Set values to 1 or 0 to create mask 
    ndwi_band = in_raster.data[:, :, 0]
    print("Using Threshold Value of:", thresh_val)
    ndwi_band[ndwi_band <= thresh_val] = -32768
    ndwi_band[ndwi_band > thresh_val] = 1
    
    ### Convert to 16S 
    ndwi_band = ndwi_band.astype(np.int16)


    mask_raster = gobs.array_to_raster(ndwi_band)
    print("Finished Creating Mask")

    ### Write to File
    with ds.new_dataset(out_img, 'PCIDSK', '') as write_dataset:
            writer = ds.BasicWriter(write_dataset)
            
            writer.create(mask_raster)
            writer.write_raster(mask_raster)

            in_aux.set_file_metadata_value('NO_DATA_VALUE', '-32768')
            mdmap = [['BandDescription', 'Water_Mask'], ['NO_DATA_VALUE', '-32768']]
            in_aux.set_chan_metadata(mdmap, 1)
            in_aux.set_chan_description("Mask using value " + str(thresh_val), 1)
            writer.aux_data = in_aux


            writer.crs = in_coord_sys
            writer.geocoding = in_geocode



    print("--Creating Mask Function Finished--")
    return




#make_sieve: applies sieve function to first channel of image to remove small polygons based on certain threshold value
### in_img: the filepath of the image to be sieved(input)
### out_img: the filepath of the sieve_function image to be created (output)
### sval: the threshold value of the minimum area polygon to keep (usually set to 100 or 150)
def make_sieve(in_img, out_img, sval):
    print("--Starting Sieve Function--")
    if os.path.isfile(out_img):
        os.remove(out_img)

    ### Specify parameters    
    dbic = [1]
    
    with ds.open_dataset(in_img) as dset:
        num_chans = dset.aux_data.chan_count
        if (num_chans < 2):
            channel = [0, 1, 0, 0]
            pcimod(file=in_img, pciop="ADD", pcival=channel)

    dboc=[2]
    sthresh = [sval]
    keepvalu =[]
    connect=[8]
    
    ### Run Sieve Function
    sieve(in_img, dbic, dboc, sthresh, keepvalu, connect)
    print("Running Sieve Function")

    ### Save Channel to Separate File for ArcPro Processing
    with ds.open_dataset(in_img) as read_dataset, ds.new_dataset(out_img, 'PCIDSK', '') as write_dataset:
        reader = ds.BasicReader(read_dataset, [2])
        in_raster = reader.read_raster(0, 0, reader.width, reader.height)
        in_aux = reader.aux_data

        writer = ds.BasicWriter(write_dataset)
            
        writer.create(in_raster)
        writer.write_raster(in_raster)

        in_aux.set_file_metadata_value('NO_DATA_VALUE', '-32768')
        in_aux.set_chan_description('Water_Mask_Sieved', 1)
        writer.aux_data = in_aux


        writer.crs = reader.crs
        writer.geocoding = reader.geocoding

    ### Deleting interim channels
    pcimod(file = in_img, pciop = "DEL", pcival=[2])

    print("--Finished Sieve Function--")
    return



# make_bandratio: Taking in multispectral image, calculates a band_ratio of (B > NIR) & (B > Red)
### in_img: the filepath of the multispectral image (input)
### out_img: the filepath of the band ratio applied image to be created(output)
### bluechan: the number of the channel in input image that is of "Blue" wavelength
### nirchan: the number of the channel in input image that is of "NIR" wavelength
### redchan: the number of the channel in input image that is of "Red" wavelength
def make_bandratio(in_img, out_img, bluechan, nirchan, redchan):
    print("--Starting Band Ratio Function--")

    if os.path.isfile(out_img):
        os.remove(out_img)
    

    ### Read in original raster 
    with ds.open_dataset(in_img) as dataset:
        #Find Blue + NIR 
        reader = ds.BasicReader(dataset, [bluechan, nirchan, redchan])
        in_aux = reader.aux_data
        in_coord_sys = reader.crs
        in_geocode = reader.geocoding
        
        in_raster = reader.read_raster(0, 0, reader.width, reader.height)
        print ("Finished Reading in raster")

    ### Give the bands an appropriate name 
    blue_band = in_raster.data[:, :, 0]
    nir_band = in_raster.data[:, :, 1]
    red_band = in_raster.data[:, :, 2]


    ### Calculate Ratio (Blue > NIR ) & ( Blue > Red )
    band_ratio_result = (blue_band > nir_band)*(blue_band > red_band)
    band_ratio_result = band_ratio_result.astype(np.int16)
    band_ratio_result[band_ratio_result == 0] = -32768

    band_ratio_raster = gobs.array_to_raster(band_ratio_result)
    
    print("Finished Calculating Band Ratio")

    ### Output
    with ds.new_dataset(out_img, 'PCIDSK', '') as write_dataset:
        writer = ds.BasicWriter(write_dataset)

        writer.create(band_ratio_raster)
        writer.write_raster(band_ratio_raster)

    
        in_aux.set_file_metadata_value('NO_DATA_VALUE', '-32768')
        mdmap = [['BandDescription', 'NDWI'], ['NO_DATA_VALUE', '-32768']]
        in_aux.set_chan_metadata(mdmap, 1)
        in_aux.set_chan_description('BandRatio', 1)
        writer.aux_data = in_aux
        
        writer.crs = in_coord_sys
        writer.geocoding = in_geocode

    print("--Band Ratio Function Finished--")
    return

#make_addras: Taking in two mask images with only one channel each, combines their values using the OR function
#             (so if either image is True, the output mask will be true)
#             As well, will take metadata of first image forreference 
### in_img1: the filepath of the first mask image (input)
### in_img2: the filepath of the second mask image (input)
### out_img: the filepath of the combined mask image (output)
def make_or(in_img1, in_img2, out_img):
    print("--Starting OR Function--")

    if os.path.isfile(out_img):
        os.remove(out_img)
    
    ### Read in original raster 
    with ds.open_dataset(in_img1) as dataset1, ds.open_dataset(in_img2) as dataset2:
    
        reader = ds.BasicReader(dataset1, [1])
        reader2 = ds.BasicReader(dataset2, [1])
        in_aux = reader.aux_data
        in_coord_sys = reader.crs
        in_geocode = reader.geocoding
        
        ras1 = reader.read_raster(0, 0, reader.width, reader.height)
        ras2 = reader2.read_raster(0, 0, reader2.width, reader2.height)
        print ("Finished Reading in rasters")

    ### Apply OR Function
    outarr = (ras1.data[:, :, 0] == 1) | (ras2.data[:, :, 0] == 1)
    outarr = outarr.astype(np.int16)
    outarr[outarr == 0] = -32768

    out_raster = gobs.array_to_raster(outarr)
    print("Finished Applying OR")
    
    ### Output
    with ds.new_dataset(out_img, 'PCIDSK', '') as write_dataset:
        writer = ds.BasicWriter(write_dataset)

        writer.create(out_raster)
        writer.write_raster(out_raster)

    
        in_aux.set_file_metadata_value('NO_DATA_VALUE', '-32768')
        mdmap = [['BandDescription', 'NDWI'], ['NO_DATA_VALUE', '-32768']]
        in_aux.set_chan_metadata(mdmap, 1)
        in_aux.set_chan_description('FinalMask', 1)
        writer.aux_data = in_aux
        
        writer.crs = in_coord_sys
        writer.geocoding = in_geocode

    print("--OR Function Finished--")
    
def make_speclass(img_in, img_out, dem_in):
    print("--Starting SPECLASS Function--")
    start = time.time()
    speclass(img_in, dem_in, [], img_out, "PIX", "")
    end = time.time()
    time_elapsed = end - start
    print("Time Elapsed:")
    print(time_elapsed)
    print("--SPECLASS Function Finished--")
    
    
#find_otsu: given an image with one channel, determines threshold for two classes using otsu's method
#           the number of bins is automatically generated
### img_in: filepath of img (input)
### binnum: number of bins for histogram
def find_otsu(img_in):
    print("--Starting Calculating Otsu Value--")

    ### Read in raster
    with ds.open_dataset(img_in) as dataset:
        reader = ds.BasicReader(dataset, [1])
        in_aux = reader.aux_data
        in_coord_sys = reader.crs
        in_geocode = reader.geocoding
        
        in_raster = reader.read_raster(0, 0, reader.width, reader.height)
        print ("Finished Reading in raster")

    raster = in_raster.data[:, :, 0]
    raster[raster == -32768] = np.nan
    
    hist, bin_edges = np.histogram(raster[np.isfinite(raster)], bins = 'auto')
    print("Finished making histograms")
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    
    mean1 = np.cumsum(hist * bin_mids) / weight1
    mean2 = (np.cumsum((hist * bin_mids)[::-1]) / weight2[::-1])[::-1]

    inter_class_variance = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    index_of_max_val = np.argmax(inter_class_variance)

    otsu = bin_mids[:-1][index_of_max_val]
    print("--OTSU Function Finished--")
    return otsu


#find_channel: given an image, returns the number of a specific channel by comparing with channel metadata 'WavelengthName' and Channel Description
#              (in edge cases where multiple channels have same keywoard in channel description, function will return the smallest channel number)
#              If no channel is found, returns -1
### in_img: the filepath of the image (input)
### channame: the keyword wavelength name or channel description to look for
def find_channel(in_img, channame):   
    with ds.open_dataset(in_img) as dataset:
        aux = dataset.aux_data
        numchans = aux.chan_count

        for i in range(1, numchans + 1):
            bandname = aux.get_chan_metadata_value('WavelengthName', i)
            banddesc = aux.get_chan_description(i)
            if ((bandname == channame) or (channame in banddesc)):
                return i
    return -1    


################################### USER INPUT RELATED FUNCTIONS ############################################################################

def run_ndwi(img_filepath, ndwi_filepath):
    ### Find NDWI Channel Names
    greenchanname = "Coastal"
    greenchan = find_channel(img_filepath, greenchanname)
    if (greenchan == -1):
        greenchanname = "Green"
        greenchan = find_channel(img_filepath, greenchanname)
    nirchanname = "NIR"
    nirchan = find_channel(img_filepath, nirchanname)
    if (nirchan == -1):
        nirchanname = "Near-IR1"
        nirchan = find_channel(img_filepath, nirchanname)
    if (nirchan == -1):
        nirchanname = "Near-IR2"
        nirchan = find_channel(img_filepath, nirchanname)

    if (nirchan == -1 or greenchan == -1):
        print("Channels not found")
        return -1

    print("Running NDWI\n", ndwi_filepath, "\n")
    print(greenchanname + " " + nirchanname)

    make_ndwi(img_filepath, ndwi_filepath, greenchan, nirchan) 

def run_watermask(ndwi_filepath, watermask_filepath):
    otsu_val = find_otsu(ndwi_filepath)
    make_mask(ndwi_filepath, watermask_filepath, otsu_val)

def run_sieve(watermask_filepath, watermasksieve_filepath, thresh):
    #thresh = int(input("Please enter threshold area (smallest area you want to keep) for sieve function\n"))
    make_sieve(watermask_filepath, watermasksieve_filepath, thresh)

def run_bandratio(img_filepath, bandratio_filepath):
    bluechanname = "Blue"
    bluechan = find_channel(img_filepath, bluechanname)
    redchanname = "Red"
    redchan = find_channel(img_filepath, redchanname)
    nirchanname = "NIR"
    nirchan = find_channel(img_filepath, nirchanname)
    if (nirchan == -1):
        nirchanname = "Near-IR1"
        nirchan = find_channel(img_filepath, nirchanname)
    if (nirchan == -1):
        nirchanname = "Near-IR2"
        nirchan = find_channel(img_filepath, nirchanname)

    if (nirchan == -1 or bluechan == -1 or redchan == -1):
        print("Channels not found")
        return -1
    
    make_bandratio(img_filepath, bandratio_filepath, bluechan, nirchan, redchan)

def run_finalfilter(watermasksieve_filepath, bandratio_filepath, finalfilter_filepath):
    make_or(watermasksieve_filepath, bandratio_filepath, finalfilter_filepath)

def run_speclass(img_in, img_out):
    print("Select DEM file")
    dem_in = filedialog.askopenfilename(
            filetypes = (
                ("All Files", "*.*"),
                ("All Files", ".")
                )
            )
    make_speclass(img_in, img_out, dem_in)

    
#console_prompt: when run, starts process to interact with user on console 
def console_prompt():
    root = tk.Tk()
    root.withdraw()

    choice = int(input("Please choose one of the following options by entering a corresponding number:\n"
                   + "1. Process all img files in a folder\n"
                   + "2. Process a single file\n"
                   + "3. Run a specific function on a single file\n"))

    ### Processing a folder
    if (choice == 1):
        dir_filepath = filedialog.askdirectory()
        item_list = os.listdir(dir_filepath)
        thresh = int(input("Please enter a threshold value for the sieve function (the smallest area you want to keep)\n"))

        for item in item_list:
            item_filepath = os.path.join(dir_filepath, item)
            if(os.path.isfile(item_filepath) and item_filepath.endswith(".pix")):
                process_file(item_filepath, thresh)
        
    ### Processing single file
    elif (choice == 2):
        print("Please choose a file to process")
        img_filepath = filedialog.askopenfilename(
            filetypes = (
                ("PCI Files", "*.pix"),
                ("All Files", ".")
                )
            )
        
        thresh = int(input("Please enter a threshold value for the sieve function (the smallest area you want to keep)\n"))
        process_file(img_filepath, thresh)

    elif (choice == 3):
        fun_num = int(input("Please select which function:\n"
                           + "1. NDWI only\n"
                           + "2. WaterMask only\n"
                           + "3. Sieve only\n"
                           + "4. Bandratio only\n"
                           + "5. Apply OR function to two raster files\n"
                           + "6. Run SPECLASS Function\n"))
        process_fun(fun_num)
        
#process_fun: runs 1 function on a file
### fun_num: the number dictating which function
def process_fun(fun_num):
    print("Please choose a file to process")
    img_filepath = filedialog.askopenfilename(
            filetypes = (
                ("PCI Files", "*.pix"),
                ("All Files", ".")
                )
            )
    img_filepath2 = ""

    if (fun_num == 5):
        img_filepath2 = filedialog.askopenfilename(
            filetypes = (
                ("PCI Files", "*.pix"),
                ("All Files", ".")
                )
            )
    suffix = input("Please choose a suffix for the output\n")
    working_dir = os.path.dirname(os.path.abspath(img_filepath))
    img_basename, img_extension = os.path.splitext(os.path.basename(img_filepath))
    output_filepath = os.path.join(working_dir, img_basename + suffix + img_extension)
        
    if (fun_num == 1):
        run_ndwi(img_filepath, output_filepath)

    elif (fun_num == 2):
        run_watermask(img_filepath, output_filepath)

    elif (fun_num == 3):
        thresh = int(input("Please enter a threshold value for the sieve function (the smallest area you want to keep)\n"))
        run_sieve(img_filepath, output_filepath, thresh)

    elif (fun_num == 4):
        run_bandratio(img_filepath, output_filepath)
        
    elif (fun_num == 5):
        run_finalfilter(img_filepath, img_filepath2, output_filepath)

    elif (fun_num == 6):
        run_speclass(img_filepath, output_filepath)
        

        
#process_file: runs all function on one file
def process_file(img_filepath, thresh):
    print("\n\n" + img_filepath)
    ### Get Working directory and basename
    working_dir = os.path.dirname(os.path.abspath(img_filepath))
    img_basename, img_extension = os.path.splitext(os.path.basename(img_filepath))
    working_dir = os.path.join(working_dir, img_basename)
    
    if (os.path.exists(working_dir) == False):
        os.makedirs(working_dir)
    
    ### Filepaths for outputs 
    ndwi_filepath = os.path.join(working_dir, img_basename + '_NDWI' + img_extension)
    watermask_filepath = os.path.join(working_dir, img_basename + '_Mask' + img_extension)
    watermasksieve_filepath = os.path.join(working_dir, img_basename + '_Sieve' + img_extension)
    bandratio_filepath = os.path.join(working_dir, img_basename + '_BandRatio' + img_extension)
    finalfilter_filepath = os.path.join(working_dir, img_basename + '_Final' + img_extension)

    ### Running NDWI Function
    run_ndwi(img_filepath, ndwi_filepath)

    ### Running Water Mask Function
    run_watermask(ndwi_filepath, watermask_filepath)
    
    ### Running Sieve Function
    run_sieve(watermask_filepath, watermasksieve_filepath, thresh)

    ### Running Band Ratio Function
    run_bandratio(img_filepath, bandratio_filepath)

    ### Make Final Filter
    run_finalfilter(watermasksieve_filepath, bandratio_filepath, finalfilter_filepath)

    ### Move file to path
    #os.rename(img_filepath, os.path.join(working_dir, img_basename + img_extension))

#################### MAIN FUNCTION ###############################
console_prompt()
