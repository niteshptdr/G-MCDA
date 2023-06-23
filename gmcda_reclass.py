# -*- coding: utf-8 -*-
"""
Standardize the thematic layers
"""

import rasterio
import numpy as np
import rasterio.mask
import fiona


in_file_path= r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\Downloadtiff\lc.tif"
out_file=r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_lc.tif"
recls_file=r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\textfiles\lc.txt"
aoi_shp=r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\shps\gy_doab.shp"
temp_out_file= r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\temp_lc.tif"



def reclassify(in_file_path, out_file, temp_out_file, recls_file, aoi_shp):
    # Read AOI
    with fiona.open(aoi_shp) as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    
    # Create dictionary from text file content
    with open(recls_file) as f:
        lines = f.readlines()
    
    dict_classes = {}
    dict_function_val = {}
    if len(lines[0].split()) > 2: # Continues data layer 
        layer_type = 'continues'
        line_no = 1
        for line in lines:
            content = line.split()
            if len(content) == 3:
                dict_classes[line_no] = [float(content[0]), float(content[1]), float(content[2])]
            elif len(content) == 4:
                dict_function_val[line_no] = [float(content[0]), float(content[1]), float(content[2]), float(content[3])]
            line_no= line_no + 1
    else: # Thematic data layer
        layer_type = 'thematic'
        for line in lines:
            content = line.split()
            dict_classes[content[0]] = float(content[1])
    
    
    with rasterio.open(in_file_path) as src:
    
        # array, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        
        array = src.read()
        profile = src.profile
        # out_meta = src.meta
        profile.update(
        {"dtype": "float32",
         # "height": array.shape[1],
         # "width": array.shape[2],
         # "transform": out_transform
         # # 'nodata': 0
         })
        array = np.float32(array)
        dict_indices = {}
        
        if layer_type == 'continues':
            for key in dict_classes.keys():
                dict_indices[key] = np.where((array >= dict_classes[key][0]) & (array < dict_classes[key][1]))
    
            for key in dict_function_val.keys():
                ymin = dict_function_val[key][2]
                ymax = dict_function_val[key][3]
                xmin = dict_function_val[key][0]
                xmax = dict_function_val[key][1]
                
                gradient = (ymax-ymin)/(xmax-xmin)
                
                bands, rows, cols = array.shape
                
                for r in range(rows):
                    for c in range(cols):
                        if (array[0,r,c] >= xmin) & (array[0,r,c] < xmax):
                            x =abs(array[0,r,c]-xmin)
                            array[0,r,c] = gradient*x + ymin
    
            for key in dict_indices.keys():
                array[dict_indices[key]] = dict_classes[key][2]     
        else:
            for key in dict_classes.keys():
                dict_indices[key] = np.where(array == int(key))
            
            for key in dict_indices.keys():
                array[dict_indices[key]] = dict_classes[key]
    
    ind_pix_not_valid = np.where((array > 1) | (array < 0))
    array[ind_pix_not_valid] = 0
    with rasterio.open(temp_out_file, 'w', **profile) as dst:
        dst.write(array)
    
    # Clip & Re-write...
    with rasterio.open(temp_out_file) as src:
    
        array, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        
        profile = src.profile
        # out_meta = src.meta
        profile.update(
        {"dtype": "float32",
          "height": array.shape[1],
          "width": array.shape[2],
          "transform": out_transform
         # # 'nodata': 0
         }  )
    with rasterio.open(out_file, 'w', **profile) as dst:
        dst.write(array)
        
        
reclassify(in_file_path, out_file, temp_out_file, recls_file, aoi_shp)
