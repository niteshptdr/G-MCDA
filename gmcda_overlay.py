# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 20:15:35 2022

@author: HP
"""

'''
1. convert to uniform resolution
2. perfrom weighted overlay
3. write clipped results
'''


import rasterio
from rasterio.enums import Resampling
import fiona
import rasterio.mask
import numpy as np
import os
import gdal


fd=r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif"
aoi_shp_path=r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\shps\gy_doab.shp"







file_dir = os.path.join(fd)
criteria_meta = {
             "soil":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_soil.tif","weight":0.140350877192982,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_soil.tif"},
             "slope":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_slope.tif","weight":0.157894736842105,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_slope.tif"},
             "aquifer_material":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_aquifer_material.tif","weight":0.12280701754386,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_aquifer_material.tif"},
             "lc":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_lc.tif","weight":0.140350877192982,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_lc.tif"},
             "gw_utilization_fraction":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_gw_utilization_fraction.tif","weight":0.0350877192982456,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_gw_utilization_fraction.tif"},
             "depth_gw_for_gmcda_2":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_depth_gw_for_gmcda_2.tif","weight":0.12280701754386,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_depth_gw_for_gmcda_2.tif"},
             "rainydays":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_rainydays.tif","weight":0.175438596491228,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_rainydays.tif"},
             "distance":  {"recl_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\reclasstif\recls_distance.tif","weight":0.105263157894737,"final_uni_file":r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\overlaytif\Uni_distance.tif"},
    
    
    
    
    
    
    
    
    
    
    
    
    
    }

def overlay_different_res(file_dir, aoi_shp_path, criteria_meta):
    
    """

    Parameters
    ----------
    file_dir : path of working directory where tif files are saved
        
    aoi_shp_path : SHAPE FILE PATH
        path of shape file of AOI
    criteria_meta : DICTIONARY
        contains previously saved TIFs and weights

    Returns
    -------
    Suitability Map
    
    Masks all pixels with value 0 (not suitable)
    """
    with fiona.open(aoi_shp_path) as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    
    
    # get the finest resolution amongst the layers
    criterion = criteria_meta.keys()
    resolutions_xy = []
    width_height = []
    for key in criteria_meta.keys():
        in_file_path = os.path.join(criteria_meta[key]['recl_file'])
        
        with rasterio.open(in_file_path) as src:
            resolutions_xy.append( [ src.profile['transform'][0],  abs(src.profile['transform'][4]) ] ) 
            width_height.append([src.profile['width'], src.profile['height']])
            # resolutions_y.append(abs(src.profile['transform'][4]))
    
            # print(key, in_file_path, src.profile['transform'][0])
    
    
    dict_resolutions = dict(zip(criterion, resolutions_xy))
    dict_width_height = dict(zip(criterion, width_height))
    
    # width = max(width_height)[0]
    # height = max(width_height)[1]
    
    resolution_x = max(resolutions_xy)[0]
    resolution_y = max(resolutions_xy)[1]
    
    des_ind = resolutions_xy.index(min(resolutions_xy))
    width = width_height[des_ind][0]
    height = width_height[des_ind][1]
    
    for key in criteria_meta.keys():
        in_file_path = os.path.join(criteria_meta[key]['recl_file'])
        with rasterio.open(in_file_path) as dataset:
             upscale_factor_x = dict_resolutions[key][0] / resolution_x
             upscale_factor_y = dict_resolutions[key][1] / resolution_y
    
             # resample data to target shape
             data = dataset.read(
                    out_shape=(
                                dataset.count,
                                # int(dataset.height * upscale_factor_y),
                                # int(dataset.width * upscale_factor_x)
                                height, # height
                                width, # width
                                
                                ),
                                resampling=Resampling.bilinear
                                )
            
             # scale image transform
             transform = dataset.transform * dataset.transform.scale(
                    (dataset.width / data.shape[-1]),
                    (dataset.height / data.shape[-2])
                )
             profile = dataset.profile
             # height =  int(dataset.height * upscale_factor_y)
             # width =   int(dataset.width * upscale_factor_x)
    
             
             profile.update(transform=transform, driver='GTiff', height=height, width=width, crs=dataset.crs)
             temp_out_file = os.path.join(file_dir, 'temp.tif')
             with rasterio.open(temp_out_file, 'w', **profile) as dst:
                 dst.write(data)
    
    
    
            # Clip to AOI
             with rasterio.open(temp_out_file) as src:
    
                array, out_transform = rasterio.mask.mask(src, shapes, crop=True)
    
                profile = src.profile
                # out_meta = src.meta
                profile.update(
                { "height": array.shape[1],
                  "width": array.shape[2],
                  "transform": out_transform
                 # # 'nodata': 0
                 }  )
             temp_out_unifile = os.path.join(file_dir, 'uni_temp.tif')
             with rasterio.open(temp_out_unifile, 'w', **profile) as dst:
                dst.write(array*100)
            
             out_unifile = os.path.join(criteria_meta[key]['final_uni_file'])
             ds = gdal.Open(temp_out_unifile)
             band = ds.GetRasterBand(1)
             arr = band.ReadAsArray()
             [rows, cols] = arr.shape
             arr_min = arr.min()
             arr_max = arr.max()
             arr_mean = int(arr.mean())
             #arr_out = np.where((arr < arr_mean), 10000, arr)
             arr_out = arr
             driver = gdal.GetDriverByName("GTiff")
             
             outdata = driver.Create(out_unifile, cols, rows, 1, gdal.GDT_UInt16)
             outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
             outdata.SetProjection(ds.GetProjection())##sets same projection as input
             outdata.GetRasterBand(1).WriteArray(arr_out)
             outdata.GetRasterBand(1).SetNoDataValue(0)##if you want these values transparent
             outdata.FlushCache() ##saves to disk!!
             outdata = None
             band=None
             ds=None
            
             os.remove(temp_out_unifile)
                
    # Perfrom Boolean & weighted overlay
    wt_sum = np.zeros((profile['height'], profile['width']))
    
    for key in criteria_meta.keys():
        
        filename = os.path.join(file_dir, criteria_meta[key]['final_uni_file'])
        with rasterio.open(filename) as src:
            array = src.read(1)
            profile = src.profile
            
            
            image_read_masked = np.ma.masked_array(array, mask=(array == 0))
            wt_sum = wt_sum + image_read_masked * criteria_meta[key]['weight']
    
    profile.update( {'nodata': 0} )
    suit_temp = os.path.join(file_dir, 'suit_temp.tif')
    wt_sum = np.reshape(wt_sum, (1,profile['height'], profile['width']))
    wt_sum = wt_sum
    wt_sum = wt_sum.astype(rasterio.uint8)
    with rasterio.Env():
        profile.update(
            dtype = rasterio.uint8,
            compress='lzw',
            )
        
        with rasterio.open(suit_temp, 'w', **profile) as dst:
            dst.write(wt_sum)
    
    
    
    ds = gdal.Open(suit_temp)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    [rows, cols] = arr.shape
    arr_min = arr.min()
    arr_max = arr.max()
    arr_mean = int(arr.mean())
    arr_out = arr #np.where((arr < arr_mean), 10000, arr)
    driver = gdal.GetDriverByName("GTiff")
    out_file = os.path.join(file_dir, 'suitability.tif')
    outdata = driver.Create(out_file, cols, rows, 1, gdal.GDT_UInt16)
    outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    outdata.SetProjection(ds.GetProjection())##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(arr_out)
    outdata.GetRasterBand(1).SetNoDataValue(0)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!
    outdata = None
    band=None
    ds=None
    
    os.remove(suit_temp)
    
overlay_different_res(file_dir, aoi_shp_path, criteria_meta)
