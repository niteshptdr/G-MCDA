# -*- coding: utf-8 -*-

shp_file_path = r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\shps\gy_doab.shp"
product_id ="Soil texture USDA"
out_dir = r"D:\Mar\MAR\Uploads\npt4@gmail.com\Ganga Yamuna 2\Downloadtiff"


import geemap
import ee
import os
import geopandas as gpd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling


#ee.Authenticate()
ee.Initialize()

#shp_file_path = r"D:\gmcda\testing\data\shps\hindon.shp"
# out_dir = os.path.join(os.path.expanduser('~'), 'Downloads')
#out_dir = r'D:\gmcda\testing\data'
max_pixels = 1000000
deg_to_meter = 111000


product_IDs = ['Land Cover Copernicus 100m', 'Precipitation CHIRPS 5km', 
               'rainydays CHIRPS 5km',
               'Slope srtm90m', 'Population GPWv411', 'TWS GRACE Monthly Mass Grids - Land',
               'Soil texture USDA', 'Significant Sens slope of ppt - CMIP5 RCP4.5',
               'Significant Sens slope of rainydays - CMIP5 RCP4.5',
               'Distance from river - Hydroshed River-based']

#product_id = 'Distance from river - Hydroshed River-based'


# MKT at GEE...
def mkt_gee(image_collection):
    
    coll = image_collection

    afterFilter = ee.Filter.lessThan(**{
      'leftField': 'system:time_start',
      'rightField': 'system:time_start'
    })

    joined = ee.ImageCollection(ee.Join.saveAll('after').apply(**{
        'primary': coll, 'secondary': coll,'condition': afterFilter
    }))


    # i and j are images
    def sign(i, j): 
        return ee.Image(j).neq(i).multiply(ee.Image(j).subtract(i).clamp(-1, 1)).int()

    def get_kendall_img(current):
        afterCollection = ee.ImageCollection.fromImages(current.get('after'))

        def map_to_sign(image):
            sign_img = sign(current, image)
            return ee.Image(sign_img).unmask(0)

        return afterCollection.map(map_to_sign)

    kendall = ee.ImageCollection(joined.map(get_kendall_img).flatten()).reduce('sum', 2) 

    # palette = ['red', 'white', 'green'];
    #Map.addLayer(kendall, {'min': -10, 'max': 10, 'palette': palette}, 'kendall');

    #Map


    def slope(i, j):  # i and j are images:
          return ee.Image(j).subtract(i) \
              .divide(ee.Image(j).date().difference(ee.Image(i).date(), 'days')) \
              .rename('slope') \
              .float()

    def func_qhx(current):
        afterCollection = ee.ImageCollection.fromImages(current.get('after'))
        def func_slope(image):
            return ee.Image(slope(current, image))
        return afterCollection.map(func_slope)

    slopes = ee.ImageCollection(joined.map(func_qhx).flatten())

    sensSlope = slopes.reduce(ee.Reducer.median(), 2); # Set parallelScale.
    #Map.addLayer(sensSlope, {'palette': palette}, 'sensSlope')

    epochDate = ee.Date('1970-01-01')

    def func_oax(image):
        epochDays = image.date().difference(epochDate, 'days').float()
        return image.subtract(sensSlope.multiply(epochDays)).float()

    # sensIntercept = coll.map(func_oax).reduce(ee.Reducer.median(), 2)

    #Map.addLayer(sensIntercept, {}, 'sensIntercept')
    #Map


    # Values that are in a group (ties).  Set all else to zero.

    def func_zbj(i):
        def i_eq(j):
            return i.eq(j)
        matches = coll.map(i_eq).sum()
        return i.multiply(matches.gt(1))

    groups = coll.map(func_zbj)

    # Compute tie group sizes in a sequence.  The first group is discarded.
    def group(array):
        length = array.arrayLength(0)
      # Array of indices.  These are 1-indexed.
        indices = ee.Image([1]) \
          .arrayRepeat(0, length) \
          .arrayAccum(0, ee.Reducer.sum()) \
          .toArray(1)
        sorted = array.arraySort()
        left = sorted.arraySlice(0, 1)
        right = sorted.arraySlice(0, 0, -1)
        # Indices of the end of runs.
        mask = left.neq(right) \
          .arrayCat(ee.Image(ee.Array([[1]])), 0)
        runIndices = indices.arrayMask(mask)
      # Subtract the indices to get run lengths.
        groupSizes = runIndices.arraySlice(0, 1) \
          .subtract(runIndices.arraySlice(0, 0, -1))
        return groupSizes

    # See equation 2.6 in Sen (1968).
    def factors(image):
        return image.expression('b() * (b() - 1) * (b() * 2 + 5)')

    groupSizes = group(groups.toArray())
    groupFactors = factors(groupSizes)
    groupFactorSum = groupFactors.arrayReduce('sum', [0]) \
          .arrayGet([0, 0])

    count = joined.count()

    kendallVariance = factors(count) \
        .subtract(groupFactorSum) \
        .divide(18) \
        .float()
    #Map.addLayer(kendallVariance, {}, 'kendallVariance')



    # Compute Z-statistics.
    zero = kendall.multiply(kendall.eq(0))
    pos = kendall.multiply(kendall.gt(0)).subtract(1)
    neg = kendall.multiply(kendall.lt(0)).add(1)

    z = zero \
        .add(pos.divide(kendallVariance.sqrt())) \
        .add(neg.divide(kendallVariance.sqrt()))
    #Map.addLayer(z, {'min': -2, 'max': 2}, 'z')

    # https:#en.wikipedia.Org/wiki/Error_function#Cumulative_distribution_function
    def eeCdf(z):
        return ee.Image(0.5) \
          .multiply(ee.Image(1).add(ee.Image(z).divide(ee.Image(2).sqrt()).erf()))

    def invCdf(p):
        return ee.Image(2).sqrt() \
          .multiply(ee.Image(p).multiply(2).subtract(1).erfInv())

    # Compute P-values.
    p = ee.Image(1).subtract(eeCdf(z.abs()))
    #Map.addLayer(p, {'min': 0, 'max': 1}, 'p')

    # Pixels that can have the None hypothesis (there is no trend) rejected.
    # Specifically, if the True trend is zero, there would be less than 5%
    # chance of randomly obtaining the observed result (that there is a trend).
    #Map.addLayer(p.lte(0.025), {'min': 0, 'max': 1}, 'significant trends')

    significant_sensSlope = sensSlope.updateMask(p.lte(0.025))
    #Map.addLayer(significant_sensSlope, {'min': 0, 'max': 1}, 'significant slopes')
    #Map
    return significant_sensSlope



def get_gee_layer(product_id, shp_file_path, out_dir, max_pixels, deg_to_meter):
    
    product_meta = {
        'Land Cover Copernicus 100m':{'name': 'Land Cover Copernicus 100m',
                                       'tif_file': 'lc.tif'},
    
        'Precipitation CHIRPS 5km':{'name': 'Precipitation CHIRPS 5km',
                                       'tif_file': 'ppt.tif'},

        'rainydays CHIRPS 5km':{'name': 'rainydays CHIRPS 5km',
                                       'tif_file': 'rainydays.tif'},    
    
    
        'Slope srtm90m': {'name': 'Slope srtm90m',
                          'tif_file': 'slope.tif'},
        
        'Population GPWv411': {'name': 'Population GPWv411',
                          'tif_file': 'population.tif'},
            
        'TWS GRACE Monthly Mass Grids - Land': {'name': 'TWS GRACE Monthly Mass Grids',
                          'tif_file': 'tws.tif'},
            
        'Soil texture USDA': {'name': 'Soil texture USDA',
                          'tif_file': 'soil.tif'},

        'Significant Sens slope of ppt - CMIP5 RCP4.5': {'name': 'Significant Sens slope of ppt - CMIP5 RCP4.5',
                          'tif_file': 'ppt_sensslope.tif'},    
        
        'Significant Sens slope of rainydays - CMIP5 RCP4.5': {'name': 'Significant Sens slope of rainydays - CMIP5 RCP4.5',
                          'tif_file': 'rainydays_sensslope.tif'},  
        
        'Distance from river - Hydroshed River-based': {'name': 'Distance from river - Hydroshed River-based',
                          'tif_file': 'distance.tif'},  
        }
    
    aoi_object = geemap.shp_to_ee(shp_file_path)
    gpd_aoi = gpd.read_file(shp_file_path)
    aoi_bounds = gpd_aoi.geometry.bounds
    
    minx = aoi_bounds.minx[0]; maxx = aoi_bounds.maxx[0]; miny = aoi_bounds.miny[0]; maxy = aoi_bounds.maxy[0]; 
    
    
    x_bound_len = (maxx-minx) * deg_to_meter
    y_bound_len = (maxy-miny) * deg_to_meter
    
    
    resolution = 50 # m
    while (y_bound_len/resolution  * x_bound_len/resolution) > max_pixels:
        resolution = resolution + 50
    
    
    if product_id == 'Land Cover Copernicus 100m':
        image = ee.ImageCollection('COPERNICUS/Landcover/100m/Proba-V-C3/Global').first()
        image = image.select('discrete_classification')
        
    elif product_id == 'Precipitation CHIRPS 5km':
        image = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate('2011-01-01','2020-12-31')
        startyear = 2010
        endyear = 2020
        years = ee.List.sequence(startyear, endyear)
        
        def get_annual_chirps(year):
            chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
            annual = chirps.filter(ee.Filter.calendarRange(year, year, 'year')).sum();
            return annual.set('year', year).set('system:time_start', ee.Date.fromYMD(year, 1, 1));
        
        annual_chirps = years.map(get_annual_chirps)
        image = ee.ImageCollection.fromImages(annual_chirps).mean()
        
    elif product_id == 'rainydays CHIRPS 5km':
        image = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate('2011-01-01','2020-12-31')
        startyear = 2010
        endyear = 2020
        years = ee.List.sequence(startyear, endyear)
        
        def get_annual_chirps(year):
            chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
            
            def count_rainyday(image):
                return image.gt(2.5)

            annual = chirps.filter(ee.Filter.calendarRange(year, year, 'year')).map(count_rainyday).sum() # IMD rainyday
            
            return annual.set('year', year).set('system:time_start', ee.Date.fromYMD(year, 1, 1));
            
        annual_chirps = years.map(get_annual_chirps)
        image = ee.ImageCollection.fromImages(annual_chirps).mean().int()


    elif product_id == 'Slope srtm90m':
        elev = ee.Image("CGIAR/SRTM90_V4")
        image = ee.Terrain.slope(elev)
    elif product_id == 'Population GPWv411':
        population = ee.ImageCollection("CIESIN/GPWv411/GPW_Basic_Demographic_Characteristics").filterDate('2010-01-01', '2020-12-31')
        image = population.select('basic_demographic_characteristics').first()
        
    elif product_id == 'TWS GRACE Monthly Mass Grids - Land':
        dataset = ee.ImageCollection('NASA/GRACE/MASS_GRIDS/LAND').filter(ee.Filter.date('2016-05-01', '2016-05-31')).first()
        image = dataset.select('lwe_thickness_csr')
    
    elif product_id == 'Soil texture USDA':
        image = ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02');
        image = image.select('b10')
        
    elif product_id == 'Distance from river - Hydroshed River-based':
        rivers = ee.FeatureCollection("WWF/HydroSHEDS/v1/FreeFlowingRivers").filterBounds(aoi_object).filter('RIV_ORD <= 6')
        
        distance = rivers.distance(**{'searchRadius': 100000, 'maxError': 50}).clip(aoi_object)
        image = distance
        
    elif product_id == 'Significant Sens slope of ppt - CMIP5 RCP4.5':
        
        Date_Start = ee.Date('2021-01-01')
        Date_End = ee.Date('2040-12-31') # currently maximum 20 years can be downloaded
        
        def get_daily(in_date):
        
            daily = ee.ImageCollection("NASA/NEX-GDDP") \
                          .filterDate(ee.Date(in_date), ee.Date(in_date).advance(1, 'day')) \
                          .filter(ee.Filter.eq('scenario', 'rcp45')) \
                          .select('pr') \
                          .mean()
        
            return daily.set('system:time_start', ee.Date(in_date).millis())
        
        # Create list of dates for time series
        n_days = Date_End.difference(Date_Start,'day').round()
        dates = ee.List.sequence(0,n_days,1)
        
        def make_datelist(n):
          return Date_Start.advance(n,'day')
        
        dates = dates.map(make_datelist)
        
        daily_fut = dates.map(get_daily)
        # daily_fut = ee.ImageCollection(daily_fut)
        
        startyear = Date_Start.get('year')
        endyear = Date_End.get('year')
        years = ee.List.sequence(startyear, endyear)
        
        fut_ppt = ee.ImageCollection(daily_fut)
        
        def get_annual(year):
            annual = fut_ppt.filter(ee.Filter.calendarRange(year, year, 'year')) \
                                .sum() \
                                .clip(aoi_object) \
                                .multiply(86400)
            return annual.set('year', year).set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis())
        
        annual_ppt = years.map(get_annual)
        annual_ppt = ee.ImageCollection(annual_ppt)
        
        coll = annual_ppt
        
        significant_sens_slope = mkt_gee(coll)
        image = significant_sens_slope


    elif product_id == 'Significant Sens slope of rainydays - CMIP5 RCP4.5':

        Date_Start = ee.Date('2021-01-01')
        Date_End = ee.Date('2040-12-31') # currently maximum 20 years can be downloaded
        
        def get_daily(in_date):
        
            daily = ee.ImageCollection("NASA/NEX-GDDP") \
                          .filterDate(ee.Date(in_date), ee.Date(in_date).advance(1, 'day')) \
                          .filter(ee.Filter.eq('scenario', 'rcp45')) \
                          .select('pr') \
                          .mean() \
                          .gt(0.000028935185185185186) # rainy day rain>2.5 mm - IMD
        
            return daily.set('system:time_start', ee.Date(in_date).millis())
        
        # Create list of dates for time series
        n_days = Date_End.difference(Date_Start,'day').round()
        dates = ee.List.sequence(0,n_days,1)
        
        def make_datelist(n):
          return Date_Start.advance(n,'day')
        
        dates = dates.map(make_datelist)
        
        daily_fut = dates.map(get_daily)
        # daily_fut = ee.ImageCollection(daily_fut)
        
        startyear = Date_Start.get('year')
        endyear = Date_End.get('year')
        years = ee.List.sequence(startyear, endyear)
        
        fut_ppt = ee.ImageCollection(daily_fut)
        
        def get_annual(year):
            annual = fut_ppt.filter(ee.Filter.calendarRange(year, year, 'year')) \
                                .sum() \
                                .clip(aoi_object)
                                # .multiply(86400)
            return annual.set('year', year).set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis())
        
        annual_ppt = years.map(get_annual)
        annual_ppt = ee.ImageCollection(annual_ppt)
        
        coll = annual_ppt
        
        significant_sens_slope = mkt_gee(coll)
        image = significant_sens_slope



    
    
    temp_filename = os.path.join(out_dir, 'temp.tif')
    roi = aoi_object.geometry()
    image = image.clip(roi).unmask()
    geemap.ee_export_image(
        image, filename=temp_filename, scale=resolution, region=roi, file_per_band=False
    )
    
    
    # Check if the layer is in WGS-84, reproject if not
    dst_crs = 'EPSG:4326'
    src_file = temp_filename
    out_tif_file = product_meta[product_id]['tif_file']
    dst_file = os.path.join(out_dir, out_tif_file)
    
    # src_file  = r'D:\gmcda\testing\data\temp_recls_lc.tif'
    # dst_file = r'D:\gmcda\testing\data\temp_recls_lc_reproj.tif'
    
    
    with rasterio.open(src_file) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
    
        with rasterio.open(dst_file, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
    

get_gee_layer(product_id, shp_file_path, out_dir, max_pixels, deg_to_meter)







