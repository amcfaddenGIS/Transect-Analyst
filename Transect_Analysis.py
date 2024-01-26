# Test bed for Thesis Tools
import os
os.environ["PROJ_LIB"] = "C:/Users/alexm/anaconda3/envs/PyQTApp/Library/share/proj"
print(os.environ["PROJ_LIB"])
import json
import glob
import math
import tempfile
from PIL import Image, ImageDraw
import shutil
import pyproj
import shapely.geometry
from jinja2 import Template
import sys
from osgeo import gdal, ogr, osr
import geopandas as gpd
from rasterio import mask, sample, plot
import matplotlib.pyplot as plt
import time
import pandas as pd
import re
import fiona
import rasterio
from rasterio import warp
import numpy as np
import folium
import pathlib as ptlib


def List_Rasters(Path, file_type="tif"):
    """
    Creates list of raster file paths

    :param Path: Directory containing a series of rasters
    :param file_type: File type for listing
    :return: List of all raster file paths
    """
    path = ptlib.Path(Path)
    All_Rasters = list(path.glob('**/*.{}'.format(file_type)))
    Raster_List = []
    for raster in All_Rasters:
        Raster_Path = str(raster).replace("/", "\\")
        Raster_List.append(Raster_Path)
    return Raster_List

def List_Geojsons(Path, file_type="geojson"):
    path = ptlib.Path(Path)
    geojsons = list(path.glob('**/*.{}'.format(file_type)))
    geojson_list = []
    for g in geojsons:
        geojson = str(g).replace("/", "\\")
        geojson_list.append(geojson)
    return geojson_list


def reproject_raster(raster_path, output_location):
    """
    Reproject rasters to CRS 4326 so they can be viewed in folium


    :param raster_path: Path for the raster. Must be a tif
    :param output_location: Output location for the reprojected raster
    :return: Reprojected raster
    """
    for ra in raster_path:
        file_name = os.path.basename(ra)
        # Grab the name of the mosaic from the list of raster files
        image_n = os.path.splitext(file_name)[0]
        with rasterio.open(ra) as r:
            out_meta = r.meta
            width, height = r.width, r.height  # 1182, 1050
            # Next get the transform information
            transform = r.transform
            # Get the coordinate reference system
            crs = r.crs
            # Load band
            array = r.read(1)
        # Use the WGS84 projection
        dst_crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs"
        # Calculate the default transform
        dst_transform, dst_width, dst_height = warp.calculate_default_transform(crs, dst_crs, width, height, *r.bounds)
        # Reproject the data to the desired CRS
        dst_data = np.zeros((dst_height, dst_width), np.float32)
        rasterio.warp.reproject(
            source=array,
            destination=dst_data,
            src_transform=transform,
            src_crs=crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=rasterio.warp.Resampling.nearest
        )
        # Save the output raster
        with rasterio.open("{}/Reproject_{}.tif".format(output_location, image_n), "w", driver="GTiff", height=dst_height, width=dst_width, count=1,
                           dtype=np.float32, crs=dst_crs, transform=dst_transform) as dst:
            dst.write(dst_data, 1)

class Sample_Transect():
    """
    The Sample Transect class is used for sampling temperatures along drawn transects from folium. Sampling transects must be initialized with
    an output location for method results, the transect itself, a conversion crs for transect conversion and the identified raster for sampling.

    """
    def __init__(self, output_location, line, conversion_crs, raster):
        """

        :param output_location: Output location for all class methods
        :param line: Identified transect
        :param conversion_crs: CRS that the transect will be converted to
        :param raster: Raster used for sampling
        """
        self.output_location = output_location
        self.line = line
        self.conversion_crs = conversion_crs
        self.raster = raster
        self.centroids = None
        self.buffer_raster = None
        self.buffer = None

    def Transform_Line_to_UTM(self):
        """
        Changes the CRS of the drawn transect from 4326 to an identified CRS by the user. CRS is meant to overlap with the crs of the
        raster being sampled

        :return: Shapefile of the reprojected transect
        """
        # Open the data source and load up the layer and layer definition (Field names, etc.)
        source = ogr.Open(self.line)
        layer = source.GetLayer()
        layer_defn = layer.GetLayerDefn()
        spatialRef = layer.GetSpatialRef()
        # Create a spatial reference object and import the UTM 11N EPSG
        UTM = osr.SpatialReference()
        UTM.ImportFromEPSG(self.conversion_crs)
        # Create a spatial reference object and import the WGS 84 Web Mercator spatial reference
        base_epsg = osr.SpatialReference()
        layer_epsg = 4326
        base_epsg.ImportFromEPSG(layer_epsg)
        # Create a transformation function to transform points and create a new line string with the coords
        transform = osr.CoordinateTransformation(base_epsg, UTM)
        # Loop through features
        # Grab the first feature in the layer
        base_feature = layer.GetNextFeature()
        x0 = base_feature.GetGeometryRef().GetX(0)
        y0 = base_feature.GetGeometryRef().GetY(0)
        # Output name for the layer, based on the title field
        name = base_feature.GetField(1)
        # Transform coordinates from WGS to UTM
        points = [
                transform.TransformPoint(base_feature.GetGeometryRef().GetY(0), base_feature.GetGeometryRef().GetX(0)),
                transform.TransformPoint(base_feature.GetGeometryRef().GetY(1), base_feature.GetGeometryRef().GetX(1))]
            # Create a new feature and add points to that feature with the transformed points
            # When creating geometry using OGR, you need to specify the type of geometry
        new_geom = ogr.Geometry(type=ogr.wkbLineString)
        new_geom.AddPoint(points[0][0], points[0][1])
        new_geom.AddPoint(points[1][0], points[1][1])
            # Create an output location for the drawn feature. This will likely be a temporary folder for the tool
        driver = ogr.GetDriverByName("ESRI Shapefile")
        datasource = driver.CreateDataSource("{}/{}.shp".format(self.output_location, name))
            # Create a layer for the transform data source
        transform_layer = datasource.CreateLayer(srs=UTM, name="Transform", geom_type=ogr.wkbLineString)
            # Grab each field name from the first layer definition
        for i in range(0, layer_defn.GetFieldCount()):
                # Grab the field definition from the layer
            field_name = layer_defn.GetFieldDefn(i)
            # Add the field to the transform layer
            # This retains attributes added by the user in the tool
            transform_layer.CreateField(field_name)
            # Add the newly created feature to the layer
        new_feature = ogr.Feature(transform_layer.GetLayerDefn())
        new_feature.SetGeometry(new_geom)
        # Grab the field values from the original feature and add it to the newly created transformed line
        for i in range(0, layer_defn.GetFieldCount()):
                field_name = layer_defn.GetFieldDefn(i)
                new_feature.SetField(field_name.name, base_feature.GetField(i))
        # Add the newly created feature to the layer
        transform_layer.CreateFeature(new_feature)
        self.line = "{}/{}.shp".format(self.output_location, name)

    def Create_Line_Buffer(self):
        """
        Creates a buffer around the transect for raster clipping


        :return: Shapefile of the buffer
        """
        # Get the line feature and load it up as a fiona collection
        with fiona.open(self.line) as src:
            # Call the first feature in the line feature
            f = src.next()
            # Call the coordinates from the geometry of the fiona feature
            cords = f["geometry"]['coordinates']
            # Convert the geometry into a shapely line string
            l = shapely.geometry.LineString(cords)
            # Draw a buffer around the
            buffer_geom = l.buffer(50, cap_style=3)
            Buffer_Geom = ogr.CreateGeometryFromWkt(buffer_geom.wkt)
        driver = ogr.GetDriverByName("ESRI Shapefile")
        # Create a data source for the buffer shapefile
        if os.path.exists("{}/Buffer.shp".format(self.output_location)):
            driver.DeleteDataSource("{}/Buffer.shp".format(self.output_location))
        buffer_data_source = driver.CreateDataSource(
            "{}/Buffer.shp".format(self.output_location))
        # Create a layer for the specified data source
        UTM = osr.SpatialReference()
        UTM.ImportFromEPSG(self.conversion_crs)
        layer = buffer_data_source.CreateLayer("Buffer.shp", srs=UTM, geom_type=ogr.wkbPolygon)
        feature_defn = layer.GetLayerDefn()
        # Create a new feature using the feature definition
        feature_1 = ogr.Feature(feature_defn)
        # Set the feature geometry from the newly created buffer
        feature_1.SetGeometry(Buffer_Geom)
        # Add the buffer feature to the created layer
        layer.CreateFeature(feature_1)
        self.buffer = "{}/Buffer.shp".format(self.output_location)

    def Clip_Raster(self):
        """
        Clips the raster using the buffer shapefile

        :return: Clipped raster
        """
        # Open a geopandas data frame
        gdf = gpd.read_file(self.buffer)
        # Convert the buffer shapefile into a geojson that can be used to clip a raster
        # Rasterio only uses geojsons for clipping
        shapes = [json.loads(gdf.to_json())['features'][0]['geometry']]
        with rasterio.open(self.raster) as r:
            out_image, out_transform = rasterio.mask.mask(r, shapes, crop=True)
            out_meta = r.meta
            out_meta.update({
                "driver": "Gtiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })
        with rasterio.open("{}/buffer.tif".format(self.output_location), 'w', **out_meta) as dst:
            dst.write(out_image)
        self.buffer_raster = "{}/buffer.tif".format(self.output_location)

    def Add_Raster_Centroid(self):
        """
        Create a point layer of raster pixel centroids. Centroids are derived from the clipped raster


        :return: Shapefile with vector points of each pixel centroid
        """
        # Open the clipped imagery
        r = gdal.Open(self.buffer_raster)
        # Open spatial reference
        UTM = osr.SpatialReference()
        UTM.ImportFromEPSG(26911)
        # Open the shapefile driver
        driver = ogr.GetDriverByName("ESRI Shapefile")
        # Get the transform information from the raster
        geom = r.GetGeoTransform()
        # Origin coord of X for the raster
        or_x = geom[0]
        # Origin coord of Y for the raster
        or_y = geom[3]
        # X width of the pixel (positive)
        pixel_width = geom[1]
        # Y width of the pixel (negative)
        pixel_height = geom[5]
        # Get the Y and X coordinates based on origin coords and pixel dimensions
        x = (or_x) + (pixel_width / 2)
        y = (or_y) + (pixel_height / 2)
        # Create point geometry
        geometry = ogr.Geometry(type=ogr.wkbPoint)
        # Create point geometry
        geometry.AddPoint(x, y)
        # Create a centroid shapefile
        Centroids_Shapefile = driver.CreateDataSource("{}/Centroids.shp".format(self.output_location))
        # Add a point layer to the shapefile
        layer = Centroids_Shapefile.CreateLayer("Centroids", srs=UTM, geom_type=ogr.wkbPoint)
        # Call the rows and columns of the raster
        cols = r.RasterXSize
        rows = r.RasterYSize
        # Get the layer definition of the centroids shapefile
        # Use this layer definition as the feature definition
        layer_defn = layer.GetLayerDefn()
        feature_1 = ogr.Feature(layer_defn)
        # Set the geometry of the feature with the point geometry
        feature_1.SetGeometry(geometry)
        # Add the feature to the centroids layer
        layer.CreateFeature(feature_1)
        # Destroy the feature so it isn't taking up memory
        feature_1.Destroy()
        for i in range(0, rows + 1):
            for j in range(0, cols + 1):
                x_cord = ((x) + (pixel_width * j))
                y_cord = ((y) + (pixel_height * i))
                # Create a new feature using the feature definition
                feature = ogr.Feature(layer_defn)
                    # Create point geometry to be added to the new feature
                point = ogr.Geometry(ogr.wkbPoint)
                    # Add x and y coords to the point geometry
                point.AddPoint(x_cord, y_cord)
                    # Set the x and y coords to the point geometry
                feature.SetGeometry(point)
                    # Add the feature to the centroids layer
                layer.CreateFeature(feature)
                feature.Destroy()
                point.Destroy()
        self.centroids = "{}/Centroids.shp".format(self.output_location)

    def Distance_To_Line(self):
        """
        Determine the distance of centroids to the transect. All centroids within 5 meters are kept, while all centroids more than 5 meters away
        are queried out and not used for sampling



        :return: Shapefile containing centroids that are within 5 meters of the transect
        """
        driver = ogr.GetDriverByName("ESRI Shapefile")
        # Call the point layer
        point_source = ogr.Open(self.centroids)
        layer = point_source.GetLayer()
        # Call the line layer
        line = ogr.Open(self.line)
        line_layer = line.GetLayer()
        # Call the line feature
        line_feature = line_layer.GetNextFeature()
        # Call the line features geometry
        line_geom = line_feature.GetGeometryRef()
        # Call the UTM srs to create a centroid data source
        UTM = osr.SpatialReference()
        UTM.ImportFromEPSG(26911)
        # Create a data source for the queried data
        query_source = driver.CreateDataSource("{}/Queried_Centroids.shp".format(self.output_location))
        # Create a layer for the query data source
        # This data source will be important for the created feature after the calculation is done
        # If the distance between the point and line is less than 5, it will be added to the data source
        query_layer = query_source.CreateLayer("Queried_Centroids", srs=UTM, geom_type=ogr.wkbPoint)
        # Get each point value from the line string
        if line_geom.GetX(0) > line_geom.GetX(1):
            x2 = line_geom.GetX(0)
            y2 = line_geom.GetY(0)
            y1 = line_geom.GetY(1)
            x1 = line_geom.GetX(1)
        else:
            x1 = line_geom.GetX(0)
            y1 = line_geom.GetY(0)
            x2 = line_geom.GetX(1)
            y2 = line_geom.GetY(1)
        # Calculate the distance between each point
        dx = x1 - x2
        dy = y1 - y2
        a = dy
        b = -dx
        c = y1 * dx - x1 * dy
        for i in range(0, layer.GetFeatureCount()):
            feature = layer.GetFeature(i)
            point_geom = feature.GetGeometryRef()
            x0 = point_geom.GetX()
            y0 = point_geom.GetY()
            d = abs(a * x0 + b * y0 + c) / math.sqrt(a * a + b * b)
            if d > 5:
                feature.Destroy()
            elif d < 5 and (x2 > x0 > x1):
                print(d)
                query_layer.CreateFeature(feature)
                feature.Destroy()
        self.centroids = "{}/Queried_Centroids.shp".format(self.output_location)

    def Sample_Pixel_Values_and_Calculate_Distance(self, dem = None):
        """
        Sample the pixel values at all centroids and determine the distance of each centroid from the transect origin.
        Distance to transect origin is for the reconstruction of the temperatures along the transect


        :param dem: Optional DEM for slop sampling of each centroid
        :return:
        """
        driver = ogr.GetDriverByName("ESRI Shapefile")
        # Call the point layer using geopandas
        source = gpd.read_file(self.centroids)
        # Create an empty list to contain each of the extracted values
        values = []
        # Create a list of x,y coordinates for distance calculations and value extraction
        coord_list = [(x, y) for x, y in zip(source['geometry'].x, source['geometry'].y)]
        # Open the clipped raster
        r = rasterio.open(self.buffer_raster)
        # For each x value in the coordinates list
        for x in coord_list:
            # Sample it's associated pixel value and append it to the value list
            for val in r.sample([x]):
                values.append(val[0])
        # Add the value list as a column in the geopandas frame
        values = [0 if i < 0 else i for i in values]
        source["Extracted_Values"] = values
        # Open the line shapefile source
        line_src = driver.Open(self.line)
        # Call the line layer
        line_layer = line_src.GetLayer()
        # Call the first feature in the line source
        feature = line_layer.GetFeature(0)
        # Get the geometry of the feature
        geom = feature.GetGeometryRef()
        x2 = geom.GetX(0)
        y2 = geom.GetY(0)
        distance = []
        for x in coord_list:
            x1 = x[0]
            y1 = x[1]
            d = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            distance.append(d)
        source["Distance"] = distance
        if dem is not None:
            dem_val = []
            with rasterio.open("C:/Users/alexm/OneDrive/Desktop/Preliminary Research/Thomas Fire Extent/ThomasDEM/ThomasDEM/slope.tif") as src:
                for x in coord_list:
                    # Sample it's associated pixel value and append it to the value list
                    for val in src.sample([x]):
                        dem_val.append(val[0])
            source["Slope"] = dem_val
        if os.path.exists("{}/Queried_Centroids.shp".format(self.output_location)):
            driver.DeleteDataSource("{}/Queried_Centroids.shp".format(self.output_location))
        source.to_file("{}/Queried_Centroids.shp".format(self.output_location))
        return source
    def Run_Transect_Sample(self):
        """
        Function that runs all methods at once

        :return:
        """
        self.Transform_Line_to_UTM()
        self.Create_Line_Buffer()
        self.Clip_Raster()
        self.Add_Raster_Centroid()
        self.Distance_To_Line()
        self.Sample_Pixel_Values_and_Calculate_Distance()

    def Output_Transect_and_Centroid(self, file_location):
        """
        Function that outputs the Transect and Centroid to a desired location

        :return: Transect and Centroids at desired file location
        """
        # Load the centroids and transect as a geopandas dataframe
        cent_gpd = gpd.read_file(self.centroids)
        transect_gpd = gpd.read_file(self.line)
        # Output the centroids and transect to desired file location
        cent_gpd.to_file("{}/Centroid.shp".format(file_location))
        transect_gpd.read_file("{}/Transect.shp".format(file_location))



# First open the transect as a


def Output_Transect_Figure(output_location, line, centroids, buffer, buffer_extent, pass_n):
    gpdf = gpd.read_file(centroids)
    sort = gpdf.sort_values(by=["Distance"])
    buf = gpd.read_file(buffer_extent)
    # Get the bounding box of the GeoDataFrame
    bbox = buf.bounds
    # Extract the bounding box coordinates
    max_value = sort["Extracted_"].max()
    min_value = sort["Extracted_"].min()
    mean_value = sort["Extracted_"].mean()
    xmin, ymin, xmax, ymax = bbox.loc[0]
    line_Df = gpd.read_file(line)
    fig, (ax1, ax2) = plt.subplots(1,2)
    ax2.scatter(x=sort["Distance"], y=sort["Extracted_"], c=sort["Extracted_"], cmap="plasma", zorder=2, s = 15, edgecolors = 'black')
    ax2.plot(sort["Distance"], sort["Extracted_"], c="red", zorder=1)
    ax2.set_title("Temperature Variations Across Transect")
    ax2.set_ylabel("Sampled Temperature (K)")
    ax2.set_xlabel("Distance From Transect Start (m)")
    fig.text(0.98, 0.75, "Max Temperature: {}".format(round(max_value,2)), ha='right', va='top', transform=plt.gca().transAxes, fontweight = 'bold', c = 'red')
    fig.text(0.98, 0.70, "Min Temperature: {}".format(round(min_value,2)), ha='right', va='top', transform=plt.gca().transAxes, fontweight = 'bold', c = 'red')
    fig.text(0.98, 0.60, "Mean Temperature: {}".format(round(mean_value,2)), ha='right', va='top', transform=plt.gca().transAxes, fontweight = 'bold', c = 'red')
    ax2.set_ylim(200,1100)
    r = rasterio.open(buffer)
    ax1.scatter(sort["geometry"].x, sort["geometry"].y, cmap="plasma", c=sort["Extracted_"], zorder = 2, edgecolors = 'black')
    rasterio.plot.show(source=r.read(1), ax=ax1, transform=r.transform, cmap="plasma", origin= 'upper')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin,ymax)
    ax1.set_title("Transect: Pass {}".format(pass_n))
    line_Df.plot(ax=ax1, color='black', linewidth=2, zorder = 1)
    plt.subplots_adjust(top=0.953,
            bottom=0.077,
            left=0.053,
            right=0.99,
            hspace=0.2,
            wspace=0.2)
    fig.set_size_inches(12.0, 9.0)
    plt.savefig("{}/{}".format(output_location, pass_n), dpi = 100)





def Transect_Animation(output_location, crs, line, raster_list):
    """
    Loops through a series of rasters, creates sample transect classes and outputs images of each sample transect class
    and sampled temperatures along the transect. Parameters are for initializing sampling transect class

    :param output_location: Location where transect and images will be saved
    :param crs: Conversion crs for transect
    :param line: Transect geojson drawn in folium
    :param raster_list: File path for rasters to be sampled
    :return: Series of images for each sample transect class
    """
    # Create a list of all raster paths
    r_list = List_Rasters(raster_list)
    # For each raster, initialize a sampling transect class for that list and sample the values from that raster
    for raster in r_list:
        Transect_Class = Sample_Transect(output_location=output_location,
                                 conversion_crs=crs,
                                 line=line,
                                 raster=raster)
        Transect_Class.Transform_Line_to_UTM()
        Transect_Class.Create_Line_Buffer()
        Transect_Class.Clip_Raster()
        Transect_Class.Add_Raster_Centroid()
        Transect_Class.Distance_To_Line()
        Transect_Class.Sample_Pixel_Values_and_Calculate_Distance()
        file_name = os.path.basename(raster)
        # Grab the pass number from the file name
        image_n = file_name[5:7]
        # Create a string that stores the pass number.
        pass_n = int(image_n)
        # Output the sampled transect to a matplotlib png
        Output_Transect_Figure(centroids=Transect_Class.centroids, buffer=Transect_Class.buffer_raster,
                               buffer_extent=Transect_Class.buffer,
                               line=Transect_Class.line, pass_n=pass_n,
                               output_location = output_location)



def Output_Transect_Figure_Slope(line, centroids, raster, buffer_extent):
    """
    Uses outputs from the Sample Transect to create figures visualizing temperature variations across transects


    :param line: Transect used for sampling
    :param centroids: Queried centroids from the sample transect class
    :param raster: Buffer clipped raster from the sampling transect class
    :param buffer_extent: Buffer extent used for clipping from the sampling transect class
    :return: Image showing variations in temperature across the transect
    """
    gpdf = gpd.read_file(centroids)
    sort = gpdf.sort_values(by=["Distance"])
    buf = gpd.read_file(buffer_extent)
    # Get the bounding box of the GeoDataFrame
    bbox = buf.bounds
    # Extract the bounding box coordinates
    xmin, ymin, xmax, ymax = bbox.loc[0]
    line_Df = gpd.read_file(line)
    fig, (ax1, ax2) = plt.subplots(1,2)
    ax2.scatter(x=sort["Distance"], y=sort["Slope"], c=sort["Slope"], cmap="turbo", zorder=2, s = 15, edgecolors = 'black')
    ax2.plot(sort["Distance"], sort["Slope"], c="red", zorder=1)
    ax2.set_title("Slope Variations Across Transect")
    ax2.set_ylabel("Slope (degrees)")
    ax2.set_xlabel("Distance From Transect Start (m)")
    ax2.set_ylim(0, 90)
    r = rasterio.open(raster)
    ax1.scatter(sort["geometry"].x, sort["geometry"].y, cmap="turbo", c=sort["Slope"], zorder = 2, edgecolors = 'black')
    rasterio.plot.show(source=r.read(1), ax=ax1, transform=r.transform, cmap="viridis", origin= 'upper')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin,ymax)
    ax1.set_title("Transect: Slope Aspect")
    line_Df.plot(ax=ax1, color='black', linewidth=2, zorder = 1)
    fig.figure.set_facecolor('#fa7a25')
    plt.subplots_adjust(top=0.953,
            bottom=0.077,
            left=0.053,
            right=0.99,
            hspace=0.1,
            wspace=0.1)
    fig.set_size_inches(18.5, 10.5)
    plt.savefig("C:/Users/alexm/OneDrive/Desktop/Figure Test/Transect 2/Slope", dpi = 100)



"""
Prototype Code below:
"""
def Output_Stacked_Plot(list_of_rasters, transect, output_image_location):
    # First, load up the list of transects into numerous file paths
    # Next, create the transect class. Use a temporary directory as the output location for the transects
    # Output each transect class to a table
    # Once the function is finished running, the temporary directory will be removed
    # Once every transect has been created, create a new function that can easily output numerous transects.
    # Same as the previous function, however the Grid Spec tool will be used to design a new output.
    # Create a temporary directory
    temp = tempfile.mkdtemp()
    try:
        # List the rasters for sampling
        rasters = List_Rasters(list_of_rasters)
        # List the transects that will be reprojected and output
        # transects = List_Geojsons(list_of_Transects) I will use this later on in development
        # Raster is called from the list of rasters
        # Either add a try statement based on the number of features in the Geojsons
        # For now, create a transect with multiple features and export each feature as an independent layer to a temp file that will be deleted
        # This list will contain temp file locations
        # Output location will always be the temporary folder
        transect_layer_list = []
        # Create a geojson driver
        driver = ogr.GetDriverByName("GeoJSON")
        # Open the transects layer
        transects = ogr.Open(transect)
        # Get the transects layer
        t_layer = transects.GetLayer()
        # Get the spatial reference for the layer so the data source has the correct srs
        srs = t_layer.GetSpatialRef()
        # Add a list to contain the bounds of all buffer layers and determine the max and min for all of them.
        buffer_list = []
        n = 0
        # Go through each feature and add a layer to a temp file
        colors = ["red", "blue", "green"]
        for i in range(0,t_layer.GetFeatureCount()):
            # Get the first feature in the transect layer
            f = t_layer.GetFeature(i)
            # Create a geojson source that will hold the features information
            # Each created layer will be converted to a transect class
            gjson_source = driver.CreateDataSource("{}/Layer_{}.geojson".format(temp, i))
            # Create a new layer using the geojson source
            new_layer = gjson_source.CreateLayer("Transect", srs=srs, geom_type =ogr.wkbLineString)
            # Create a new feature using the new layer
            new_layer.CreateFeature(f)
            # Append the file path to the geojson to the transect layer list
            # These paths will be used to create the transect classes
            transect_layer_list.append("{}/Layer_{}.geojson".format(temp, i))
            # Destroy the gjson_source so memory isn't constantly being used
            gjson_source.Destroy()
        for ra in rasters:
            n += 1
            # A list is created to house each of the geopandas dataframes that will contain the sampled centroids from the transect
            gpd_list = []
            # A list is created to house each of the geopandas dataframes that will contain the transects used for sampling
            gpd_transect_list = []
            # Create a list of point cords to determine the highest x and y coordinates for the output transect extent
            buffer_coords = [[]]
            for t in transect_layer_list:
                # Create a Transect class
                Transect_Class = Sample_Transect(line = t, raster=ra, conversion_crs=26911, output_location=temp)
                # Develop function that runs every single method and samples values around the Transect
                Transect_Class.Run_Transect_Sample()
                # Take the centroids created by the sampling and output them to a geopandas data frame.
                # This goepandas dataframe will be used for the subsequent transect graphs
                cent_gpd = gpd.read_file(Transect_Class.centroids)
                # Append the geopandas dataframe to the gpd_list
                gpd_list.append(cent_gpd)
                # Load the transect as a geopandas data frame
                trans_gpd = gpd.read_file(Transect_Class.line)
                # Append the geopandas dataframe to the transect list
                gpd_transect_list.append(trans_gpd)
                # Open the buffer and get the bounds of the buffer feature and append it to a list
                buf_df = gpd.read_file(Transect_Class.buffer)
                buffer_list.append(buf_df)
                # Use each of the geopandas dataframes and plot them using a line graph and scatter plot
                # First, need to determine a method to create grid specs the size of the number of transects
                # Ideally, the left side will have all the mapped transects and the right side will have all the extracted values
                # The colors of the mapped transects on the right side will coincide with the colors of the transects on the left side
            fig = plt.figure(constrained_layout = True)
                # Each Transect Plot needs 2 grid spec rows and 6 columns
                # Therefore, the length of the rows is determined by number of transects  * 2
                # In addition, the raster plot takes up 4 columns, while the
            gs = fig.add_gridspec(nrows = len(gpd_transect_list) * 2, ncols=6, left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
            cn = 0
            i = 0
            j = 2
            for d in gpd_list:
                ax = fig.add_subplot(gs[i:j, 3:])
                ax.set_ylim(0, 500)
                sort = d.sort_values(by=["Distance"])
                #ax.scatter(x=sort["Distance"], y=sort["Extracted_"], c=sort["Extracted_"], cmap="plasma", zorder=2,
                 #               s=15, edgecolors='black')
                ax.plot(sort["Distance"], sort["Extracted_"], color=colors[cn], zorder=1)
                ax.set_title("Temperature Variations Across Transect")
                ax.set_ylabel("Sampled Temperature (K)")
                ax.set_xlabel("Distance From Transect Start (m)")
                i += 2
                j += 2
                cn += 1
            # Now create the min max bounds of the map figure
            # Combine the geodataframes in the list and get the bounds for the min max values
            buffer_df = pd.concat(buffer_list, ignore_index=True)
            max_y = max(buffer_df["geometry"].bounds.maxy.values.tolist())
            min_y = min(buffer_df["geometry"].bounds.miny.values.tolist())
            max_x = max(buffer_df["geometry"].bounds.maxx.values.tolist())
            min_x = min(buffer_df["geometry"].bounds.minx.values.tolist())
            print(max_y)
            print(min_y)
            print(max_x)
            print(min_x)
            # [minx, miny, maxx, maxy]
            r = rasterio.open(ra)
            ax1 = fig.add_subplot(gs[:j, :3])
            rasterio.plot.show(source=r.read(1), ax=ax1, transform=r.transform, cmap="plasma", origin='upper')
            plt.subplots_adjust(top=0.953,
                                bottom=0.077,
                                left=0.053,
                                right=0.99,
                                hspace=0.1,
                                wspace=0.1)
            fig.set_size_inches(18.5, 10.5)
            cn_2 = 0
            for d in gpd_transect_list:
                d.plot(ax=ax1, color=colors[cn_2], linewidth=2, zorder=1)
                cn_2 += 1
            ax1.set_xlim(min_x, max_x)
            ax1.set_ylim(min_y, max_y)
            plt.savefig("{}/No_Scatter_Figure_{}.png".format(output_image_location, n))
        shutil.rmtree(temp)
    except Exception as e:
        shutil.rmtree(temp)
        print(e)

# Outputs drawn transect and centroids to file location








