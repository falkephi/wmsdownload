#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Download WMS data for offline use.

"""

import os
import tempfile
import glob
import textwrap

from osgeo import gdal
from osgeo import osr
from osgeo import ogr

from owslib.wms import WebMapService

import argparse as ap


class WmsDataset(object):

    def __init__(self, servername):
        self.wms = WebMapService(servername)
        self.layers = self.wms.contents
        #self.epsg = 'EPSG:4326'
        self.epsg = 'EPSG:21781'
        self.datafiles = []
        self.boundingbox = None
        self.resolution = None
        self.nx = 1024
        self.ny = 1024
        self.has_tmp = False

        self.proj = get_projection(self.epsg)

        self.progress = gdal.TermProgress_nocb

    def list_layers(self, layer=None, out_format='long'):
        """List all available layers. If a layer is specified it prints all
        information available."""
        if layer is None:
            for key in self.layers.keys():
                print u'{:<35} {:<43}'.format(key, self.layers[key].title)
        else:
            if out_format == 'long':
                print layer
                print
                print self.layers[layer].title
                print
                print textwrap.fill(self.layers[layer].abstract, 80)
                print
            elif out_format == 'short':
                print u'{:<35} {:<43}'.format(layer, self.layers[layer].title)

    def search_layers(self, string):
        """Search for layers containing 'string' in the layer name, title or
        abstract."""
        for key in self.layers.keys():
            if string in key or string in self.layers[key].title \
                    or string in self.layers[key].abstract:
                self.list_layers(layer=key, out_format='short')

    def create_tmp_folder(self):
        """Create a temp folder only when we are downloading something."""
        self.temp_folder = tempfile.mkdtemp(prefix='WMSdata_', dir='/tmp/')
        self.has_tmp = True

    def set_projection(self, epsg_code):
        self.epsg = epsg_code
        self.proj = get_projection(self.epsg)

    def download_data(self, layer, data_format='raster', **kwargs):
        """Download tiles of the specified layer for the full extent."""
        self.create_tmp_folder()
        if data_format == 'raster':
            if 'output' in kwargs.keys() and kwargs['output'] is not None:
                outfile = kwargs['output']
                del kwargs['output']
            #if 'reclass' in kwargs.keys():
            #    reclass = kwargs['reclass']
            #    del kwargs['reclass']
            #else:
            #    reclass = False
            if 'overwrite' in kwargs.keys():
                overwrite = kwargs['overwrite']
                del kwargs['overwrite']
            else:
                outfile = None
            self.download_raster(layer, **kwargs)
            if outfile is not None:
                self.merge_tiles(outfile, overwrite=overwrite)
                #self.merge_tiles(outfile, overwrite=overwrite,
                #                 reclass=reclass)
        elif data_format == 'vector':
            return self.download_vector(layer, **kwargs)

    def download_raster(self, layer, resolution=1, nx=None, ny=None,
                        boundingbox=None):
        #                reclass=False):
        """Download a raster dataset by downloading different tiles and merging
        them. Reclassify unique RGB value if necessary."""
        print "Downloading data ..."
        if nx is None:
            nx = self.nx
        else:
            self.nx = int(nx)
        if ny is None:
            ny = self.ny
        else:
            self.ny = int(ny)
        # calculate the number and bounds of the tiles to download
        self.resolution = resolution
        if boundingbox is None:
            self.boundingbox = self.layers[layer].boundingBox
            self.boundingbox = reproject_boundingbox(self.boundingbox,
                                                     self.proj)
        else:
            self.boundingbox = boundingbox
        minx = self.boundingbox[0]
        miny = self.boundingbox[1]
        maxx = self.boundingbox[2]
        maxy = self.boundingbox[3]
        self.x_extent = maxx - minx
        self.y_extent = maxy - miny
        self.x_npix = self.x_extent / resolution
        self.y_npix = self.y_extent / resolution
        x_ntiles = int(self.x_npix // nx)  # number of tiles in x direction
        if self.x_npix % nx > 0:
            x_ntiles += 1
        y_ntiles = int(self.y_npix // ny)  # number of tiles in y direction
        if self.y_npix % ny > 0:
            y_ntiles += 1

        x_tile_bounds = []
        for xt in range(x_ntiles):
            x_tile_bounds.append(minx + xt * nx * resolution)
        x_tile_bounds.append(maxx)

        y_tile_bounds = []
        for yt in range(y_ntiles):
            y_tile_bounds.append(miny + yt * ny * resolution)
        y_tile_bounds.append(maxy)

        self.progress(0.)
        tottiles = x_ntiles * y_ntiles
        n = 0
        for i in range(x_ntiles):
            xbounds = (x_tile_bounds[i], x_tile_bounds[i + 1])
            x_size = (xbounds[1] - xbounds[0]) / self.resolution
            for j in range(y_ntiles):
                ybounds = (y_tile_bounds[j], y_tile_bounds[j + 1])
                y_size = (ybounds[1] - ybounds[0]) / self.resolution
                filename = self.download_tile(layer,
                                              (xbounds[0], ybounds[0],
                                               xbounds[1], ybounds[1]),
                                              (x_size, y_size))
                self.datafiles.append(filename)
                n += 1
                self.progress(n / float(tottiles))
        self.progress(1.0)

    def download_vector(self, layer, **kwargs):
        output = kwargs['output']
        overwrite = kwargs['overwrite']

        output = os.path.expanduser(output)
        if os.path.exists(output):
            if overwrite:
                os.remove(output)
            else:
                raise IOError(
                    'File exists! Please choose a different file name.')

        print "Downloading data ..."

        vecmap = self.layers[layer]
        self.boundingbox = vecmap.boundingBox
        self.boundingbox = reproject_boundingbox(self.boundingbox, self.proj)

        kmldata = self.wms.getmap(
            layers=[vecmap.name],
            bbox=self.boundingbox,
            srs=self.epsg,
            format='application/vnd.google-earth.kml+xml',
            size=(1000, 1000),
            transparent=True)

        # KML files support lat / lon coordinates only:
        # https://developers.google.com/kml/documentation/kmlreference
        # We convert the file to an ESRI shapefile to have the correct
        # projection.
        fhandle, fname = tempfile.mkstemp(suffix='.kml', dir=self.temp_folder)
        os.write(fhandle, kmldata.read())
        os.close(fhandle)

        print "Reprojecting data ..."
        # Input file
        kml_driver = ogr.GetDriverByName('KML')
        kml_datasource = kml_driver.Open(fname)
        kml_layer = kml_datasource.GetLayer()

        # Output file
        shp_driver = ogr.GetDriverByName('ESRI Shapefile')
        shp_datasource = shp_driver.CreateDataSource(output)
        shp_layer = shp_datasource.CreateLayer(
            kml_layer.GetName(), geom_type=kml_layer.GetGeomType())

        # Copy fields
        kml_layerdefn = kml_layer.GetLayerDefn()
        for i in range(kml_layerdefn.GetFieldCount()):
            field_defn = kml_layerdefn.GetFieldDefn(i)
            shp_layer.CreateField(field_defn)

        shp_layerdefn = shp_layer.GetLayerDefn()

        # Define geometric transformation
        spatial_ref = kml_layer.GetSpatialRef()
        transformation = osr.CoordinateTransformation(spatial_ref, self.proj)

        # Copy and reproject all features
        kml_feature = kml_layer.GetNextFeature()
        while kml_feature:
            geom = kml_feature.GetGeometryRef()
            geom.Transform(transformation)
            shp_feature = ogr.Feature(shp_layerdefn)
            shp_feature.SetGeometry(geom)
            # Copy attributes
            for i in range(shp_layerdefn.GetFieldCount()):
                shp_feature.SetField(
                    shp_layerdefn.GetFieldDefn(i).GetNameRef(),
                    kml_feature.GetField(i))
            shp_layer.CreateFeature(shp_feature)
            kml_feature.Destroy()
            shp_feature.Destroy()
            kml_feature = kml_layer.GetNextFeature()

        kml_datasource.Destroy()
        shp_datasource.Destroy()

    def download_tile(self, layer, bbox, size):
        img = self.wms.getmap(layers=[layer],
                              styles=['default'],
                              srs=self.epsg,
                              bbox=bbox,
                              size=size,
                              format='image/tiff',
                              transparent=True)

        fhandle, fname = tempfile.mkstemp(suffix='.tif', dir=self.temp_folder)
        os.write(fhandle, img.read())
        os.close(fhandle)
        return fname

    #def merge_tiles(self, output, overwrite=True, reclass=False):
    def merge_tiles(self, output, overwrite=True):
        """Merge tiles downloaded to the tmp directory into one single file at
        the specified location (``output``). Overwrite file by default (can be
        changed by setting ``overwrite`` to ``False``."""

        print "Merging data ..."
        output = os.path.expanduser(output)
        if os.path.exists(output):
            if overwrite:
                os.remove(output)
            else:
                raise IOError(
                    'File exists! Please choose a different file name.')

        # open the first valid file
        source_file = gdal.Open(self.datafiles[0])
        i = 1
        while source_file is None:
            source_file = gdal.Open(self.datafiles[i])
            i += 1
        bands = source_file.RasterCount
        #if reclass:
        #    band_type = gdal.GDT_UInt32
        #    out_bands = 1
        #else:
        #    band_type = source_file.GetRasterBand(1).DataType
        #    out_bands = bands
        band_type = source_file.GetRasterBand(1).DataType
        out_bands = bands

        # Create a new file:
        driver = gdal.GetDriverByName('GTiff')
        target_file = driver.Create(output,
                                    int(self.x_npix),
                                    int(self.y_npix),
                                    out_bands, band_type)
        geotransform = [self.boundingbox[0], self.resolution, 0,
                        self.boundingbox[3], 0, -self.resolution]
        target_file.SetGeoTransform(geotransform)
        target_file.SetProjection(self.proj.ExportToWkt())

        self.progress(0)
        for i, datfile in enumerate(self.datafiles):
            ff = gdal.Open(datfile)
            if ff is not None:
                # Transform coordinates to pixel coordinates
                geotransform = ff.GetGeoTransform()
                ulx = geotransform[0]
                uly = geotransform[3]
                ulbbx = self.boundingbox[0]
                ulbby = self.boundingbox[3]
                x_off = int((ulbbx - ulx) / -self.resolution)
                y_off = int((ulbby - uly) / self.resolution)
                x_size = ff.RasterXSize
                y_size = ff.RasterYSize
                #if reclass is True:
                #    # Do the reclassification
                #    allbands = ff.ReadAsArray()
                #    sumband = allbands[0, :, :].flatten()
                #    for band in range(1, bands):
                #        sumband = sumband + allbands[band, :, :].flatten() \
                #            * 256 ** band
                #    # Retrieve unique values and map them to integer values
                #    unique_values = np.unique(sumband)
                #    mapping_values = range(len(unique_values))
                #    for i, val in enumerate(unique_values):
                #        # replace val in sumband with mapping_values[i]
                #        sumband[sumband == val] = mapping_values[i]
                #    sumband = tuple(sumband)
                #    write_band = struct.pack('I' * x_size * y_size, *sumband)
                #    target_band = target_file.GetRasterBand(1)
                #    target_band.WriteRaster(x_off, y_off, x_size, y_size,
                #                            write_band, x_size, y_size,
                #                            gdal.GDT_UInt32)
                #else:
                #    # Proceed with copying single bands
                for band in range(bands):
                    source_band = ff.GetRasterBand(band + 1)
                    target_band = target_file.GetRasterBand(band + 1)
                    source_data = \
                        source_band.ReadRaster(0, 0, x_size, y_size,
                                               x_size, y_size,
                                               target_band.DataType)
                    target_band.WriteRaster(x_off, y_off, x_size, y_size,
                                            source_data, x_size, y_size,
                                            target_band.DataType)
            # Close input file
            ff = None
            self.progress(float(i) / len(self.datafiles))

        # Close target file
        target_file = None

        self.progress(1.0)
        print "Layer saved to %s." % output

    def destroy(self):
        """Empty and delete the temporary directory."""
        if self.has_tmp:
            files = glob.glob(self.temp_folder + '/*')
            for f in files:
                os.remove(f)

            os.removedirs(self.temp_folder)


def get_projection(epsgstring):
    _, epsg = epsgstring.split(':', 1)
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(int(epsg))
    return proj


def reproject_boundingbox(bbox, target_proj):
    if len(bbox) == 5:
        source_proj = get_projection(bbox[4])
        transform = osr.CoordinateTransformation(source_proj, target_proj)
        point1 = ogr.CreateGeometryFromWkt("POINT (%s %s)" %
                                           (bbox[0], bbox[1]))
        point2 = ogr.CreateGeometryFromWkt("POINT (%s %s)" %
                                           (bbox[2], bbox[3]))
        point1.Transform(transform)
        point2.Transform(transform)

        return (point1.GetX(), point1.GetY(), point2.GetX(), point2.GetY())


def commandline_parser():
    parser = ap.ArgumentParser(
        description="""Download Web Map Service (WMS) data from for offline
use. Be aware that some WMS services have licensing restrictions which prohibit
downloading the data for offline use.

(c) Copyright 2014, Philipp Meier <philipp@diemeiers.ch>
        """, formatter_class=ap.RawDescriptionHelpFormatter)
    parser.add_argument('-u', '--url',
                        help='URL of the WMS server.')
    parser.add_argument('-l', '--layer',
                        help='Name of the layer that should be downloaded.')
    parser.add_argument('-f', '--data-format',
                        help='''Format of the data to be downloaded (raster or
                        vector).''')
    parser.add_argument('-r', '--resolution',
                        help='''Set resolution (in Meters) of the data that
                        will be downloaded.''')
    parser.add_argument('-o', '--output',
                        help='File to which data will be saved.')
    parser.add_argument('-x', '--overwrite', action='store_true',
                        help='Overwrite existing file.')

    parser.add_argument('-a', '--list-all', action='store_true',
                        help='Print a list of available layers.')
    parser.add_argument('-i', '--info',
                        help='Get more info about a specific layer.')
    parser.add_argument('-s', '--search',
                        help='''Find a keyword in the available layers.
                        Searches in name, title and abstract.''')
    parser.add_argument('-p', '--projection',
                        help='''Set the projection of the output file as
                        EPSG code (e.g. WGS84 --> EPSG:4326; CH1903/LV09 -->
                        EPSG:21781).''')
    parser.add_argument('-bb', '--boundingbox',
                        help='''Specify the area of interest for downloading
                        data. Format: minx,miny,maxx,maxy (no spaces between
                        numbers and comma!).''')
    #parser.add_argument('-c', '--reclassify', action='store_true',
    #                    help='''Reclassify unique RGB values to integer
    #                    classes [1, 2, ..., z].''')

    return parser


if __name__ == '__main__':
    parser = commandline_parser()
    args = parser.parse_args()

    dd = WmsDataset(args.url)

    if args.list_all:
        dd.list_layers(out_format='short')
    elif args.info is not None:
        dd.list_layers(layer=args.info, out_format='long')
    elif args.search is not None:
        dd.search_layers(args.search)
    else:
        layer = args.layer
        if args.data_format is not None:
            data_format = args.data_format
        else:
            data_format = 'raster'

        #if args.reclassify is not None:
        #    reclassify = args.reclassify
        #else:
        #    reclassify = False

        outfile = args.output

        if args.resolution is None:
            res = 25
        else:
            res = float(args.resolution)

        if args.projection is not None:
            dd.set_projection(args.projection)

        if args.boundingbox is not None:
            boundingbox = \
                tuple([float(i) for i in args.boundingbox.split(',')])
        else:
            boundingbox = None

        try:
            dd.download_data(layer,
                             data_format=data_format,
                             resolution=res,
                             output=args.output,
                             overwrite=args.overwrite,
                             boundingbox=boundingbox)
                             #reclass=reclassify)
        except KeyboardInterrupt:
            print "\nCaught keyboard interrupt."
            print "Cleaning up ..."
        except:
            print "An error occurred!"
            print "Cleaning up ..."
            dd.destroy()
            raise
    dd.destroy()
