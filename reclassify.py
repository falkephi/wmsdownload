#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Reclassify RGB GeoTIFF layer into one songle band.

"""
import os
import struct

from osgeo import gdal

import argparse as ap
import numpy as np


class RGBlayer(object):

    def __init__(self):
        self.outfile = None
        self.inputdata = None
        self.outputdata = None

    def load(self, filename):
        filename = os.path.expanduser(filename)
        self.infile = gdal.Open(filename)
        if self.infile is not None:
            self.bands = self.infile.RasterCount
            self.geotransform = self.infile.GetGeoTransform()
            self.proj = self.infile.GetProjection()
            self.x_size = self.infile.RasterXSize
            self.y_size = self.infile.RasterYSize
            self.inputdata = self.infile.ReadAsArray()
            # Close input file
            self.infile = None

    def check_outfile(self, outfile, overwrite=False):
        outfile = os.path.expanduser(outfile)
        if os.path.exists(outfile) and overwrite is False:
            raise IOError(
                'File exists! Please choose a different file name.')

    def save(self, outfile, overwrite=False):
        if self.outputdata is not None:
            self.outfile = os.path.expanduser(outfile)
            self.check_outfile(self.outfile, overwrite=overwrite)
            if overwrite:
                os.remove(outfile)

            driver = gdal.GetDriverByName('GTiff')
            target_file = driver.Create(self.outfile,
                                        int(self.x_npix),
                                        int(self.y_npix),
                                        1, gdal.GDT_UInt32)
            target_file.SetGeoTransform(self.geotransform)
            target_file.SetProjection(self.proj)

            write_band = struct.pack('I' * self.x_size * self.y_size,
                                     *self.outputdata)
            target_band = target_file.GetRasterBand(1)
            target_band.WriteRaster(0, 0, self.x_size, self.y_size,
                                    write_band, self.x_size, self.y_size,
                                    gdal.GDT_UInt32)

    def reclassify(self):

        # Create unique values from all bands
        sumband = self.inputdata[0, :, :].flatten()
        for band in range(1, self.bands):
            sumband = sumband + self.inputdata[band, :, :].flatten() \
                * 256 ** band
        # Retrieve unique values and map them to integer values
        unique_values = np.unique(sumband)
        mapping_values = range(len(unique_values))
        for i, val in enumerate(unique_values):
            # replace val in sumband with mapping_values[i]
            sumband[sumband == val] = mapping_values[i]
        self.outputdata = tuple(sumband)


def commandline_parser():
    parser = ap.ArgumentParser(
        description="""Reclassify a multiband GeoTIFF into a single band
GeoTIFF with unique classes."

(c) Copyright 2014, Philipp Meier <philipp@diemeiers.ch>
        """, formatter_class=ap.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-file',
                        help='Input file.')
    parser.add_argument('-o', '--output-file',
                        help='File to which data will be saved.')
    parser.add_argument('-x', '--overwrite', action='store_true',
                        help='Overwrite existing file.')

    return parser

if __name__ == '__main__':
    parser = commandline_parser()
    args = parser.parse_args()

    layer = RGBlayer()
    # Throw an error if file exists before heavy lifting is done
    layer.check_outfile(args.output_file, overwrite=args.overwrite)
    layer.load(args.input_file)
    layer.reclassify()
    layer.save(args.output_file, overwrite=args.overwrite)
