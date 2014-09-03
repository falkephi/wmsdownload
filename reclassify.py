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
        self.progress = gdal.TermProgress_nocb

    def load(self, filename):
        print "Loading file ..."
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
        exists = os.path.exists(outfile)
        if exists and not overwrite:
            raise IOError(
                'File exists! Please choose a different file name.')
        return exists

    def save(self, outfile, overwrite=False):
        print "Saving data ..."
        if self.outputdata is not None:
            self.outfile = os.path.expanduser(outfile)
            exists = self.check_outfile(self.outfile, overwrite=overwrite)
            if exists and overwrite:
                os.remove(outfile)

            driver = gdal.GetDriverByName('GTiff')
            target_file = driver.Create(self.outfile,
                                        self.x_size,
                                        self.y_size,
                                        1, gdal.GDT_UInt32)
            target_file.SetGeoTransform(self.geotransform)
            target_file.SetProjection(self.proj)

            write_band = struct.pack('I' * self.x_size * self.y_size,
                                     *self.outputdata)
            target_band = target_file.GetRasterBand(1)
            target_band.WriteRaster(0, 0, self.x_size, self.y_size,
                                    write_band, self.x_size, self.y_size,
                                    gdal.GDT_UInt32)

    def reclassify(self, n):

        print "Classify ..."
        # Create unique values from all bands
        sumband = self.inputdata[0, :, :].flatten()
        for band in range(1, self.bands):
            sumband = sumband + self.inputdata[band, :, :].flatten() \
                * 256 ** band
        # Retrieve unique values and map them to integer values
        unique_values, idx = np.unique(sumband, return_inverse=True)
        count = np.bincount(idx)
        frequencies = zip(unique_values, count)
        frequencies.sort(key=lambda tup: tup[1], reverse=True)
        # The first n most frequent classes are the classes we need:
        pixclasses = frequencies[0:n]
        classified = np.zeros(sumband.shape, np.uint8)
        for i, pixclass in enumerate(pixclasses):
            classified[sumband == pixclass[0]] = i
        # Now it's time to handle ambiguous cases
        # will be added later

        self.outputdata = tuple(classified)


def commandline_parser():
    parser = ap.ArgumentParser(
        description="""Reclassify a multiband GeoTIFF into a single band
GeoTIFF with unique classes.

(c) Copyright 2014, Philipp Meier <philipp@diemeiers.ch>
        """, formatter_class=ap.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-file',
                        help='Input file.')
    parser.add_argument('-o', '--output-file',
                        help='File to which data will be saved.')
    parser.add_argument('-x', '--overwrite', action='store_true',
                        help='Overwrite existing file.')
    parser.add_argument('-c', '--classes',
                        help='Number of classes.')

    return parser


if __name__ == '__main__':
    parser = commandline_parser()
    args = parser.parse_args()

    layer = RGBlayer()
    # Throw an error if file exists before heavy lifting is done
    layer.check_outfile(args.output_file, overwrite=args.overwrite)
    layer.load(args.input_file)
    layer.reclassify(int(args.classes))
    layer.save(args.output_file, overwrite=args.overwrite)
