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

    def reclassify(self, n, classify_ambiguous=False):

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
        pixclasses = frequencies[:n]
        classified = np.zeros(sumband.shape, np.uint8)

        # get colors
        if classify_ambiguous:
            red = self.inputdata[0, :, :].flatten()
            green = self.inputdata[1, :, :].flatten()
            blue = self.inputdata[2, :, :].flatten()
            img_class_colors = dict()

        for i, pixclass in enumerate(pixclasses):
            mask = sumband == pixclass[0]
            classified[mask] = i
            if classify_ambiguous:
                img_class_colors[i] = dict(r=red[mask][0],
                                           g=green[mask][0],
                                           b=blue[mask][0])

        if classify_ambiguous:
            print 'Classify ambiguous pixels (go get a cup of coffee) ...'
            # Now it's time to handle ambiguous cases
            # first convert RGB to HSV:
            self.progress(0.0)
            for i in img_class_colors.keys():
                color = img_class_colors[i]
                r = float(color['r']) / 255
                g = float(color['g']) / 255
                b = float(color['b']) / 255
                H, S, V = RGBtoHSV(r, g, b)
                img_class_colors[i].update(H=H, S=S, V=V)

            it = 0
            for val, count in frequencies[n:]:
                mask = sumband == val
                r = float(red[mask][0]) / 255
                g = float(green[mask][0]) / 255
                b = float(blue[mask][0]) / 255
                H, S, V = RGBtoHSV(r, g, b)
                old_match = 280
                match_idx = -1
                for i in img_class_colors.keys():
                    match = get_color_proximity((H, S, V),
                                                (img_class_colors[i]['H'],
                                                 img_class_colors[i]['S'],
                                                 img_class_colors[i]['V']))
                    if match < old_match:
                        old_match = match
                        match_idx = i
                classified[mask] = match_idx
                self.progress(float(it) / len(frequencies))
                it += 1
            self.progress(1.0)

        self.outputdata = tuple(classified)


def RGBtoHSV(r, g, b):
    cmax = max([r, g, b])
    cmin = min([r, g, b])
    delta = cmax - cmin
    if delta == 0:
        S = 0
        H = 0
    else:
        S = delta / cmax
        if r >= g and r >= b:
            H = 60 * ((g - b) / delta % 6)
        elif g >= r and g >= b:
            H = 60 * ((b - r) / delta + 2)
        elif b >= r and b >= g:
            H = 60 * ((r - g) / delta + 4)

    V = cmax
    return H, S, V


def get_color_proximity(color1, color2):
    diff1 = min(abs(color1[0] - color2[0]), abs(360 - color2[0] + color1[0]))
    diff2 = abs(color1[1] - color2[1]) * 10
    diff3 = abs(color1[2] - color2[2]) * 10

    return diff1 + diff2 + diff3


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
    parser.add_argument('-a', '--reclass-ambiguous', action='store_true',
                        help='Classify ambiguous pixels (might take a while).')

    return parser


if __name__ == '__main__':
    parser = commandline_parser()
    args = parser.parse_args()

    layer = RGBlayer()
    # Throw an error if file exists before heavy lifting is done
    layer.check_outfile(args.output_file, overwrite=args.overwrite)
    layer.load(args.input_file)
    layer.reclassify(int(args.classes),
                     classify_ambiguous=args.reclass_ambiguous)
    layer.save(args.output_file, overwrite=args.overwrite)
