wmsdownload
===========

Python tool for downloading WMS data for offline use.

DISCLAIMER
----------

Be aware that some WMS services have licensing restrictions which prohibit
downloading the data for offline use.

How to use it
-------------

```
python wmsdownload.py [-h] -u URL [-l LAYER [-f DATA_FORMAT]
                      -r RESOLUTION -o OUTPUT [-p PROJECTION] [-x] |
                      [-a] [-i INFO] [-s SEARCH]]

  -h, --help            show this help message and exit
  -u URL, --url URL     URL of the WMS server.
  -l LAYER, --layer LAYER
                        Name of the layer that should be downloaded.
  -f DATA_FORMAT, --data-format DATA_FORMAT
                        Format of the data to be downloaded (raster or
                        vector).
  -r RESOLUTION, --resolution RESOLUTION
                        Set resolution (in Meters) of the data that will be
                        downloaded.
  -o OUTPUT, --output OUTPUT
                        File to which data will be saved.
  -x, --overwrite       Overwrite existing file.
  -a, --list-all        Print a list of available layers.
  -i INFO, --info INFO  Get more info about a specific layer.
  -s SEARCH, --search SEARCH
                        Find a keyword in the available layers. Searches in
                        name, title and abstract.
  -p PROJECTION, --projection PROJECTION
                        Set the projection of the output file as EPSG code
                        (e.g. WGS84 --> EPSG:4326; CH1903/LV09 -->
                        EPSG:21781).
  -bb BOUNDINGBOX, --boundingbox BOUNDINGBOX
                        Specify the area of interest for downloading data.
                        Format: minx,miny,maxx,maxy (no spaces between numbers
                        and comma!).

```

LICENSE
-------

This software can be distributed freely under the GPL v2 license. Please see the
LICENSE file for more information.

(c) Copyright 2014, Philipp Meier.

