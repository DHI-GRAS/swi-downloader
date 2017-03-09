#Definition of inputs and outputs
#==================================
##GWA scripts=group
##Download Soil Water Index=name
##ParameterSelection|product_param|Soil Water Product|Daily;10 Day Composite
##ParameterExtent|extent|Extent to subset after downloading (geographic coordinates)|''|True
##ParameterString|startdate|Start date (YYYYMMDD) - if left empty the default start date will be used|
##ParameterString|enddate|End date (YYYYMMDD) - if left empty today will be used as end date|
##ParameterString|username|User name|
##ParameterString|password|Password|
##ParameterBoolean|delete_downloaded|Delete downloaded files|True
##ParameterBoolean|load_to_canvas|Open output file after running algorithm|True
##ParameterFile|outdir|Output directory|True|False|
import os
import sys
from processing.tools import dataobjects
here = os.path.dirname(scriptDescriptionFile)
if here not in sys.path:
    sys.path.append(here)

import download_swi
from logs import set_qgis_logger

set_qgis_logger(progress)

product = ['SWI', 'SWI10'][product_param]

outfiles = download_swi.download(
        outdir=outdir, progress=progress,
        product=product,
        username=username, password=password,
        startdate=startdate, enddate=enddate, extent=extent,
        delete_downloaded=delete_downloaded)

if load_to_canvas:
    for outfile in outfiles:
        dataobjects.load(outfile, os.path.basename(outfile), style=None)

# Release file pointers held by QGIS Processing
dataobjects.resetLoadedLayers()
