# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ by Gonzalo Simarro and Daniel Calvete in 2025
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
import ee
import json
import os
import sys
import warnings
warnings.simplefilter("error", RuntimeWarning)
#
sys.path.insert(0, 'uortos')
import uortos as uortos
#
#
pathFldMain = 'example'
#
assert os.path.exists(pathFldMain)
#
#
pathJson = os.path.join(pathFldMain, 'data', 'parameters.json')
with open(pathJson, 'r') as f:
    par = json.load(f)
#
pathJson = os.path.join(pathFldMain, 'data', 'satellites.json')
with open(pathJson, 'r') as f:
    par.update(json.load(f))
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ Google Earth Engine initialization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
print('GEE initialization')
ee.Authenticate()
ee.Initialize(project=par['eeProject'])
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ download of zip, tif and png
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
print('Download of zip, tif and png files')
uortos.DownloadImages(pathFldMain, par)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ obtain FM pairs and purge
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
print('Obtain FM pairs')
uortos.FMPairs(pathFldMain, par, par['limFM'] > 0)
#
print('Purge FM pairs')
uortos.PurgePairs(pathFldMain, par, 'fm')
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ obtain cor pairs and purge
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
print('Obtain corr pairs')
uortos.CorPairs(pathFldMain, par)
#
print('Purge corr pairs')
uortos.PurgePairs(pathFldMain, par, 'cor')
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ obtain affine transformations and final products
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
print('Obtain affine transformations')
uortos.Affines(pathFldMain, par)
#
print('Obtain final outputs')
uortos.Output(pathFldMain, par)
#
