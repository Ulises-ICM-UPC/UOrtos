# UOrtos

`UOrtos` is an open source software to obtain a set of outlier-free co-registered and georeferenced satellite images with high alignment quality.

### Description
The algorithm for co-registration and georeferentation of satellite images assumes that images can be realigned through an isometric transformation. The image alignment is achieved by minimising the distances between common image features. Clustering techniques are used to obtain a high quality set of aligned images without outliers. As a result of the process, the metadata information related to the geotransformation of the satellite images (GeoTIFF) is corrected and co-registered PNG images are derived. The development of this software is aimed for obtaining images to accurately determine shoreline displacements. In general, it can be used for the precise detection of features at different times. Details about the algorithm and methodology are described in
> *Simarro, G.; Calvete, D.; Ribas, F.; Castillo Y.; Puig-Polo, C. UOrtos: A Python Tool for Subpixel co-registration of Landsat and Sentinel 2 Imagery. Remote Sens. 2025, xx, xxxx. https://doi.org/10.3390/rs12345678*

The co-registration process consists of the following steps:

1. [Image downloading](#image-download)
2. [Feature pairing](#feature-pairing)
3. [Image clustering](#clustering)
4. [Outputs](#output)

### Requirements and project structure
To run the software it is necessary to have Python (3.12) and install the following dependencies:
- cv2 (4.1.0)
- numpy (2.1.0)
- scipy (1.14.1)
- osgeo (3.9.1)
- ee (1.418)

In parenthesis we indicate the version with which the software has been tested. It is possible that it works with older versions.

The structure of the project is the following:
* `example.py`
* `example_notebook.py`
* **`uortos`**
  * `uortos.py`
  * `ulises_uortos.py`
* **`example`**
  * **`data`**
    * `parameters.json`
    * `satellites.json`
    * `<res01>_<domain01>_reference.tif`
  * **`output`**
    * **`<parameter_set_01>`**
      * **`<res01>_<domain01>`**
        * `<time-stamp>.log`
        * `<ref_image>.ref`
        * **`geopng`**
          * `<image01>.png`
          * . . .
        * **`tif`**
          * **`<image01>`**
            * `<image01>_<band01>.tif`
            * . . .
          * . . .
      * . . .
    * . . .
  * **`scratch`**
    * . . .

The local modules of `UOrtos` are located in the **`uortos`** folder. Folder **`scratch`** contains auxiliary files for the coregistration process. Once completed, the files can be deleted. If the process is rerun, the auxiliary files are automatically deleted if they are no longer needed and the `clean` parameter is set to `True` in the `parameters.json` file. To make use of the GDAL (https://gdal.org) utilities, here the software developed by OSGEO (https://www.osgeo.org) has been used. It can be downloaded from https://github.com/OSGeo/gdal. In case of using another implementation of GDAL, `ulises_uortos.py` must be modified accordingly. 

To run a demo in folder **`example`** experienced users can run the `example.py` file in a terminal. Alternatively we provide the file `example_notebook.ipynb` to be used in a Jupyter Notebook. In that case, import modules and set the main path of the example:


```python
import ee
import json
import os
import sys
sys.path.insert(0, 'uortos')
import uortos as uortos
pathFldMain = 'example'
```

Set also the data files in the folder **`data`** containing the parameters for the coregistration process (`parameters.json`) and the information related to the downloading of satellite images (`satellites.json`). This second file should not be modified unless you have expertise in the storage of satellite images. These files are read below.


```python
pathJson = os.path.join(pathFldMain, 'data', 'parameters.json')
with open(pathJson, 'r') as f:
    par = json.load(f)
#
pathJson = os.path.join(pathFldMain, 'data', 'satellites.json')
with open(pathJson, 'r') as f:
    par.update(json.load(f))
```

The parameters of these files will be described in the different steps of the process. If any of the values are modified, the process must be restarted from the beginning.

## Image downloading
The images are downloaded through the Google Earth Engine (GEE). The first time images are downloaded, it will be necessary to authenticate in the GEE. You need to link the image downloading to a Google Cloud Project and provide the project name in `eeProject` in the file `parameters.json`.

| Object-name | Description | Example value | 
|:--|:--:|:--:|
| `eeProject` | Google Cloud Project  | `"ee-yourname"` |

Sustitute `"ee-yourname"` for the actual name of the project. Authenticate and initializate GEE:


```python
ee.Authenticate()
ee.Initialize(project=par['eeProject'])
```

The information needed to download the required images must be specified in the `satellites.json` and `parameters.json` files. The satellites, spectral bands and resolutions to be used must be set in the file `satellites.json`. This file contains the following parameters:

| Object-name | Description | Example value | 
|:--|:--:|:--:|
| `satCode2ICode` | Image collection for each satellite code | `"S2": "COPERNICUS/S2_HARMONIZED"` |
| `satCode2RGBBandsCodes` | Indices specifying the RGB-bands for each satellite code | `"S2": ["B4", "B3", "B2"]` |
| `satCode2DownloadBandsCodes` |  Downloaded/output band for each satellite code | `"S2": "all"` |
| `satCode2ResCodes` | List of resolutions codes for each satellite code | `"S2": ["R30", "R10"]` |
| `resCode2Res` | Values (m) for each resolution code | `"R30": 30.0` |

Note that the algorithm processes simultaneously collections of images from different satellites. The different resolutions are processed and clustered separately.

The spatial domain and the temporal period of the images to be co-registered, as well as the amount of clouds in the images and the resolutions to work with, are set in the file `parameters.json`. This file contains the following parameters:

| Object-name | Description | Example value | 
|:--|:--:|:--:|
| `time0` | Initial time period  | `"2024-01-01"` |
| `time1` | Final time period | `"2025-01-01"` |
| `domainCode2Limits` | Coordinates for each domain code | `"CFA": [ 1.94, 2.03, 41.26, 41.30]` |
| `domainCode2Projection` | Filter in the metadata projection for each domain code | `"CFA": "UTM zone 31N"` |
| `cloudThres` | Maximum cloud threshold on the images | `25` |
| `resCodes` | Resolution codes to use | `["R10", "R30"]` |

The image domains are rectangles oriented in the direction of the cardinal directions and defined by the longitude and latitude of each of the edges in the following order `[long-E, long-W, lat-S, lat-N]`. For   each domain code, images whose projection metadata does not contain the characters set in `domainCode2Projection` are discarded. In case `domainCode2Projection` for a given domain codes is `“”`, no filtering is performed. Download the images:


```python
uortos.DownloadImages(pathFldMain, par)
```

## Feature pairing
In a first phase, pairs of common features between the images are obtained using ORB or SIFT method. These pairs are purged by performing a RANSAC (RANdom SAmple Consensus) which evaluates the quality of the projected pairs from image to image using transformation of rotation and translation. The points fit the transformations when the distance between projected pairs is less than `errorsC1`. To ensure the quality of the points used for the transformation, it is also required that the distance of the farthest points verify a minimum ratio, `fractions`, of the horizontal and vertical dimension of the image. In addition, the points must be well distributed across the images. This is set by the parameters `nsPairs` and `nsBoxes`. In some cases it is not necessary to correct the rotation between images, adjust the `rotations` parameter accordingly. The parameters for adjusting this process can be found in the file `parameters.json`:

| Object-name | Description | Example value | 
|:--|:--:|:--:|
| `methodFM` | Feature Matching method | `"SIFT"` |
| `limFM` | Maximum number of images to match each image | `15` |
| `nsFM` | List of maximum number of features | `[1000]` |
| `errorsC1` | Errors in the RANSAC transformations (in pixels) | `[1.0]` |
| `fractions` | Minimum distance between the furthest points (fraction of the image size) | `[0.6]` |
| `nsPairs` | Minumin number of boxes with points | `[5]` |
| `nsBoxes` | Minumun number of boxes in the grid images | `[50]` |
| `rotations` | Option to select whether transformations include rotations | `["free", "null"]` |

Set `limFM` to `0` for avoiding limits. Execute the codes to obtain the first selection of matching pairs:


```python
uortos.FMPairs(pathFldMain, par, par['limFM'] > 0)
uortos.PurgePairs(pathFldMain, par, 'fm')
```

Next, the position of the matched pairs is adjusted by cross-correlation of the images around the points. This process increases the quality of the matches. The error tolerance is set with the parameter `errorsC2` found in the file `parameters.json`:

| Object-name | Description | Example value | 
|:--|:--:|:--:|
| `errorsC2` | Errors in the RANSAC transformations (in pixels)  | `[0.2]` |

Run the code to obtain a final selection of pairs:


```python
uortos.CorPairs(pathFldMain, par)
uortos.PurgePairs(pathFldMain, par, 'cor')
```

## Image clustering
This step takes advantage of the fact that certain images are linked through transformations. Using these connections, the largest cluster of related images is created. Within a cluster, any image can be transformed into another image using these relationships. To ensure high quality connections and filter out outliers, a minimum number of connections per image is required, set by the parameter `conDegrees`. In addition, the quality of the cluster can be further improved by increassing the number of iterations `nSets`. Both parameters are defined in the `parameters.json` file:

| Object-name | Description | Example value | 
|:--|:--:|:--:|
| `conDegrees` | Degree of connections | `[2]` |
| `nSets` | Iterations in the refining process  | `500` |

Execute the code to get the cluster of images:


```python
uortos.Affines(pathFldMain, par)
```

## Outputs
As a result of the process, the metadata information related to the geotransformation of the satellite images (GeoTIFF) is corrected and co-registered PNG images are derived. 
      
Execute the code to get the outputs:


```python
uortos.Output(pathFldMain, par)
```

Output images are located in the folder `output`. For each set of parameters, a folder labelled with the parameter values collects the coregistered images. The label of this folder is coded as follows:

`<parameter_set_01>`=`"nsFM"_"errorsC1"_"fractions"_"nsPairs"_"nsBoxes"_"rotations"_"errorsC2"_"conDegrees"`

Where the numeric values have been written without comma. Each folder contains subfolders named `<res01>_<domain01>`, where `<res01>` is the image resolution code and `<domain01>` is the spatial domain code. Inside these subfolders, you will find the `geopng` and `tif` directories. Additionally, within the `<res01>_<domain01>` folder the georeference information, extracted from the metadata of the most _centered_ image of the output images, is stored in the file `<ref_image>.ref`. The structure of this file is the following:
* One line for each corner pixel
>`pixel-column`, `pixel-row`, `x-coordinate`, `y-coordinate`

If the `<res01>_<domain01>_reference.tif>` image, which can be optionally provided in the `data` folder, is part of the output images, it will be used for `<ref_image>.ref`. Otherwise, the most _centered_ image will be selected. The `scratch` folder contains the original images downloaded from GEE, organized following the structure of `output`. Once the georeferencing process has been completed, it is recommended to delete the `scratch` folder.

## Contact us

Are you experiencing problems? Do you want to give us a comment? Do you need to get in touch with us? Please contact us!

To do so, we ask you to use the [Issues section](https://github.com/Ulises-ICM-UPC/UOrtos/issues) instead of emailing us.

## Contributions

Contributions to this project are welcome. To do a clean pull request, please follow these [guidelines](https://github.com/MarcDiethelm/contributing/blob/master/README.md).

## License

UOrtos is released under a [GPLv3 license](https://github.com/Ulises-ICM-UPC/UOrtos/blob/main/LICENSE). If you use UOrtos in an academic work, please cite:

    @Article{rs12345678,
      AUTHOR = {Simarro, Gonzalo and Calvete, Daniel and Ribas, Francesca and Castillo, Yerai and {Puig-Polo}, Carol},
      TITLE = {UOrtos: A Python Tool for Subpixel co-registration of Landsat and Sentinel 2 Imagery},
      JOURNAL = {Remote Sensing},
      VOLUME = {xx},
      YEAR = {2025},
      NUMBER = {xx},
      ARTICLE-NUMBER = {xxxx},
      URL = {https://www.mdpi.com/2072-4292/xx/xx/xxxx},
      ISSN = {2072-4292},
      DOI = {10.3390/rs12345678}
      }

    @Online{ulisesortos, 
      author = {Simarro, Gonzalo and Calvete, Daniel},
      title = {UOrtos},
      year = 2025,
      url = {https://github.com/Ulises-ICM-UPC/UOrtos}
      }
