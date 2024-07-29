# UOrtOs

`UOrtos` is an open source software to obtain a set of co-registered satellite images with high alignment quality and outlier-free.

### Description
The algorithm for co-registration of satellite images assumes that images can be realigned through an isometric transformation. The image alignment is achieved by minimising the distances between common image features. Clustering techniques are used to obtain a high quality set of aligned images without outliers. As a result of the process, the metadata information related to the geotransformation of the satellite images (GeoTIFF) is corrected and co-registered PNG images are derived. The development of this software is aimed for obtaining images to accurately determine shoreline displacements. In general, it can be used for the precise detection of features at different times. Details about the algorithm and methodology are described in
> *Simarro, G.; Calvete, D.; Puig, C.; Ribas, F. UOrtos: Accurate Co-Registration of Satellite Images. Remote Sens. 2024, xx, xxxx. https://doi.org/10.3390/rs12345678*

The co-registration process consists of the following steps:

1. [Image downloading](#image-download)
2. [Feature pairing](#feature-pairing)
3. [Image clustering](#clustering)
4. [Outputs](#output)

### Requirements and project structure
To run the software it is necessary to have Python (3.9) and install the following dependencies:
- cv2 (4.2.0)?
- numpy (1.19.5)?
- scipy (1.6.0)?
- osgeo (3.9.1)????
- ee ()????

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
  * **`output`**
    * **`<parameter_set_01>`**
      * **`<code01>_<domain01>`**
        * `<ref_image>.log`
        * `<ref_image>.ref`
        * **`geopng`**
          * `<image01>.png`
          * . . .
        * **`tif`**
          * **`<image01>`**
            * `<image01>_<band01>.png`
            * . . .
          * . . .
      * . . .
    * . . .
  * **`scratch`**
    * . . .

The local modules of `UCalib` are located in the **`ucalib`** folder. Folder **`scratch`** contains auxiliary files for the coregistration process. Once completed, the files can be deleted. In case of rerunning the process the auxiliary files are automatically purged if no longer needed and paramater `clean` is set to `True` in the parameter file `parameters.json`. To make use of the GDAL (https://gdal.org) utilities, here the software developed by OSGEO (https://www.osgeo.org) has been used. It can be downloaded from https://github.com/OSGeo/gdal. In case of using another implementation of GDAL, `ulises_uortos.py` must be modified accordingly. 

To run a demo in folder **`example`** experienced users can run the `example.py` file in a terminal. Alternatively we provide the file `example_notebook.ipynb` to be used in a Jupyter Notebook. In that case, import modules and set the main path of the example:


```python
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
The images are downloaded through the Google Earth Engine (GEE). The satellites, spectral bands and resolutions to be used must be set in the file `satellites.json`. This file contains the following parameters:

| Object-name | Description | Value sample | 
|:--|:--:|:--:|
| `satCode2ICode` | List of codes and corresponding satellite image collection | `"S2": "COPERNICUS/S2_HARMONIZED"` |
| `satCode2RGBBandsCodes` | List of indices specifying the bands to select for each collection | `"S2": ["B4", "B3", "B2"]` |
| `satCode2DownloadBandsCodes` |  List of new names for the output bands | `"S2": "RGB"` |
| `satCode2ResCodes` | List of resolutions for each collection  | `"S2": ["R30", "R10"]` |
| `resCode2Res` | Value of the resolution sacale (m) | `"R30": 30.0` |

Note that the algorithm processes simultaneously collections of images from different satellites. The different resolutions are processed and clustered separately.

The spatial domain and the temporal period of the images to be co-registered, as well as the amount of clouds in the images and the resolutions to work with, are set in the file `parameters.json`. This file contains the following parameters:

| Object-name | Description | Value sample | 
|:--|:--:|:--:|
| `time0` | Initial time period  | `"2022-01-01"` |
| `time1` | Final time period | `"2022-07-01"` |
| `domainCode2Limits` | List of codes and corresponding coordenates | `"NAR": [ 151.280, 151.320, -33.750, -33.690]` |
| `cloudThres` | Maximum cloud threshold on the images | `5` |
| `resCodes` | Codes of the resolutions to use | `["R10", "R30"]` |

The image domains are rectangles oriented in the direction of the cardinal directions and defined by the longitude and latitude of each of the edges in the following order `[long-E, long-W, lat-S, lat-N]`.

The first time images are downloaded, it will be necessary to authenticate in the GEE. To run the code to generate the meshes:


```python
uortos.DownloadImages(pathFldMain, par)
```

## Feature pairing
In a first phase, pairs of common features between the images are obtained using ORB. These pairs are purged by performing a RANSAC (RANdom SAmple Consensus) which evaluates the quality of the projected pairs from image to image using transformation of rotation and translation. The points fit the transformations when the distance between projected pairs is less than `errorC1s`. To ensure the quality of the points used for the transformation, it is also required that the diatance of the furthest points verifies a minimum ratio, `fractions`, of the diagonal of the image. In addition, the points must be well distributed across the images. This is set by the parameters `nFillBoxes` and `nBoxes`. In some cases it is not necessary to correct the rotation between images, adjust the `rotations` parameter accordingly. The parameters for adjusting this process can be found in the file `parameters.json`:

| Object-name | Description | Value sample | 
|:--|:--:|:--:|
| `limORB` | _ni puta idea_  | `0` |
| `ORBPairs` | Number of ORB features | `[1000]` |
| `errorC1s` | Errors in the RANSAC transformations (in pixels) | `[1.0]` |
| `fractions` | Minimum distance between the furthest points (fraction of the image size) | `[0.5]` |
| `nFillBoxes` | Minumin number of boxes with points | `[5]` |
| `nBoxes` | Minumun number of boxes in the grid images | `[50]` |
| `rotations` | Option to select whether transformations include rotations | `["free", "null"]` |

Other parameters that can be adjusted are in the file `satellites.json`:

| Object-name | Description | Value sample | 
|:--|:--:|:--:|
| `resCode2ScaleOfORB` | _multiplicador para orb?_ | `"R30": 3` |
| `resCode2SemiSize` | _reductor para comparar?_  | `"R30": 5` |

Execute the codes to obtain the first selection of ORB pairs:


```python
uortos.ORBPairs(pathFldMain, par, par['limORB'] > 0)
uortos.CleanPairs(pathFldMain, par, '02', '03')
```

Next, the position of the ORB pairs is adjusted by cross-correlation of the images surrounding the points. This process increases the quality of the ORB points. The error tolerance is set with the parameter `errorC2s` found in the file `parameters.json`:

| Object-name | Description | Value sample | 
|:--|:--:|:--:|
| `errorC2s` | Errors in the RANSAC transformations (in pixels)  | `[0.2]` |

Run the code to obtain a final selection of pairs:


```python
uortos.CorPairs(pathFldMain, par)
uortos.CleanPairs(pathFldMain, par, '04', '05')
```

Together with the resulting set of ORB pairs between images, rotations and translations are obtained, which allow to transform the position and orientation from image to image. 

## Image clustering
This step make use of the fact that certain images can be connected to other images by transformations. Through these connections the largest cluster of connected images is constructed. Using the transformations of the connected images all images in a cluster can be transformed to any one of them. To guarantee that connections between images in the cluster are of high quality and to avoid the presence of outlayers, a minimum number of connections can be required for each image in the cluster. This is the cluster degree set with the `conDegrees` parameter. The quality of the final cluster can also be improved by performing `nIter` iterations of this process. These two parameters can be set with the file `parameters.json`:

| Object-name | Description | Value sample | 
|:--|:--:|:--:|
| `conDegrees` | Degree of connections | `[2]` |
| `nIter` | Iterations in the clustering process  | `6` |

Execute the code to get the cluster of images:


```python
uortos.Affines(pathFldMain, par, '05', '06')
```

## Outputs
As a result of the process, the metadata information related to the geotransformation of the satellite images (GeoTIFF) is corrected and co-registered PNG images are derived. 
      
Execute the code to get the outputs:


```python
uortos.Output(pathFldMain, par, '06')
```

Output images are located in the folder `output`. For each set of parameters a folder labelled with the parameter values collects the coregistered images. The label of this folder is coded as follows:

`<parameter_set_01>`=`"nORBs"_"errorC1s"_"fractions"_"nFillBoxes"_"nBoxes"_"rotations"_"errorC2s"_"conDegrees"`

Where the muneric values have been written without comma. Within each folder are other folders labelled with the resolution of the images and the spatial domain code.


## Contact us

Are you experiencing problems? Do you want to give us a comment? Do you need to get in touch with us? Please contact us!

To do so, we ask you to use the [Issues section](https://github.com/Ulises-ICM-UPC/UOrtos/issues) instead of emailing us.

## Contributions

Contributions to this project are welcome. To do a clean pull request, please follow these [guidelines](https://github.com/MarcDiethelm/contributing/blob/master/README.md).

## License

UCalib is released under a [GPLv3 license](https://github.com/Ulises-ICM-UPC/UOrtos/blob/main/LICENSE). If you use UCalib in an academic work, please cite:

    @Article{rs13142795,
      AUTHOR = {Simarro, Gonzalo and Calvete, Daniel and Souto, Paola},
      TITLE = {UCalib: Cameras Autocalibration on Coastal Video Monitoring Systems},
      JOURNAL = {Remote Sensing},
      VOLUME = {13},
      YEAR = {2021},
      NUMBER = {14},
      ARTICLE-NUMBER = {2795},
      URL = {https://www.mdpi.com/2072-4292/13/14/2795},
      ISSN = {2072-4292},
      DOI = {10.3390/rs13142795}
      }

    @Online{ulisesdrone, 
      author = {Simarro, Gonzalo and Calvete, Daniel},
      title = {UCalib},
      year = 2021,
      url = {https://github.com/Ulises-ICM-UPC/UCalib}
      }
