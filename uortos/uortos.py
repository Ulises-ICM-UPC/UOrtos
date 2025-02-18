# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~ by Gonzalo Simarro and Daniel Calvete in 2024
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
import cv2  # type: ignore
import ee  # type: ignore
import itertools
import json
import numpy as np  # type: ignore
import os
import random
import requests  # type: ignore
import shutil
import sys
#
import ulises_uortos as ulises
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def DownloadImages(pathFldMain, par):  # generates 00 and 01
    #
    # obtain pathFldScratch and clean
    print('  ... preliminary checkings...')
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    ulises.CleanScratch240227(pathFldMain, par)
    #
    # obtain and check yyyymmdd0 and yyyymmdd1
    yyyymmdd0, yyyymmdd1 = [par[item].replace('-', '') for item in ['time0', 'time1']]
    if yyyymmdd1 <= yyyymmdd0:
        print('    ... make sure time1 > time0 in parameters.json')
        sys.exit()
    #
    # download images
    viewCodes = ulises.Par2ViewCodes240227(par)
    for viewCode in viewCodes:
        #
        # obtain resCode, domainCode and res
        resCode, domainCode = viewCode.split('_')
        res = par['resCode2Res'][resCode]
        #
        # obtain AoI
        E0, E1, N0, N1 = par['domainCode2Limits'][domainCode]  # E = lon, N = lat
        if not (E1 > E0 and N1 > N0):
            print('    ... make sure that E1 > E0 and N1 > N0 in domainCode2Limits in parameter.json')
            sys.exit()
        AoI = [[[E0, N0], [E0, N1], [E1, N1], [E1, N0], [E0, N0]]]
        #
        # obtain paths of folders and inform
        pathFld00View = os.path.join(pathFldScratch, '00_zip_tif', viewCode)
        pathFld00ViewZip = os.path.join(pathFld00View, 'zip')
        pathFld00ViewTif = os.path.join(pathFld00View, 'tif')
        pathFld01View = os.path.join(pathFldScratch, '01_images', viewCode)
        print('  ... working on {:} (domain {:}, resolution {:}m); '.format(viewCode, domainCode, int(res)), end='', flush=True)
        #
        # check internet
        try:
            # try:
            #     ee.Initialize()
            # except Exception:
            ee.Initialize(project=par['eeProject'])
            internet = True
        except Exception:
            print('... authentication or initialization not successful in GEE, images will not be downloaded', end='', flush=True)
            internet = False
        #
        # download pathsZips
        if internet:
            for satCode in sorted(par['satCode2ICode']):
                #
                # disregard
                if resCode not in par['satCode2ResCodes'][satCode]:
                    continue
                #
                # try:
                #     ee.Initialize()
                # except Exception:
                ee.Initialize(project=par['eeProject'])
                #
                # obtain lResults
                eeCollection = ee.ImageCollection(par['satCode2ICode'][satCode])
                eeCollection = eeCollection.filterBounds(ee.Geometry.Polygon(AoI)).filterDate(par['time0'], par['time1'])
                lResults = eeCollection.getInfo().get('features')
                if satCode.startswith('L'):  # WATCH OUT
                    lResults = [item for item in lResults if item['properties']['CLOUD_COVER'] < par['cloudThres']]
                elif satCode.startswith('S'):  # WATCH OUT
                    lResults = [item for item in lResults if item['properties']['CLOUDY_PIXEL_PERCENTAGE'] < par['cloudThres']]
                else:  # WATCH OUT
                    continue
                #
                # download zip files (with RGB or all bands)
                for result in lResults:  # random.sample(lResults, len(lResults)):
                    try:
                        # obtain pathZip
                        name = '{:}_{:}'.format(satCode, result['id'].split('/')[-1])  # IMP*; '{:}_'.format(satCode) in zip-name
                        pathZip = os.path.join(pathFld00ViewZip, '{:}.zip'.format(name))
                        if not os.path.exists(pathZip):  # IMP*; does not overwrite
                            # write pathZip
                            img = ee.Image(result['id'])
                            if par['satCode2DownloadBandsCodes'][satCode] == 'all':
                                bandsCodes = img.bandNames().getInfo()
                            elif par['satCode2DownloadBandsCodes'][satCode] == 'RGB':
                                bandsCodes = par['satCode2RGBBandsCodes'][satCode]
                            else:
                                assert False
                            dTMP = {'image': img, 'bands': bandsCodes, 'region': AoI, 'scale': res, 'fileFormat': 'GeoTIFF'}
                            url = ee.data.makeDownloadUrl(ee.data.getDownloadId(dTMP))
                            response = requests.get(url)
                            os.makedirs(os.path.dirname(pathZip), exist_ok=True)
                            with open(pathZip, 'wb') as fd:
                                fd.write(response.content)
                            ulises.WriteTimeStampLog(pathFld00View)
                    except Exception:
                        pass
        #
        # inform
        fnsZips = ulises.PathFld2Fns240227(pathFld00ViewZip, ext='.zip')
        print('{:} zip-files available; '.format(len(fnsZips)), end='', flush=True)
        #
        # obtain tif for zips; checks, and creates no pathFldBands if it does not succeed
        for fnZip in fnsZips:  # reference does not show up in zips
            #
            # obtain pathZip, name, yyyymmdd and disregard
            pathZip = os.path.join(pathFld00ViewZip, fnZip)
            name = os.path.splitext(fnZip)[0]
            yyyymmdd = ulises.Name2YYYYMMDD240227(name)
            if not (int(yyyymmdd0) <= int(yyyymmdd) <= int(yyyymmdd1)):
                continue
            #
            # write pathFldBands
            pathFldBands = os.path.join(pathFld00ViewTif, name)
            if not os.path.exists(pathFldBands):  # IMP*; does not overwrite
                wellUnzipped = ulises.PathZip2PathFldBands240227(pathZip, pathFldBands, par)  # IMP*
                if wellUnzipped:
                    ulises.WriteTimeStampLog(pathFld00View)
                else:
                    assert not os.path.exists(pathFldBands)  # WATCH OUT; could be rmtree if exists
            if not os.path.exists(pathFldBands):
                continue
            #
            # check projection
            desiredProjection = par['domainCode2Projection'][viewCode.split('_')[1]]
            if desiredProjection.strip() == '':
                continue
            for fnTif in os.listdir(pathFldBands):  # are to be .tif
                assert fnTif.endswith('.tif')
                metadata = ulises.PathTif2Metadata240227(os.path.join(pathFldBands, fnTif))
                if desiredProjection not in metadata['projection']:
                    shutil.rmtree(pathFldBands)
                    ulises.WriteTimeStampLog(pathFld00View)
                    break
        #
        # obtain tif for reference; checks, and no pathFldBands if it does not succeed
        pathAuxRef = os.path.join(pathFldMain, 'data', '{:}_reference.tif'.format(viewCode))  # IMP*
        if os.path.exists(pathAuxRef):
            pathFldBands = os.path.join(pathFld00ViewTif, '00_reference')  # IMP*
            wellUnzipped = ulises.PathGivenRefTif2PathFldBands240227(pathAuxRef, pathFldBands)
            if wellUnzipped:
                ulises.WriteTimeStampLog(pathFld00View)
        #
        # obtain png for tifs
        names = ulises.PathFld2Flds240227(pathFld00ViewTif)
        for name in names:
            pathFldBands = os.path.join(pathFld00ViewTif, name)
            pathPng = os.path.join(pathFld01View, '{:}.png'.format(name))
            if os.path.exists(pathPng):  # WATCH OUT: does not overwrite
                continue
            img = ulises.PathFldBands2ImgRGBRaw240227(pathFldBands, par['satCode2RGBBandsCodes'][name.split('_')[0]])
            if img is None:
                continue
            os.makedirs(os.path.dirname(pathPng), exist_ok=True)
            cv2.imwrite(pathPng, img)
            ulises.WriteTimeStampLog(pathFld01View)
        #
        fnsPngs = ulises.PathFld2Fns240227(pathFld01View, ext='.png')
        print('n = {:} png-files to analyze'.format(len(fnsPngs)))
    #
    return None
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def FMPairs(pathFldMain, par, isLimFM):  # generates 02 from 01
    #
    # set hidden parameters (used to be in par)
    methodFM = par['methodFM']
    if methodFM == 'ORB':
        resCode2ScaFM = {'R30': 3, 'R10': 2}
    elif methodFM == 'SIFT':
        resCode2ScaFM = {'R30': 2, 'R10': 1}  # enough
    else:
        assert False
    #
    # obtain pathFldScratch and clean
    print('  ... preliminary checkings...')
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    ulises.CleanScratch240227(pathFldMain, par)
    #
    # obtain viewCodes and flds02
    viewCodes = ulises.Par2ViewCodes240227(par)
    flds02 = ulises.ParAndXX2FldsXX240227(par, '02')
    #
    # run for viewCodes x flds02
    for viewCode, fld02 in itertools.product(viewCodes, flds02):
        #
        # obtain paths
        pathFld01View = os.path.join(pathFldScratch, '01_images', viewCode)
        pathFld02View = os.path.join(pathFldScratch, fld02, viewCode)
        if not os.path.exists(pathFld01View):
            ulises.RmTreeIfExists(pathFld02View)
            continue
        #
        # obtain nFM and scaFM
        dInfoFld02 = ulises.FldXX2DInfo240227(fld02)
        nFM = dInfoFld02['nFM']
        scaFM = resCode2ScaFM[viewCode.split('_')[0]]
        #
        # obtain fnsImgs and inform
        fnsImgs = ulises.PathFld2Fns240227(pathFld01View, ext='.png')
        names = [os.path.splitext(item)[0] for item in fnsImgs]
        print('  ... working on {:} to get {:}: n = {:04d} initial png-files; '.format(viewCode, fld02, len(names)), end='', flush=True)
        #
        # manage isLimFM
        if not isLimFM:  # manage if the (exhaustive) work is already done
            # obtain namesDone
            if os.path.exists(pathFld02View):
                namesDone = ulises.PathFld2Flds240227(pathFld02View)
            else:
                namesDone = []
            assert set(namesDone) <= set(names)
            if len(namesDone) == len(names):
                nOfNpzs02 = len([os.path.join(x[0], fn) for x in os.walk(pathFld02View) for fn in x[2] if fn.endswith('.npz')])
                N = ulises.ViewCode2N240227(pathFldScratch, viewCode)
                print('{:06d} npz-files: {:5.1f}% of N = n * (n - 1) / 2'.format(nOfNpzs02, nOfNpzs02 / N * 100))
                continue
        #
        # make all folders
        for name in names:
            os.makedirs(os.path.join(pathFld02View, name), exist_ok=True)
        ulises.WriteTimeStampLog(pathFld02View)
        #
        # disregard
        if len(names) <= 1:
            nOfNpzs02, N = 0, 1  # fake
            print('{:06d} npz-files: {:5.1f}% of N = n * (n - 1) / 2'.format(nOfNpzs02, nOfNpzs02 / N * 100))
            return None
        #
        # obtain precomputations
        ncsOri, nrsOri, ncsSca, nrsSca, kpssSca, dessSca = ulises.FMPre240227(pathFld01View, fnsImgs, methodFM, scaFM, nFM)
        assert all(np.isclose(ncsSca[pos] / ncsOri[pos], scaFM) for pos in range(len(ncsOri)) if kpssSca[pos] is not None)
        assert all(np.isclose(nrsSca[pos] / nrsOri[pos], scaFM) for pos in range(len(nrsOri)) if kpssSca[pos] is not None)
        #
        # obtain npz-files
        if not isLimFM:
            # run exhaustively pos0 < pos1
            for pos0 in range(len(names)):
                # obtain name0, pathFld02ViewName0 and disregard
                name0 = names[pos0]
                pathFld02ViewName0 = os.path.join(pathFld02View, name0)
                if kpssSca[pos0] is None:
                    continue
                for pos1 in range(len(names)):
                    # check and disregard
                    if not pos0 < pos1:
                        continue
                    # obtain name1 and pathNpz02 and diregard
                    name1 = names[pos1]
                    pathNpz02 = os.path.join(pathFld02ViewName0, '{:}.npz'.format(name1))
                    if kpssSca[pos1] is None or os.path.exists(pathNpz02) or all(item in namesDone for item in [name0, name1]):
                        continue
                    # obtain and write pathNpz02 (timestamps)
                    kps0Sca, des0Sca, kps1Sca, des1Sca = kpssSca[pos0], dessSca[pos0], kpssSca[pos1], dessSca[pos1]
                    nc0Ori, nr0Ori, nc1Ori, nr1Ori = ncsOri[pos0], nrsOri[pos0], ncsOri[pos1], nrsOri[pos1]
                    ulises.WritePathNpz02240227(pathNpz02, methodFM, nc0Ori, nr0Ori, nc1Ori, nr1Ori, kps0Sca, des0Sca, kps1Sca, des1Sca, scaFM)
        else:
            # run to obtain a subset of useful npzs
            for _ in range(10):  # WATCH OUT; could be while True
                workDone0 = ulises.MFMConnsAddPairs240227(pathFld02View, methodFM, names, ncsOri, nrsOri, kpssSca, dessSca, scaFM, par['limFM'])
                workDone1 = ulises.MFMConnsDeletePairs240227(pathFld02View, names, par['limFM'])
                if not workDone0 and not workDone1:
                    break
            ulises.WriteTimeStampLog(pathFld02View)
        #
        nOfNpzs02 = len([os.path.join(x[0], fn) for x in os.walk(pathFld02View) for fn in x[2] if fn.endswith('.npz')])
        N = ulises.ViewCode2N240227(pathFldScratch, viewCode)
        print('{:06d} npz-files: {:5.1f}% of N = n * (n - 1) / 2'.format(nOfNpzs02, nOfNpzs02 / N * 100))
        #
    return None
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def PurgePairs(pathFldMain, par, case):  # generates yy = 03 from xx = 02 and yy = 05 from xx = 04
    #
    # manage case (xx < yy)
    if case == 'fm':
        xx, yy = '02', '03'
    elif case == 'cor':
        xx, yy = '04', '05'
    else:
        assert False
    #
    # obtain pathFldScratch and clean
    print('  ... preliminary checkings...')
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    ulises.CleanScratch240227(pathFldMain, par)
    #
    # obtain viewCodes and fldsYY (to do)
    viewCodes = ulises.Par2ViewCodes240227(par)
    fldsYY = ulises.ParAndXX2FldsXX240227(par, yy)
    #
    # run for viewCodes x fldsYY
    for viewCode, fldYY in itertools.product(viewCodes, fldsYY):
        #
        # obtain paths
        fldXX = ulises.FldXX2FldYY240227(fldYY, xx, options={})
        pathFldXXView = os.path.join(pathFldScratch, fldXX, viewCode)
        pathFldYYView = os.path.join(pathFldScratch, fldYY, viewCode)
        if not os.path.exists(pathFldXXView):
            continue
        #
        # obtain dInfoFldYY and relevant values
        dInfoFldYY = ulises.FldXX2DInfo240227(fldYY)
        fraction, nPairs, nBoxes, rotation = [dInfoFldYY[item] for item in ['fraction', 'nPairs', 'nBoxes', 'rotation']]
        if case == 'fm':
            errorC = dInfoFldYY['errorC1']
        elif case == 'cor':
            errorC = dInfoFldYY['errorC2']
        else:
            assert False
        #
        # obtain pathsNpzsXX and inform
        if os.path.exists(pathFldXXView):
            pathsNpzsXX = [os.path.join(x[0], fn) for x in os.walk(pathFldXXView) for fn in x[2] if fn.endswith('.npz')]
        else:
            pathsNpzsXX = []
        print('  ... working on {:} to get {:}: {:06d} initial npz-files; '.format(viewCode, fldYY, len(pathsNpzsXX)), end='', flush=True)
        #
        # run through pathNpzXX
        for pathNpzXX in random.sample(pathsNpzsXX, len(pathsNpzsXX)):
            #
            # obtain dInfoXX
            dInfoXX = ulises.PathNpzOfPairs2DInfo240227(pathNpzXX, options={'nBoxes': nBoxes, 'par': par})
            #
            # load useful variables from dInfoXX
            nc0, nr0, nc1, nr1, cs0, rs0, cs1, rs1, ers = [dInfoXX[item] for item in ['nc0', 'nr0', 'nc1', 'nr1', 'cs0', 'rs0', 'cs1', 'rs1', 'ers']]
            nOfCBands, nOfRBands = [dInfoXX[item] for item in ['nOfCBands', 'nOfRBands']]
            #
            # check and disregard
            if dInfoXX['nOfFBoxes'] < nPairs:
                continue
            if not (np.max(cs0) - np.min(cs0) >= fraction * nc0 and np.max(rs0) - np.min(rs0) >= fraction * nr0):
                continue
            #
            # obtain pathNpzYY and nOfFBoxes0, check and disregard
            pathNpzYY = os.path.join(pathFldYYView, dInfoXX['name0'], '{:}.npz'.format(dInfoXX['name1']))
            if os.path.exists(pathNpzYY):  # check and continue
                dInfoYY = ulises.PathNpzOfPairs2DInfo240227(pathNpzYY, options={'nBoxes': nBoxes, 'par': par})
                if dInfoYY['nOfFBoxes'] < nPairs or dInfoYY['datetime'] <= dInfoXX['datetime']:
                    os.remove(pathNpzYY)
                    nOfFBoxes0 = 0  # IMP*
                else:
                    nOfFBoxes0 = dInfoYY['nOfFBoxes']
                    continue  # WATCH OUT; continue or not (also possible)
            else:
                nOfFBoxes0 = 0  # IMP*
            #
            # obtain RANSAC (PurgePixelsRANSAC240227, if not providing null, ensures fraction but not nOfFBoxes) and check nOfFBoxes
            cs0, rs0, cs1, rs1, ers = ulises.PurgePixelsRANSAC240227(1, cs0, rs0, cs1, rs1, ers, nc0, nr0, fraction, errorC, rotation, nOfCBands, nOfRBands, nOfIter=10)
            if not ulises.EnoughImprovingFBoxes240227(nc0, nr0, nOfCBands, nOfRBands, cs0, rs0, nPairs, nOfFBoxes0):
                continue
            #
            # write pathNpzYY
            os.makedirs(os.path.dirname(pathNpzYY), exist_ok=True)
            np.savez(pathNpzYY, nc0=nc0, nr0=nr0, nc1=nc1, nr1=nr1, cs0=cs0, rs0=rs0, cs1=cs1, rs1=rs1, ers=ers)
            ulises.WriteTimeStampLog(pathFldYYView)
            #
        #
        # obtain nOfNpzsYY, N and inform
        textAux = {'03': 'm1', '05': 'm2'}[yy]
        nOfNpzYY = len([os.path.join(x[0], fn) for x in os.walk(pathFldYYView) for fn in x[2] if fn.endswith('.npz')])
        N = ulises.ViewCode2N240227(pathFldScratch, viewCode)
        if np.isclose(N, 0):
            print('{:06d} npz-files: {:} = {:5.1f}% of N'.format(0, textAux, 0))
        else:
            print('{:06d} npz-files: {:} = {:5.1f}% of N'.format(nOfNpzYY, textAux, nOfNpzYY / N * 100))
    #
    return None
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def CorPairs(pathFldMain, par):  # generates 04 from 03
    #
    # set hidden parameters
    resCode2SemiSize = {'R30': 5, 'R10': 10}
    #
    # obtain pathFldScratch and clean
    print('  ... preliminary checkings...')
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    ulises.CleanScratch240227(pathFldMain, par)
    #
    # obtain viewCodes and flds04
    viewCodes = ulises.Par2ViewCodes240227(par)
    flds04 = ulises.ParAndXX2FldsXX240227(par, '04')
    #
    # run for viewCodes x flds04
    for viewCode, fld04 in itertools.product(viewCodes, flds04):
        #
        # obtain paths
        fld03 = ulises.FldXX2FldYY240227(fld04, '03', options={})  # options={} since '03' < '04'
        pathFld01View = os.path.join(pathFldScratch, '01_images', viewCode)
        pathFld03View = os.path.join(pathFldScratch, fld03, viewCode)
        pathFld04View = os.path.join(pathFldScratch, fld04, viewCode)
        if not os.path.exists(pathFld03View):
            continue
        #
        # obtain pathsNpzs03 and inform
        if os.path.exists(pathFld03View):
            pathsNpzs03 = [os.path.join(x[0], fn) for x in os.walk(pathFld03View) for fn in x[2] if fn.endswith('.npz')]
        else:
            pathsNpzs03 = []
        print('  ... working on {:} to get {:}: {:06d} initial npz-files; '.format(viewCode, fld04, len(pathsNpzs03)), end='', flush=True)
        #
        # run through all files
        for pathNpz03 in pathsNpzs03:  # random.sample(pathsNpzs03, len(pathsNpzs03)) 
            #
            # obtain dInfo03
            dInfo03 = ulises.PathNpzOfPairs2DInfo240227(pathNpz03, options={})
            #
            # obtain pathNpz04, check and disregard
            pathNpz04 = os.path.join(pathFld04View, dInfo03['name0'], '{:}.npz'.format(dInfo03['name1']))
            if os.path.exists(pathNpz04):
                dInfo04 = ulises.PathNpzOfPairs2DInfo240227(pathNpz04, options={})
                if dInfo04['datetime'] <= dInfo03['datetime']:
                    os.remove(pathNpz04)
                else:
                    continue
            else:
                pass
            #
            # load useful variables from dInfo03
            nc0, nr0, nc1, nr1, cs0, rs0, cs1, rs1 = [dInfo03[item] for item in ['nc0', 'nr0', 'nc1', 'nr1', 'cs0', 'rs0', 'cs1', 'rs1']]
            sameSize = nc0 == nc1 and nr0 == nr1
            #
            # obtain img0 and img1
            img0 = cv2.imread(os.path.join(pathFld01View, '{:}.png'.format(dInfo03['name0'])))
            img1 = cv2.imread(os.path.join(pathFld01View, '{:}.png'.format(dInfo03['name1'])))
            #
            # obtain cs0, rs0, cs1 and rs1
            semiSize = resCode2SemiSize[viewCode.split('_')[0]]
            cs0, rs0, cs1, rs1 = ulises.ImprovePairsThroughCorrelation240227(cs0, rs0, cs1, rs1, img0, img1, sameSize, semiSize=semiSize, scale0=5)  # WATCH OUT
            #
            # write pathNpz04
            os.makedirs(os.path.dirname(pathNpz04), exist_ok=True)
            np.savez(pathNpz04, nc0=nc0, nr0=nr0, nc1=nc1, nr1=nr1, cs0=cs0, rs0=rs0, cs1=cs1, rs1=rs1, ers=np.random.random(len(cs0)))  # WATCH OUT; ers
            ulises.WriteTimeStampLog(pathFld04View)
            #
        #
        # obtain nOfNpzs04, N and inform
        nOfNpzs04 = len([os.path.join(x[0], fn) for x in os.walk(pathFld04View) for fn in x[2] if fn.endswith('.npz')])
        N = ulises.ViewCode2N240227(pathFldScratch, viewCode)
        print('{:06d} npz-files (as above, {:5.1f}% of N)'.format(nOfNpzs04, nOfNpzs04 / N * 100))
        #
    return None
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def Affines(pathFldMain, par, case='default'):
    #
    # manage case
    if case == 'default':
        xx, yy = '05', '06'
    elif case == 'bypass':
        xx, yy = '15', '16'
    else:
        assert False
    #
    # obtain pathFldScratch and clean yy
    print('  ... preliminary checkings...')
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    ulises.CleanScratch240227(pathFldMain, par)
    #
    # obtain viewCodes and flds06
    viewCodes = ulises.Par2ViewCodes240227(par)
    flds06 = ulises.ParAndXX2FldsXX240227(par, yy)  # since yy = 06 typically
    #
    # run for viewCodes x flds06
    for viewCode, fld06 in itertools.product(viewCodes, flds06):
        #
        # obtain paths
        fld05 = ulises.FldXX2FldYY240227(fld06, xx, options={})  # since xx = 05 typically
        pathFld01View = os.path.join(pathFldScratch, '01_images', viewCode)
        pathFld05View = os.path.join(pathFldScratch, fld05, viewCode)
        pathFld06View = os.path.join(pathFldScratch, fld06, viewCode)
        if not os.path.exists(pathFld05View):
            continue
        #
        # obtain dInfoFld06 and relevant values
        dInfoFld06 = ulises.FldXX2DInfo240227(fld06)
        fraction, nPairs, rotation, errorC2, conDegree = [dInfoFld06[item] for item in ['fraction', 'nPairs', 'rotation', 'errorC2', 'conDegree']]
        #
        # obtain pathsNpzs05, disregard and inform
        if os.path.exists(pathFld05View):
            pathsNpzs05 = [os.path.join(x[0], fn) for x in os.walk(pathFld05View) for fn in x[2] if fn.endswith('.npz')]
        else:
            pathsNpzs05 = []
        if len(pathsNpzs05) == 0:
            continue
        print('  ... working on {:} to get {:}: {:06d} initial npz-files; '.format(viewCode, fld06, len(pathsNpzs05)), end='', flush=True)
        #
        # obtain size of the images and precomputations
        fnsPng = ulises.PathFld2Fns240227(pathFld01View, ext='.png')
        nc4SatCode, nr4SatCode = [{} for _ in range(2)]
        for fnPng in fnsPng:
            satCode = fnPng.split('_')[0]  # can be 00
            if fnPng.split('_')[0] in nc4SatCode.keys():
                continue
            nrImg, ncImg = cv2.imread(os.path.join(pathFld01View, fnPng)).shape[0:2]
            nc4SatCode[satCode] = ncImg
            nr4SatCode[satCode] = nrImg
        #
        # obtain names and inform ('global' poss are referred to names hereafter)
        names = ulises.ViewCode2Names240227(pathFldScratch, viewCode)
        #
        # obtain mConns and dPixelsConns (pos0 < pos1 in dPixelsConns)
        mConns, dPixelsConns = ulises.MConnsAndDPixelsConns240227(names, pathFld05View)
        #
        # obtain groups  (group = list of global poss) from mConns
        groups = ulises.MConns2Groups240227(mConns, conDegree)  # has repetitivity
        assert sorted(list([item0 for item1 in groups for item0 in item1])) == list(range(len(names)))  # check for readability
        sizesOfGroups = [len(item) for item in groups]
        print(' s1 = {:5.1f}% of n;'.format(max(sizesOfGroups) / len(names) * 100), end='', flush=True)
        #
        # obtain date0
        date0 = ulises.ReadTimeStampLog(pathFld06View)
        #
        # analyse each group (group = list of global poss)
        maxSizeOfSubgroup = 0
        for posGroup, group in enumerate(groups):
            #
            # obtain and make pathFldGroup
            assert group == sorted(group)  # avoidable check
            groupCode = ulises.PosGroup2GroupCode240227(posGroup)
            pathFldGroup = os.path.join(pathFld06View, groupCode)
            os.makedirs(pathFldGroup, exist_ok=True)
            #
            # obtain namesGroup, mConnsGroup, nameRefGroup (not very relevant) and write names.txt, connections_FM.json and *.ref (timestamps)
            namesGroup, mConnsGroup, nameRefGroup = ulises.WriteGroups240227(pathFldGroup, names, mConns, group)
            assert len(namesGroup) == mConnsGroup.shape[0] == mConnsGroup.shape[1]
            #
            # scape
            if len(group) <= conDegree+1:  # VERY IMP*; scape
                continue
            if not len(ulises.MConns2Groups240227(mConnsGroup, conDegree)) == 1:  # VERY IMP*; scape
                continue
            #
            # obtain subgroups for non-trivial groups
            nIter = int((par['nSets'] + 1) ** (1 / 3)) + 1
            for _ in range(nIter):
                #
                # obtain a lighter valid mConnsGroup (only to find DRSSTsRefIsOK240227 pseudowalks)
                mConnsGroupLight = ulises.MConnsLight240227(mConnsGroup, conDegree)
                if mConnsGroupLight is None:
                    continue
                #
                # run through randomPWalkAsPairs
                for _ in range(nIter):
                    #
                    # obtain random pseudowalk with (group) local numbering
                    randomPWalkAsPairsLocal = ulises.MConnsConvex2RandomPWalkAsPairs240227(mConnsGroupLight)
                    assert all([item[0] < item[1] for item in randomPWalkAsPairsLocal])  # check for readability
                    #
                    # obtain random pseudowalk with (global) names
                    randomPWalkAsPairs = []  # as [name0, name1]
                    for pairLocal in randomPWalkAsPairsLocal:
                        pairGlobalNames = [names[group[pairLocal[item]]] for item in range(2)]
                        randomPWalkAsPairs.append(pairGlobalNames)
                    #
                    # check for readability
                    assert all([item[0] < item[1] for item in randomPWalkAsPairs])
                    assert set(ulises.PWalkAsPairs2PWalkAsList240227(randomPWalkAsPairs)) == set(namesGroup)
                    assert all(ulises.MergeTwoStringsWithTo240227(item[0], item[1]) in dPixelsConns.keys() for item in randomPWalkAsPairs)
                    #
                    # run through pairs
                    cntRnd2 = 0
                    while cntRnd2 < nIter:
                        #
                        # update cntRnd2
                        cntRnd2 = cntRnd2 + 1
                        #
                        # obtain dHsRef and dRSSTsRef
                        dHsRef = ulises.PWalkAsPairs2DHsRef240227(randomPWalkAsPairs, dPixelsConns, nameRefGroup, nPairs, rotation)
                        assert np.isclose(2*len(namesGroup)-1, len(dHsRef.keys()))  # check for readability
                        dRSSTsRef = {key: ulises.H2RSST240227(dHsRef[key]) for key in dHsRef.keys()}
                        if not ulises.DRSSTsRefIsOK240227(dRSSTsRef):
                            continue
                        #
                        # obtain dConnsRecov and quality
                        dConnsRecov = ulises.DConnsRecov240227(namesGroup, nameRefGroup, dPixelsConns, dHsRef, errorC2, nPairs, fraction)
                        quality = ulises.DConns2NOfConns240227(dConnsRecov)
                        #
                        # obtain toSave and update RSST_quality.json and connections_quality.json in pathFldGroup
                        fnsRSSTjson = sorted([item for item in os.listdir(pathFldGroup) if item.startswith('RSST_') and item.endswith('.json')])
                        if len(fnsRSSTjson) == 0:
                            toSave = True
                            cntRnd2 = max([1, cntRnd2 - int(nIter / 2) - 1])
                        else:
                            assert len(fnsRSSTjson) == 1
                            if quality > int(os.path.splitext(fnsRSSTjson[0])[0].split('_')[1]):
                                pathTMP = os.path.join(pathFldGroup, fnsRSSTjson[0])
                                os.remove(pathTMP)
                                os.remove(pathTMP.replace('RSST_', 'connections_'))
                                toSave = True
                                cntRnd2 = max([1, cntRnd2 - int(nIter / 2) - 1])
                            else:
                                toSave = False
                        if toSave:  # update RSST_quality.json and connections_quality.json
                            auxStrH = str(quality).zfill(8)
                            pathRSST = os.path.join(pathFldGroup, 'RSST_{:}.json'.format(auxStrH))
                            fileout = open(pathRSST, 'w')
                            json.dump(dRSSTsRef, fileout, indent=2)
                            fileout.close()
                            pathConns = os.path.join(pathFldGroup, 'connections_{:}.json'.format(auxStrH))
                            fileout = open(pathConns, 'w')
                            json.dump(dConnsRecov, fileout, indent=2)  # IMP*
                            fileout.close()
                            ulises.WriteTimeStampLog(pathFld06View)
                        #
            #
            if ulises.ReadTimeStampLog(pathFld06View) > date0:
                ulises.DeleteSubgroupsFlds240227(pathFldGroup)
                ulises.WriteSubgroupsFlds240227(pathFldGroup, nc4SatCode, nr4SatCode)
            #
            # update maxSizeOfSubgroup
            maxSizeOfSubgroup = max([maxSizeOfSubgroup, ulises.PathFldGroup2MaxSizeOfSubgroup240227(pathFldGroup)])
            #
        #
        print(' s2 = {:5.1f}% of n'.format(maxSizeOfSubgroup / len(names) * 100))
        #
    return None
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
def Output(pathFldMain, par, case='default'):
    #
    # manage case
    if case == 'default':
        xx = '06'
    elif case == 'bypass':
        xx = '16'
    else:
        assert False
    #
    # obtain pathFldScratch and clean
    print('  ... preliminary checkings...')
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    ulises.CleanScratch240227(pathFldMain, par)
    #
    # obtain viewCodes and flds06
    viewCodes = ulises.Par2ViewCodes240227(par)
    flds06 = ulises.ParAndXX2FldsXX240227(par, xx)
    #
    # run for viewCodes x flds06
    for viewCode, fld06 in itertools.product(viewCodes, flds06):
        #
        # obtain paths
        pathFld06View = os.path.join(pathFldScratch, fld06, viewCode)
        if not os.path.exists(pathFld06View):
            continue
        #
        # inform
        print('  ... working on {:} and {:} (generating output)'.format(viewCode, fld06))
        #
        # write output
        ulises.PathFld06View2Out(pathFld06View, par)
        #
    return None
#
