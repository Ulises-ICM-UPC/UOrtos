#
# Mon Feb 17 14:04:42 2025, extract from Ulises by Gonzalo Simarro and Daniel Calvete
#
import copy
import cv2  # type: ignore
import datetime
import itertools
import json
import numpy as np  # type: ignore
import os
from osgeo import gdal, osr  # type: ignore
import osgeo.gdalnumeric as gdn  # type: ignore
import random
import shutil
import zipfile
#
def ApplyH01240227(H01, cs0, rs0):  # lm:2024-04-20; lr:2025-02-08
    RSST01 = H2RSST240227(H01)
    cs1, rs1 = ApplyRSST01240227(RSST01, cs0, rs0)
    return cs1, rs1
def ApplyRSST01240227(RSST01, cs0, rs0):  # lm:2024-03-11; lr:2025-02-08
    cs1 = +RSST01[0] * cs0 + RSST01[1] * rs0 + RSST01[2]
    rs1 = -RSST01[1] * cs0 + RSST01[0] * rs0 + RSST01[3]
    return cs1, rs1
def AreAllMetadataInFldEqual240227(pathFld):  # lm:2025-01-08; lr:2025-02-10
    fnsTifs = [item.name for item in os.scandir(pathFld) if item.is_file() and item.name.endswith('.tif')]
    if len(fnsTifs) == 0:  # IMP*
        return False
    if len(fnsTifs) == 1:
        return True
    for posFnTif, fnTif in enumerate(fnsTifs):
        metadata = PathTif2Metadata240227(os.path.join(pathFld, fnTif))
        if posFnTif == 0:
            metadata0 = copy.deepcopy(metadata)
        elif not AreTwoMetadataEqual240227(metadata, metadata0):
            return False
    return True
def AreTwoMetadataEqual240227(metadata0, metadata1):  # lm:2023-05-31; lr:2025-02-11
    areTwoMetadataEqual = True
    if not set(['width', 'height', 'projection', 'geotransform']) <= set(metadata0.keys()):
        return False
    if not set(['width', 'height', 'projection', 'geotransform']) <= set(metadata1.keys()):
        return False
    if not np.isclose(metadata0['width'], metadata1['width']):
        return False
    if not np.isclose(metadata0['height'], metadata1['height']):
        return False
    if metadata0['projection'] != metadata1['projection']:
        return False
    if not np.allclose(np.asarray(metadata0['geotransform']), np.asarray(metadata1['geotransform'])):
        return False
    return areTwoMetadataEqual
def CR2CRIntegerAroundAndWeights(cs, rs): # lm:2021-09-13; lr:2025-01-08
    csFloor, rsFloor = np.floor(cs).astype(int), np.floor(rs).astype(int)
    csDelta, rsDelta = cs - csFloor, rs - rsFloor
    mZeros = np.zeros((len(cs), 4))
    csIAround, rsIAround, wsAround = mZeros.astype(int), mZeros.astype(int), mZeros
    csIAround[:, 0], rsIAround[:, 0], wsAround[:, 0] = csFloor + 0, rsFloor + 0, (1. - csDelta) * (1. - rsDelta)
    csIAround[:, 1], rsIAround[:, 1], wsAround[:, 1] = csFloor + 1, rsFloor + 0, (0. + csDelta) * (1. - rsDelta)
    csIAround[:, 2], rsIAround[:, 2], wsAround[:, 2] = csFloor + 0, rsFloor + 1, (1. - csDelta) * (0. + rsDelta)
    csIAround[:, 3], rsIAround[:, 3], wsAround[:, 3] = csFloor + 1, rsFloor + 1, (0. + csDelta) * (0. + rsDelta)
    assert np.allclose(wsAround.sum(axis=1), 1)
    possCs0, possCs1 = [np.where(np.abs(csDelta - item) < 1.e-8)[0] for item in [0, 1]]
    possRs0, possRs1 = [np.where(np.abs(rsDelta - item) < 1.e-8)[0] for item in [0, 1]]
    for pos in range(4):  # all pos (corners) are given the same value!
        if len(possCs0) > 0:
            csIAround[possCs0, pos] = csFloor[possCs0]
        if len(possCs1) > 0:
            csIAround[possCs1, pos] = csFloor[possCs1] + 1
        if len(possRs0) > 0:
            rsIAround[possRs0, pos] = rsFloor[possRs0]
        if len(possRs1) > 0:
            rsIAround[possRs1, pos] = rsFloor[possRs1] + 1
    return csIAround, rsIAround, wsAround
def CRXY2Geotranform240227(cs, rs, xs, ys):  # lm:2024-06-07; lr:2025-02-10
    if not (len(cs) == len(rs) == len(xs) == len(ys) and len(cs) >= 3):
        return None
    nOfPoints = len(cs)
    A, b = np.zeros((2 * nOfPoints, 6)), np.zeros(2 * nOfPoints)
    poss0, poss1 = Poss0AndPoss1InFind2DTransform(nOfPoints)
    A[poss0, 0], A[poss0, 1], A[poss0, 2], b[poss0] = np.ones(nOfPoints), cs, rs, xs
    A[poss1, 3], A[poss1, 4], A[poss1, 5], b[poss1] = np.ones(nOfPoints), cs, rs, ys
    try:
        geotransform = np.linalg.lstsq(A, b, rcond=None)[0]
        geotransform = tuple(geotransform)
    except Exception:  # aligned points
        geotransform = None
    return geotransform
def CenteredSubgroupName(dRSSTs, nameRefGroup, namesSubgroup, nc4SatCode, nr4SatCode):  # lm:2024-12-29; lr:2025-02-03
    ms = np.zeros((len(namesSubgroup), len(namesSubgroup)))
    for pos0, pos1 in itertools.product(range(len(namesSubgroup)), range(len(namesSubgroup))):
        name0, name1 = [namesSubgroup[item] for item in [pos0, pos1]]
        RSST0Ref = dRSSTs[MergeTwoStringsWithTo240227(name0, nameRefGroup)]
        RSSTRef1 = dRSSTs[MergeTwoStringsWithTo240227(nameRefGroup, name1)]
        RSST01 = RSST01AndRSST12ToRSST02240227(RSST0Ref, RSSTRef1)
        satCode0 = name0.split('_')[0]
        aux1 = RSSTNcNr2RMSE2240227(RSST01, nc4SatCode[satCode0], nr4SatCode[satCode0])
        ms[pos0, pos1] = np.sqrt(aux1)
    ms = (ms + np.transpose(ms)) / 2
    posCentered = np.argmin(np.sum(ms, axis=0))
    nameCenteredSubgroup = namesSubgroup[posCentered]
    return nameCenteredSubgroup
def CleanScratch240227(pathFldMain, par):  # lm:2025-02-12; lm:2025-02-12
    CleanScratchGivenReference240227(pathFldMain)
    CleanScratchCases240227(pathFldMain, par)
    CleanScratchDates240227(pathFldMain, par)
    return None
def CleanScratchCases240227(pathFldMain, par):  # DANIs cleaner; lm:2024-06-18; lr:2025-02-10
    if not par['clean']:
        return None
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    if not os.path.exists(pathFldScratch):
        return None
    dValidFldsXX = {xx: ParAndXX2FldsXX240227(par, xx) for xx in ['00', '01', '02', '03', '04', '05', '06', '15', '16']}
    viewCodes = Par2ViewCodes240227(par)
    for fldXX in PathFld2Flds240227(pathFldScratch):
        pathFldXX = os.path.join(pathFldScratch, fldXX)
        xx = fldXX.split('_')[0]
        if xx not in dValidFldsXX or fldXX not in dValidFldsXX[xx]:
            shutil.rmtree(pathFldXX)
            continue
        for viewCode in PathFld2Flds240227(pathFldXX):
            if viewCode not in viewCodes:
                shutil.rmtree(os.path.join(pathFldXX, viewCode))
    RmEmptyFldsRecursively(pathFldMain, ['02_pairs'])  # WATCH OUT; 02_pairs is retained
    return None
def CleanScratchDates240227(pathFldMain, par):  # GONZALOs cleaner; lm:2025-01-08; lr:2025-02-10
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    if not os.path.exists(pathFldScratch):
        return None
    viewCodes = ViewCodesInScratch240227(pathFldScratch)
    if len(viewCodes) == 0:
        shutil.rmtree(pathFldScratch)
        return None
    yyyymmdd0, yyyymmdd1 = [par[item].replace('-', '') for item in ['time0', 'time1']]
    for viewCode, fldXX in itertools.product(viewCodes, PathFld2Flds240227(pathFldScratch)):
        pathFldXXView = os.path.join(pathFldScratch, fldXX, viewCode)
        if not os.path.exists(pathFldXXView):
            continue
        if fldXX.split('_')[0] == '00':  # zip files (and dates) rule
            assert fldXX == '00_zip_tif'  # avoidable check for readability
            pathFldXXViewZip = os.path.join(pathFldXXView, 'zip')
            pathFldXXViewTif = os.path.join(pathFldXXView, 'tif')
            if not os.path.exists(pathFldXXViewZip):
                RmTreeIfExists(pathFldXXView)  # including tif and zip
            else:
                namesZips = [os.path.splitext(item)[0] for item in PathFld2Fns240227(pathFldXXViewZip, ext='.zip')]
                if len(namesZips) == 0:
                    shutil.rmtree(pathFldXXView)  # including tif and zip
                else:
                    for fnZip in PathFld2Fns240227(pathFldXXViewZip, ext='.zip'):  # zip files
                        nameZip = os.path.splitext(fnZip)[0]
                        assert nameZip.split('_')[0] != '00'  # WATCH OUT; forbidden satellite code
                        yyyymmddNZip = Name2YYYYMMDD240227(nameZip)
                        if not (nameZip in namesZips and int(yyyymmdd0) <= int(yyyymmddNZip) <= int(yyyymmdd1)):
                            os.remove(os.path.join(pathFldXXViewZip, fnZip))
                            WriteTimeStampLog(pathFldXXView)
                    if os.path.exists(pathFldXXViewTif):
                        for name0 in PathFld2Flds240227(pathFldXXViewTif):  # tif folders
                            if name0.split('_')[0] == '00':  # WATCH OUT; handled by CleanScratchGivenReference240227
                                continue
                            yyyymmddN0 = Name2YYYYMMDD240227(name0)
                            if not (name0 in namesZips and int(yyyymmdd0) <= int(yyyymmddN0) <= int(yyyymmdd1)):
                                shutil.rmtree(os.path.join(pathFldXXViewTif, name0))
                                WriteTimeStampLog(pathFldXXView)
        elif fldXX.split('_')[0] == '01':  # tif files (and dates) rule
            assert fldXX == '01_images'  # avoidable check for readability
            pathFld00ViewTif = os.path.join(pathFldScratch, '00_zip_tif', viewCode, 'tif')
            if not os.path.exists(pathFld00ViewTif):
                shutil.rmtree(pathFldXXView)
            else:
                namesTifs = PathFld2Flds240227(pathFld00ViewTif)
                for fnPng in PathFld2Fns240227(pathFldXXView, ext='.png'):  # IMP*; png
                    namePng = os.path.splitext(fnPng)[0]
                    yyyymmddNPng = Name2YYYYMMDD240227(namePng, yyyymmddG=yyyymmdd1)  # yyyymmddG is only used if namePng.split('_')[0] == '00'
                    if not (namePng in namesTifs and int(yyyymmdd0) <= int(yyyymmddNPng) <= int(yyyymmdd1)):
                        os.remove(os.path.join(pathFldXXView, fnPng))
                        WriteTimeStampLog(pathFldXXView)
        elif fldXX.split('_')[0] == '02':  # png files (and dates) rule
            pathFld01View = os.path.join(pathFldScratch, '01_images', viewCode)
            if not os.path.exists(pathFld01View):
                shutil.rmtree(pathFldXXView)
            else:
                namesPngs = [os.path.splitext(item)[0] for item in PathFld2Fns240227(pathFld01View, ext='.png')]  # IMP*
                for name0 in PathFld2Flds240227(pathFldXXView):  # folders
                    pathFldXXViewName0 = os.path.join(pathFldXXView, name0)
                    yyyymmddN0 = Name2YYYYMMDD240227(name0, yyyymmddG=yyyymmdd1)  # yyyymmddG is only used if name0.split('_')[0] == '00'
                    if not (name0 in namesPngs and int(yyyymmdd0) <= int(yyyymmddN0) <= int(yyyymmdd1)):
                        shutil.rmtree(pathFldXXViewName0)
                        WriteTimeStampLog(pathFldXXView)
                    else:  # look inside pathFldXXViewName0
                        for fnNpz in PathFld2Fns240227(pathFldXXViewName0, ext='.npz'):
                            nameNpz = os.path.splitext(fnNpz)[0]
                            yyyymmddNNpz = Name2YYYYMMDD240227(nameNpz, yyyymmddG=yyyymmdd1)  # yyyymmddG is only used if nameNpz.split('_')[0] == '00'
                            if not (nameNpz in namesPngs and int(yyyymmdd0) <= int(yyyymmddNNpz) <= int(yyyymmdd1)):
                                os.remove(os.path.join(pathFldXXViewName0, fnNpz))
                                WriteTimeStampLog(pathFldXXView)
        elif fldXX.split('_')[0] in ['03', '04', '05', '15']:  # previous npz (and dates) rule
            yy = {'03': '02', '04': '03', '05': '04', '15': '02'}[fldXX.split('_')[0]]  # yy < xx = fldXX.split('_')[0]
            fldYY = FldXX2FldYY240227(fldXX, yy, options={})  # options = {} since yy < xx = fldXX.split('_')[0]
            for name0 in PathFld2Flds240227(pathFldXXView):
                pathFldXXViewName0 = os.path.join(pathFldXXView, name0)
                pathFldYYViewName0 = os.path.join(pathFldScratch, fldYY, viewCode, name0)
                yyyymmddN0 = Name2YYYYMMDD240227(name0, yyyymmddG=yyyymmdd1)  # yyyymmddG is only used if name0.split('_')[0] == '00'
                if not (os.path.exists(pathFldYYViewName0) and int(yyyymmdd0) <= int(yyyymmddN0) <= int(yyyymmdd1)):
                    shutil.rmtree(pathFldXXViewName0)
                    WriteTimeStampLog(pathFldXXView)
                else:  # look inside pathFldXXViewName0
                    for fnNpz in PathFld2Fns240227(pathFldXXViewName0, ext='.npz'):
                        pathXXViewName0Fn = os.path.join(pathFldXXViewName0, fnNpz)
                        pathYYViewName0Fn = os.path.join(pathFldYYViewName0, fnNpz)
                        nameNpz = os.path.splitext(fnNpz)[0]
                        yyyymmddNNpz = Name2YYYYMMDD240227(nameNpz, yyyymmddG=yyyymmdd1)  # yyyymmddG is only used if nameNpz.split('_')[0] == '00'
                        if not (os.path.exists(pathYYViewName0Fn) and int(yyyymmdd0) <= int(yyyymmddNNpz) <= int(yyyymmdd1)):
                            os.remove(pathXXViewName0Fn)
                            WriteTimeStampLog(pathFldXXView)
        elif fldXX.split('_')[0] in ['06', '16']:  # via logs
            yy = {'06': '05', '16': '15'}[fldXX.split('_')[0]]  # yy < xx = fldXX.split('_')[0]
            fldYY = FldXX2FldYY240227(fldXX, yy, options={})  # options = {} since yy < xx = fldXX.split('_')[0]
            pathFldYYView = os.path.join(pathFldScratch, fldYY, viewCode)
            if not os.path.exists(pathFldYYView):
                shutil.rmtree(pathFldXXView)
            elif int(ReadTimeStampLog(pathFldYYView)) > int(ReadTimeStampLog(pathFldXXView)):
                shutil.rmtree(pathFldXXView)
        else:
            assert False
    RmEmptyFldsRecursively(pathFldMain, ['02_pairs'])
    return None
def CleanScratchGivenReference240227(pathFldMain):  # lm:2025-01-08; lr:2025-02-10
    pathFldData = os.path.join(pathFldMain, 'data')
    pathFldScratch = os.path.join(pathFldMain, 'scratch')
    if not os.path.exists(pathFldScratch):
        return None
    for viewCode in ViewCodesInScratch240227(pathFldScratch):
        if os.path.exists(os.path.join(pathFldData, '{:}_reference.tif'.format(viewCode))):
            continue  # IMP*; no work to do ---to change the given reference tif: delete the old + run the code + add the new---
        pathFld = os.path.join(pathFldScratch, '00_zip_tif', viewCode, 'tif', '00_reference')  # IMP*
        if os.path.exists(pathFld):
            shutil.rmtree(pathFld)
            WriteTimeStampLog(os.path.join(pathFldScratch, '00_zip_tif', viewCode))
        pathFile = os.path.join(pathFldScratch, '01_images', viewCode, '00_reference.png')  # IMP*
        if os.path.exists(pathFile):
            os.remove(pathFile)
            WriteTimeStampLog(os.path.join(pathFldScratch, '01_images', viewCode))
        for fldXX in [item for item in PathFld2Flds240227(pathFldScratch) if item[:2] in ['02', '03', '04', '05', '15']]:
            pathFld = os.path.join(pathFldScratch, fldXX, viewCode, '00_reference')
            if os.path.exists(pathFld):
                shutil.rmtree(pathFld)
                WriteTimeStampLog(os.path.join(pathFldScratch, fldXX, viewCode))
            pathFld = os.path.join(pathFldScratch, fldXX, viewCode)
            for root, _, fns in os.walk(pathFld):  # 
                for fn in fns:
                    if '00_reference' in fn:
                        os.remove(os.path.join(root, fn))
                        WriteTimeStampLog(os.path.join(pathFldScratch, fldXX, viewCode))
        for fld06 in [item for item in PathFld2Flds240227(pathFldScratch) if item[:2] in ['06', '16']]:
            pathFld06View = os.path.join(pathFldScratch, fld06, viewCode)
            names = PathFld06View2Names240227(pathFld06View)
            if '00_reference' in names:
                RmTreeIfExists(pathFld06View)
    RmEmptyFldsRecursively(pathFldMain, ['02_pairs'])  # WATCH OUT; 02_pairs is retained; the work has been done
    return None
def CompleteADictionary(theDictionary, keys, defaultValues): # lm:2021-09-10; lr:2025-02-06
    if set(keys) <= set(theDictionary.keys()):
        pass
    else:
        if isinstance(defaultValues, list):
            assert len(keys) == len(defaultValues)
            for posKey, key in enumerate(keys):
                if key not in theDictionary.keys():  # only assigns if there is no key
                    theDictionary[key] = defaultValues[posKey]
        else:  # defaultValues is a single value
            for key in keys:
                if key not in theDictionary.keys():  # only assigns if there is no key
                    theDictionary[key] = defaultValues
    return theDictionary
def DConns2NOfConns240227(dConns):  # lm:2024-05-10; lr:2025-01-08
    mConns = DConns2NamesAndMConns240227(dConns)[1]
    nOfConns = MConns2NOfConns240227(mConns)
    return nOfConns
def DConns2NameRef240227(dConns):  # lm:2024-06-06; lr:2025-02-11
    names, mConns = DConns2NamesAndMConns240227(dConns)
    if len(names) == 1:
        nameRef = names[0]
    else:
        posRef = np.argmax(np.sum(mConns, axis=0))  # WATCH OUT; np.argmax has repetitivity; cheap
        nameRef = names[posRef]
    return nameRef
def DConns2NamesAndMConns240227(dConns):  # lm:2024-05-10; lr:2025-02-11
    names = sorted(dConns)  # sorted(dConns.keys())
    mConns = np.zeros((len(names), len(names)), dtype=int)
    for pos0, name0 in enumerate(names):
        for name1 in dConns[name0]:
            pos1 = names.index(name1)
            mConns[pos0, pos1] = 1
    assert np.allclose(mConns, np.transpose(mConns))
    assert np.isclose(np.sum(np.abs(np.diagonal(mConns))), 0)
    return names, mConns
def DConnsRecov240227(names, nameRef, dPixelsConns, dHsRef, errorC, nPairs, fraction):  # lm:2024-05-07; lr:2025-02-11
    assert nameRef in names
    dConnsRecov = {key: [] for key in names}
    for name0, name1 in itertools.product(names, names):
        str01 = MergeTwoStringsWithTo240227(name0, name1)
        if str01 not in dPixelsConns.keys():
            continue
        assert name0 < name1  # check for readability (MConnsAndDPixelsConns240227)
        H0ToRef = dHsRef[MergeTwoStringsWithTo240227(name0, nameRef)]
        HRefTo1 = dHsRef[MergeTwoStringsWithTo240227(nameRef, name1)]
        H01 = H01AndH12ToH02240227(H0ToRef, HRefTo1)
        nc0, nr0, cs0, rs0 = [dPixelsConns[str01][item] for item in ['nc0', 'nr0', 'cs0', 'rs0']]
        _, _, cs1, rs1 = [dPixelsConns[str01][item] for item in ['nc1', 'nr1', 'cs1', 'rs1']]
        cs1R, rs1R = ApplyH01240227(H01, cs0, rs0)
        possGood = np.where(np.sqrt((cs1R - cs1) ** 2 + (rs1R - rs1) ** 2) <= errorC)[0]
        if len(possGood) < nPairs:
            continue
        if not (np.max(cs0[possGood]) - np.min(cs0[possGood]) >= fraction * nc0 and np.max(rs0[possGood]) - np.min(rs0[possGood]) >= fraction * nr0):
            continue
        listTMP = dConnsRecov[name0] + [name1]
        dConnsRecov[name0] = sorted(listTMP)
        listTMP = dConnsRecov[name1] + [name0]
        dConnsRecov[name1] = sorted(listTMP)
    return dConnsRecov
def DRSSTsRef2NameRef240227(dRSSTsRef):  # lm:2024-06-07; lr:2025-02-08
    DRSSTsRefCheck240227(dRSSTsRef)
    keys = list(dRSSTsRef)
    nameRef = SplitTwoStringsWithTo240227(keys[0])[0]
    return nameRef
def DRSSTsRef2Names240227(dRSSTsRef):  # lm:2025-02-08; lr:2025-02-08
    DRSSTsRefCheck240227(dRSSTsRef)
    keys = list(dRSSTsRef)
    names0 = sorted(set([SplitTwoStringsWithTo240227(item)[0] for item in keys]))
    names1 = sorted(set([SplitTwoStringsWithTo240227(item)[1] for item in keys]))
    assert names0 == names1
    names = names0  # IMP*; sorted
    return names
def DRSSTsRefCheck240227(dRSSTsRef):  # lm:2025-01-02; lr:2025-02-08
    keys = list(dRSSTsRef.keys())
    strA, strB = SplitTwoStringsWithTo240227(keys[0])  # first key is nameRef_to_nameRef
    assert strA == strB
    assert all(strA in key for key in keys)  # nameRef is in all keys
    return None
def DRSSTsRefIsOK240227(dRSSTsRef):  # lm:2024-03-13; lr:2025-02-08
    DRSSTsRefCheck240227(dRSSTsRef)
    isOK = True
    for str01 in dRSSTsRef.keys():
        str10 = ReverseTwoStringsWithTo240227(str01)
        if str01 > str10:
            continue
        H01 = RSST2H240227(dRSSTsRef[str01])
        H10 = RSST2H240227(dRSSTsRef[str10])
        if not np.allclose(np.dot(H01, H10), np.eye(3)):
            isOK = False
            break
    return isOK
def DRSSTsRefNew240227(dRSSTsRef, namesNew, nameRefNew):  # lm:2024-05-23; lr:2025-02-08
    DRSSTsRefCheck240227(dRSSTsRef)
    assert nameRefNew in namesNew
    nameRefOld = DRSSTsRef2NameRef240227(dRSSTsRef)
    namesOld = DRSSTsRef2Names240227(dRSSTsRef)
    assert set(namesNew) <= set(namesOld)
    strRefNewRefNew = MergeTwoStringsWithTo240227(nameRefNew, nameRefNew)
    dRSSTsRefNew = {strRefNewRefNew: H2RSST240227(np.eye(3))}
    RSSTRefOldRefNew = dRSSTsRef[MergeTwoStringsWithTo240227(nameRefOld, nameRefNew)]
    RSSTRefNewRefOld = dRSSTsRef[MergeTwoStringsWithTo240227(nameRefNew, nameRefOld)]
    for nameNew in [item for item in namesNew if item != nameRefNew]:
        strNewRefOld = MergeTwoStringsWithTo240227(nameNew, nameRefOld)
        strRefOldNew = MergeTwoStringsWithTo240227(nameRefOld, nameNew)
        strNewRefNew = MergeTwoStringsWithTo240227(nameNew, nameRefNew)
        strRefNewNew = MergeTwoStringsWithTo240227(nameRefNew, nameNew)
        dRSSTsRefNew[strNewRefNew] = RSST01AndRSST12ToRSST02240227(dRSSTsRef[strNewRefOld], RSSTRefOldRefNew)
        dRSSTsRefNew[strRefNewNew] = RSST01AndRSST12ToRSST02240227(RSSTRefNewRefOld, dRSSTsRef[strRefOldNew])
        assert np.allclose(np.dot(RSST2H240227(dRSSTsRefNew[strNewRefNew]), RSST2H240227(dRSSTsRefNew[strRefNewNew])), np.eye(3))
    return dRSSTsRefNew
def Datetime2Date17(theDatetime):  # lm:2025-02-01; lr:2025-02-12
    year, month, day = str(theDatetime.year).zfill(4), str(theDatetime.month).zfill(2), str(theDatetime.day).zfill(2)
    hour, minute = str(theDatetime.hour).zfill(2), str(theDatetime.minute).zfill(2)
    second, millisecond = str(theDatetime.second).zfill(2), str(int(theDatetime.microsecond / 1.e+3)).zfill(3)
    date17 = '{:}{:}{:}{:}{:}{:}{:}'.format(year, month, day, hour, minute, second, millisecond)
    assert len(date17) == 17
    return date17
def DeleteSubgroupsFlds240227(pathFldGroup):  # lm:2024-11-27; lr:2025-02-11
    for pathFldSubgroup in [item.path for item in os.scandir(pathFldGroup) if item.is_dir() and item.name.startswith('subgroup_')]:
        shutil.rmtree(pathFldSubgroup)
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
    return None
def DoTwoPWalksAsPairsTouch240227(pWalkAsPairs0, pWalkAsPairs1):  # share a node, not an segment; lm:2024-03-03; lr:2025-02-08  
    pWalkAsList0, pWalkAsList1 = [PWalkAsPairs2PWalkAsList240227(item) for item in [pWalkAsPairs0, pWalkAsPairs1]]
    doTwoPWalksAsPairsTouch = bool(set(pWalkAsList0) & set(pWalkAsList1))
    return doTwoPWalksAsPairsTouch
def EnoughImprovingFBoxes240227(nc, nr, nOfCBands, nOfRBands, cs, rs, nPairs, nOfFBoxes0):  # lm:2024-04-16; lr:2025-02-08
    nOfFBoxes = Pixels2NOfFBoxes240227(nc, nr, nOfCBands, nOfRBands, cs, rs)
    improvingFBoxes = nOfFBoxes >= nPairs and nOfFBoxes > nOfFBoxes0
    return improvingFBoxes
def FMPre240227(pathFldImgs, fnsImgs, methodFM, scale, nOfFeatures):  # lm:2025-01-03; lr:2025-02-12
    assert isinstance(scale, int)
    optionsTMP = {'nOfFeatures': nOfFeatures}
    ncsOri, nrsOri, ncsSca, nrsSca, kpssSca, dessSca = [[] for _ in range(6)]
    for fnImg in fnsImgs:
        imgOri = cv2.imread(os.path.join(pathFldImgs, fnImg))
        nrOri, ncOri = imgOri.shape[0:2]
        ncsOri.append(ncOri), nrsOri.append(nrOri)
        imgSca = cv2.resize(imgOri, (int(np.round(scale*ncOri)), int(np.round(scale*nrOri))))
        assert np.isclose(imgSca.shape[1], scale*ncOri) and np.isclose(imgSca.shape[0], scale*nrOri)  # check for readability
        if methodFM == 'ORB':
            ncSca, nrSca, kpsSca, desSca, ctrl = ORBKeypoints(imgSca, options=optionsTMP)
        elif methodFM == 'SIFT':
            ncSca, nrSca, kpsSca, desSca, ctrl = SIFTKeypoints(imgSca, options=optionsTMP)
        else:
            assert False
        if ctrl:
            ncsSca.append(ncSca), nrsSca.append(nrSca), kpssSca.append(kpsSca), dessSca.append(desSca)
            assert np.isclose(ncSca, scale*ncOri) and np.isclose(nrSca, scale*nrOri)  # check for readability
        else:
            ncsSca.append(None), nrsSca.append(None), kpssSca.append(None), dessSca.append(None)  # IMP*
    return ncsOri, nrsOri, ncsSca, nrsSca, kpssSca, dessSca
def FindForcedRSST01240227(dt01, cs0, rs0, cs1, rs1, nOfIter=10):  # lm:2024-04-17; lr:2025-02-08
    assert len(cs0) == len(rs0) == len(cs1) == len(rs1) >= 2
    nOfPoints = len(cs0)
    RSST01 = FindRSST01240227(cs0, rs0, cs1, rs1)
    if RSST01 is None:
        return None
    if False:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        cs1R, rs1R = ApplyRSST01240227(RSST01, cs0, rs0)
        print('{:8.5f} (dt01 = {:8.5f}, {:} points)'.format(np.abs(np.sqrt(np.mean((cs1R-cs1) ** 2 + (rs1R-rs1) ** 2))), np.sqrt(RSST01[0] ** 2 + RSST01[1] ** 2), len(cs0)))
    ng0, dc0, dr0 = np.angle(RSST01[0] + 1j * RSST01[1]), RSST01[2], RSST01[3]
    hasConverged = False
    for _ in range(nOfIter):
        aux1, cs1p = cs0 * dt01 * np.cos(ng0) + rs0 * dt01 * np.sin(ng0), cs1 - dc0
        aux2, rs1p = cs0 * dt01 * np.sin(ng0) - rs0 * dt01 * np.cos(ng0), rs1 - dr0
        e1dng, e1ddx, e1ddy, e10 = aux2, -np.ones(nOfPoints), np.zeros(nOfPoints), cs1p - aux1
        e2dng, e2ddx, e2ddy, e20 = aux1, np.zeros(nOfPoints), -np.ones(nOfPoints), rs1p + aux2
        e3dng, e3ddx, e3ddy, e30 = aux1 * cs1p - aux2 * rs1p, - aux2, - aux1, aux2 * cs1p + aux1 * rs1p
        A, b = np.zeros((3, 3)), np.zeros(3)
        A[0, 0], A[0, 1], A[0, 2], b[0] = np.sum(e1dng), np.sum(e1ddx), np.sum(e1ddy), -np.sum(e10)
        A[1, 0], A[1, 1], A[1, 2], b[1] = np.sum(e2dng), np.sum(e2ddx), np.sum(e2ddy), -np.sum(e20)
        A[2, 0], A[2, 1], A[2, 2], b[2] = np.sum(e3dng), np.sum(e3ddx), np.sum(e3ddy), -np.sum(e30)
        try:
            sol = np.linalg.solve(A, b)
        except Exception:
            hasConverged = False
            break
        ng0 = ng0 + sol[0]
        dc0 = dc0 + sol[1]
        dr0 = dr0 + sol[2]
        RSST01 = [dt01*np.cos(ng0), dt01*np.sin(ng0), dc0, dr0]
        if False:
            cs1R, rs1R = ApplyRSST01240227(RSST01, cs0, rs0)
            print('{:8.5f} (dt01 = {:8.5f}, {:} points)'.format(np.abs(np.sqrt(np.mean((cs1R-cs1) ** 2 + (rs1R-rs1) ** 2))), np.sqrt(RSST01[0] ** 2 + RSST01[1] ** 2), len(cs0)))
        de1 = np.sum((cs1-cs0*dt01*np.cos(ng0)-rs0*dt01*np.sin(ng0)-dc0)*(+cs0*dt01*np.sin(ng0)-rs0*dt01*np.cos(ng0)))
        de1 = np.sum((rs1+cs0*dt01*np.sin(ng0)-rs0*dt01*np.cos(ng0)-dr0)*(+cs0*dt01*np.cos(ng0)+rs0*dt01*np.sin(ng0))) + de1
        de2 = np.sum(cs1-cs0*dt01*np.cos(ng0)-rs0*dt01*np.sin(ng0)-dc0)
        de3 = np.sum(rs1+cs0*dt01*np.sin(ng0)-rs0*dt01*np.cos(ng0)-dr0)
        if max([np.abs(item) for item in [de1, de2, de3]]) < 1.e-8:
            hasConverged = True
            break
    if hasConverged:
        RSST01 = [dt01*np.cos(ng0), dt01*np.sin(ng0), dc0, dr0]  # list
    else:
        RSST01 = None
    return RSST01
def FindGoodPossRANSAC240227(dt01, cs0, rs0, cs1, rs1, parRANSAC, nc0, nr0, fraction, rotation, nOfCBands, nOfRBands, nOfIter=10):  # lm:2025-02-08; lr:2025-02-08
    assert np.max(cs0) - np.min(cs0) >= fraction * nc0 and np.max(rs0) - np.min(rs0) >= fraction * nr0  # check for readability
    assert len(cs0) == len(rs0) == len(cs1) == len(rs1) >= 2  # (comes from "if len(ers) < 4:  # IMP*; BLOCKED")
    NOfRANSAC = NForRANSAC(parRANSAC['pOutlier'], parRANSAC['pDesired'], parRANSAC['nOfPoints'])
    possGood, counter = [], 0
    while counter < NOfRANSAC:
        poss2 = random.sample(range(0, len(cs0)), 2)
        cs0H, rs0H, cs1H, rs1H = [item[poss2] for item in [cs0, rs0, cs1, rs1]]
        if rotation == 'free':
            RSST01 = FindForcedRSST01240227(dt01, cs0H, rs0H, cs1H, rs1H, nOfIter=nOfIter)
        elif rotation == 'null':
            assert np.isclose(dt01, 1)  # IMP*; no dilatation
            RSST01 = FindRSST01Translation240227(cs0, rs0, cs1, rs1)
        else:
            assert False
        if RSST01 is None:
            continue
        cs1R, rs1R = ApplyRSST01240227(RSST01, cs0, rs0)
        ers1 = np.sqrt((cs1R - cs1) ** 2 + (rs1R - rs1) ** 2)
        possGoodH = np.where(ers1 < parRANSAC['errorC'])[0]  # in "1"
        pOutlier = min(0.5, 1 - len(possGoodH) / len(cs0))
        NOfRANSAC = NForRANSAC(pOutlier, parRANSAC['pDesired'], parRANSAC['nOfPoints'])
        counter = counter + 1
        if len(possGoodH) == 0:
            continue
        cs0H, rs0H, cs1H, rs1H, ers1H = [item[possGoodH] for item in [cs0, rs0, cs1, rs1, ers1]]
        if not (np.max(cs0H) - np.min(cs0H) >= fraction * nc0 and np.max(rs0H) - np.min(rs0H) >= fraction * nr0):
            continue
        possSelectInGoodH = SelectPossInGrid240227(nc0, nr0, nOfCBands, nOfRBands, cs0H, rs0H, ers1H, fraction)
        if len(possSelectInGoodH) == 0:
            continue
        cs0H, rs0H, cs1H, rs1H, ers1H = [item[possSelectInGoodH] for item in [cs0H, rs0H, cs1H, rs1H, ers1H]]  # IMP*
        assert np.max(cs0H) - np.min(cs0H) >= fraction * nc0 and np.max(rs0H) - np.min(rs0H) >= fraction * nr0  # check for readability
        if len(possSelectInGoodH) > len(possGood):  # IMP*
            possGood = [possGoodH[item] for item in possSelectInGoodH]  # IMP*
    return possGood
def FindRSST01240227(cs0, rs0, cs1, rs1):  # lm:2024-03-11; lr:2025-02-08
    assert len(cs0) == len(rs0) == len(cs1) == len(rs1) >= 2
    nOfPoints = len(cs0)
    A, b = np.zeros((2 * nOfPoints, 4)), np.zeros(2 * nOfPoints)
    poss0, poss1 = Poss0AndPoss1InFind2DTransform(nOfPoints)
    A[poss0, 0], A[poss0, 1], A[poss0, 2], b[poss0] = cs0,  rs0, np.ones(nOfPoints), cs1
    A[poss1, 0], A[poss1, 1], A[poss1, 3], b[poss1] = rs0, -cs0, np.ones(nOfPoints), rs1
    try:
        RSST01 = list(np.linalg.lstsq(A, b, rcond=None)[0])
    except Exception:  # aligned points
        RSST01 = None
    return RSST01
def FindRSST01Translation240227(cs0, rs0, cs1, rs1):  # lm:2024-04-26; lr:2025-02-08
    assert len(cs0) == len(rs0) == len(cs1) == len(rs1) >= 1
    RSST01 = [1, 0, np.mean(cs1 - cs0), np.mean(rs1 - rs0)]
    return RSST01
def FldXX240227(xx, options={}):  # lm:2024-06-02; lr:2025-02-10
    if xx == '00':
        fldXX = '00_zip_tif'
    elif xx == '01':
        fldXX = '01_images'
    elif xx == '02':  # nFM
        nFMStr = Integer2Str240227(options['nFM'])
        fldXX = '{:}_pairs_{:}'.format(xx, nFMStr)
        assert len(fldXX.split('_')) == 3  # avoidable check for readability
    elif xx in ['03', '04']:  # nFM, errorC1, fraction, nPairs, nBoxes, rotation
        nFMStr = Integer2Str240227(options['nFM'])
        errorC1Str = Float2Str240227(options['errorC1'])
        fractionStr = Float2Str240227(options['fraction'])
        nPairsStr = Integer2Str240227(options['nPairs'])
        nBoxesStr = Integer2Str240227(options['nBoxes'])
        rotation = options['rotation']
        aux = '_'.join([nFMStr, errorC1Str, fractionStr, nPairsStr, nBoxesStr, rotation])
        fldXX = '{:}_pairs_{:}'.format(xx, aux)
        assert len(fldXX.split('_')) == 8  # avoidable check for readability
    elif xx in ['05', '15']:  # nFM, errorC1, fraction, nPairs, nBoxes, rotation, errorC2
        nFMStr = Integer2Str240227(options['nFM'])
        errorC1Str = Float2Str240227(options['errorC1'])
        fractionStr = Float2Str240227(options['fraction'])
        nPairsStr = Integer2Str240227(options['nPairs'])
        nBoxesStr = Integer2Str240227(options['nBoxes'])
        rotation = options['rotation']
        errorC2Str = Float2Str240227(options['errorC2'])
        aux = '_'.join([nFMStr, errorC1Str, fractionStr, nPairsStr, nBoxesStr, rotation, errorC2Str])
        fldXX = '{:}_pairs_{:}'.format(xx, aux)
        assert len(fldXX.split('_')) == 9  # avoidable check for readability
    elif xx in ['06', '16']:  # nFM, errorC1, fraction, nPairs, nBoxes, rotation, errorC2, conDegree
        nFMStr = Integer2Str240227(options['nFM'])
        errorC1Str = Float2Str240227(options['errorC1'])
        fractionStr = Float2Str240227(options['fraction'])
        nPairsStr = Integer2Str240227(options['nPairs'])
        nBoxesStr = Integer2Str240227(options['nBoxes'])
        rotation = options['rotation']
        errorC2Str = Float2Str240227(options['errorC2'])
        conDegreeStr = Integer2Str240227(options['conDegree'])
        aux = '_'.join([nFMStr, errorC1Str, fractionStr, nPairsStr, nBoxesStr, rotation, errorC2Str, conDegreeStr])
        fldXX = '{:}_affines_{:}'.format(xx, aux)
        assert len(fldXX.split('_')) == 10  # avoidable check for readability
    return fldXX
def FldXX2DInfo240227(fldXX):  # lm:2024-06-02; lr:2025-02-10
    fldXXSplit = fldXX.split('_')
    dInfo = {'fldXX': fldXX, 'xx': fldXXSplit[0]}
    if fldXXSplit[0] == '00':
        assert fldXX == '00_zip_tif'  # check for readability
    elif fldXXSplit[0] == '01':
        assert fldXX == '01_images'  # check for readability
    elif fldXXSplit[0] == '02':
        assert len(fldXXSplit) == 3 and fldXXSplit[1] == 'pairs'  # check for readability
        dInfo['nFM'] = Str2Integer240227(fldXXSplit[2])
    elif fldXXSplit[0] in ['03', '04']:
        assert len(fldXXSplit) == 8 and fldXXSplit[1] == 'pairs'  # check for readability
        dInfo['nFM'] = Str2Integer240227(fldXXSplit[2])
        dInfo['errorC1'] = Str2Float240227(fldXXSplit[3])
        dInfo['fraction'] = Str2Float240227(fldXXSplit[4])
        dInfo['nPairs'] = Str2Integer240227(fldXXSplit[5])
        dInfo['nBoxes'] = Str2Integer240227(fldXXSplit[6])
        dInfo['rotation'] = fldXXSplit[7]
    elif fldXXSplit[0] in ['05', '15']:
        assert len(fldXXSplit) == 9 and fldXXSplit[1] == 'pairs'  # check for readability
        dInfo['nFM'] = Str2Integer240227(fldXXSplit[2])
        dInfo['errorC1'] = Str2Float240227(fldXXSplit[3])
        dInfo['fraction'] = Str2Float240227(fldXXSplit[4])
        dInfo['nPairs'] = Str2Integer240227(fldXXSplit[5])
        dInfo['nBoxes'] = Str2Integer240227(fldXXSplit[6])
        dInfo['rotation'] = fldXXSplit[7]
        dInfo['errorC2'] = Str2Float240227(fldXXSplit[8])
    elif fldXXSplit[0] in ['06', '16']:
        assert len(fldXXSplit) == 10 and fldXXSplit[1] == 'affines'  # check for readability
        dInfo['nFM'] = Str2Integer240227(fldXXSplit[2])
        dInfo['errorC1'] = Str2Float240227(fldXXSplit[3])
        dInfo['fraction'] = Str2Float240227(fldXXSplit[4])
        dInfo['nPairs'] = Str2Integer240227(fldXXSplit[5])
        dInfo['nBoxes'] = Str2Integer240227(fldXXSplit[6])
        dInfo['rotation'] = fldXXSplit[7]
        dInfo['errorC2'] = Str2Float240227(fldXXSplit[8])
        dInfo['conDegree'] = Str2Integer240227(fldXXSplit[9])
    assert fldXX == FldXX240227(fldXXSplit[0], options=dInfo)  # check for readability
    return dInfo
def FldXX2FldYY240227(fldXX, yy, options={}):  # lm:2024-12-29; lr:2025-02-10
    dInfoXX = FldXX2DInfo240227(fldXX)
    keysCommon = set(options.keys()) & set(dInfoXX.keys())
    assert all(options[key] == dInfoXX[key] for key in keysCommon)  # IMP*
    options = options | dInfoXX
    fldYY = FldXX240227(yy, options=options)
    return fldYY
def Float2Str240227(flt, length=4):  # lm:2024-03-07; lr:2025-02-10
    fltx100 = int(np.round(flt * 100))
    assert np.isclose(fltx100 / 100, flt)  # IMP*: 2 non-null decimals at most
    string = str(fltx100).zfill(length)  # IMP*
    return string
def GeotranformCR2XY240227(geotransform, cs, rs):  # lm:2024-06-07; lr:2025-02-10
    xs = geotransform[0] + cs * geotransform[1] + rs * geotransform[2]
    ys = geotransform[3] + cs * geotransform[4] + rs * geotransform[5]
    return xs, ys
def GeotranformXY2CR240227(geotransform, xs, ys):  # lm:2024-06-07; lr:2025-02-10
    m = np.zeros((2, 2))
    m[0, 0], m[0, 1] = geotransform[1], geotransform[2]
    m[1, 0], m[1, 1] = geotransform[4], geotransform[5]
    mInv = np.linalg.inv(m)
    xsd = xs - geotransform[0]
    ysd = ys - geotransform[3]
    cs = mInv[0, 0] * xsd + mInv[0, 1] * ysd
    rs = mInv[1, 0] * xsd + mInv[1, 1] * ysd
    xsR, ysR = GeotranformCR2XY240227(geotransform, cs, rs)
    assert np.allclose(xs, xsR) and np.allclose(ys, ysR)
    return cs, rs
def H01AndH12ToH02240227(H01, H12):  # lm:2024-02-29; lr:2025-02-08
    H02 = np.dot(H12, H01)
    return H02
def H2RSST240227(H):  # lm:2025-02-08; lr:2025-02-08
    assert all(np.isclose(item, 0) for item in [H[0, 0] - H[1, 1], H[0, 1] + H[1, 0], H[2, 0], H[2, 1], H[2, 2] - 1])
    RSST = [H[0, 0], H[0, 1], H[0, 2], H[1, 2]]  # IMP*; list
    assert np.allclose(RSST2H240227(RSST), H)
    return RSST
def Img2ImgColorClipped240227(img):  # lm:2023-06-07; lr:2025-02-10
    imgHSV = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    hs, ss, vs = cv2.split(imgHSV)
    v02, v98 = np.percentile(vs, 2), np.percentile(vs, 98)
    vs = np.clip(((vs - v02) / (v98 - v02) * 255), 1, 254).astype(np.uint8)
    imgHSV = cv2.merge((hs, ss, vs))
    img = cv2.cvtColor(imgHSV, cv2.COLOR_HSV2BGR)
    return img
def ImgA2ImgAAsInB240227(imgA, RSSTImgBToImgA, ncB=None, nrB=None, margin=0):  # lm:2025-02-10; lr:2025-02-10
    nrA, ncA = imgA.shape[0:2]
    if ncB is None or nrB is None:
        nrB, ncB = imgA.shape[0:2]
    csB, rsB = np.meshgrid(np.arange(ncB), np.arange(nrB))
    csB, rsB = csB.ravel(), rsB.ravel()
    if margin > 0:
        color = [0, 0, 0]
        imgA[:margin, :, :], imgA[-margin:, :, :], imgA[:, :margin, :], imgA[:, -margin:, :] = color, color, color, color
    imgAAsInB = np.zeros((nrB, ncB, 3))
    csA, rsA = ApplyRSST01240227(RSSTImgBToImgA, csB, rsB)
    csAIA, rsAIA, wsAA = CR2CRIntegerAroundAndWeights(csA, rsA)
    csAIA, rsAIA = np.clip(csAIA, 0, ncA-1), np.clip(rsAIA, 0, nrA-1)
    wsAAAux = np.asarray([wsAA, wsAA, wsAA])
    for corner in range(4):
        auxs = imgA[rsAIA[:, corner], csAIA[:, corner], :] * np.transpose(wsAAAux[:, :, corner])
        imgAAsInB[rsB, csB, :] = imgAAsInB[rsB, csB, :] + auxs
    imgAAsInB = imgAAsInB.astype(np.uint8)
    return imgAAsInB
def ImgA2ImgAAsInB250210(imgA, RSSTImgAToImgB):  # lm:2025-02-10; lr:2025-02-10
    nr, nc = imgA.shape[0:2]
    imgAAsInB = cv2.warpAffine(imgA, RSST2H240227(RSSTImgAToImgB)[:2, :], (nc, nr))
    return imgAAsInB
def ImprovePairsThroughCorrelation240227(cs0, rs0, cs1, rs1, img0, img1, sameSize, semiSize=15, scale0=10):  # lm:2024-06-11; lr:2025-02-11
    csDom, rsDom, csPat, rsPat = cs0, rs0, cs1, rs1
    RSSTDomPat = FindForcedRSST01240227(1, csDom, rsDom, csPat, rsPat, nOfIter=10)  # angle always free, OK
    if RSSTDomPat is None:  # not converged
        cs0, rs0, cs1, rs1 = [np.asarray([]) for _ in range(4)]
        return cs0, rs0, cs1, rs1
    RSSTPatDom = InvertRSST240227(RSSTDomPat)
    imgDom, imgPat = img0, img1
    nrDom, ncDom = imgDom.shape[0:2]
    if sameSize:
        imgPot = ImgA2ImgAAsInB250210(imgPat, RSSTPatDom)  # works if images have the same image
    else:
        imgPot = ImgA2ImgAAsInB240227(imgPat, RSSTDomPat, ncB=ncDom, nrB=nrDom, margin=0)  # WATCH OUT; can be expensive
    assert imgPot.shape == imgDom.shape
    assert np.isclose(scale0, int(scale0))
    scale = int(scale0)  # IMP*; WATCH OUT <= 10 (integer)
    nrS, ncS, semiSizeS = int(imgDom.shape[0] * scale), int(imgDom.shape[1] * scale), int(semiSize * scale)
    imgPotS, imgDomS = [cv2.resize(item, (ncS, nrS)) for item in [imgPot, imgDom]]
    imgPotS, imgDomS = [np.mean(item, axis=2) for item in [imgPotS, imgDomS]]
    csDomInt, rsDomInt, csPot, rsPot = [[] for _ in range(4)]
    for posH in range(len(csDom)):
        cDomHIntS, rDomHIntS = int(csDom[posH]) * scale, int(rsDom[posH]) * scale
        condCDom = (cDomHIntS > semiSizeS + 3*scale) and (cDomHIntS < ncS - 1 - semiSizeS - 3*scale)
        condRDom = (rDomHIntS > semiSizeS + 3*scale) and (rDomHIntS < nrS - 1 - semiSizeS - 3*scale)
        if not (condCDom and condRDom):
            continue
        imgDomSCropH = imgDomS[rDomHIntS-semiSizeS:rDomHIntS+semiSizeS, :][:, cDomHIntS-semiSizeS:cDomHIntS+semiSizeS]
        sigmaDomH = np.sqrt(TwoMatrices2Corr240227(imgDomSCropH, imgDomSCropH))
        rhoOpt, cPotHIntSOpt, rPotHIntSOpt = 0, None, None
        for iteration in range(1):  # WATCH OUT; 1 or 2  (Montecarlo)
            if iteration == 0:
                cPotHIntSBest, rPotHIntSBest = 1 * cDomHIntS, 1 * rDomHIntS
            else:
                auxC, auxR = [int(4 * (np.random.random() - 0.5) * scale) for _ in range(2)]
                cPotHIntSBest, rPotHIntSBest = 1 * cDomHIntS + auxC, 1 * rDomHIntS + auxR
            rhoBest = 0
            while True:  # in each loop we analyze a 3x3 square around the current best
                rhoBestSq, posCBestSq, posRBestSq = 0, None, None  # Sq refer to 3 x 3 square
                for posC, posR in itertools.product(range(-1, 2), range(-1, 2)):
                    cSH, rSH = cPotHIntSBest + posC, rPotHIntSBest + posR
                    imgPotSCropH = imgPotS[rSH-semiSizeS:rSH+semiSizeS, :][:, cSH-semiSizeS:cSH+semiSizeS]
                    sigmaPotH = np.sqrt(TwoMatrices2Corr240227(imgPotSCropH, imgPotSCropH))
                    sigmaPotDomH = TwoMatrices2Corr240227(imgPotSCropH, imgDomSCropH)
                    rhoH = sigmaPotDomH / (sigmaPotH * sigmaDomH)
                    if rhoH > rhoBestSq:
                        rhoBestSq, posCBestSq, posRBestSq = [copy.deepcopy(item) for item in [rhoH, posC, posR]]
                if any(item is None for item in [posCBestSq, posRBestSq]):
                    isOK = False
                    break  # gradient
                if rhoBestSq > rhoBest:
                    rhoBest, cPotHIntSBest, rPotHIntSBest = 1 * rhoBestSq, cPotHIntSBest + posCBestSq, rPotHIntSBest + posRBestSq
                distance = np.sqrt((cPotHIntSBest - cDomHIntS) ** 2 + (rPotHIntSBest - rDomHIntS) ** 2) / scale
                if posCBestSq == 0 and posRBestSq == 0 or distance > 3:
                    isOK = distance < 3 and rhoBest > 0.5  # WATCH OUT; 3 could be larger
                    break  # while True, i.e., the kind of gradient
            if isOK and rhoBest > rhoOpt:
                rhoOpt, cPotHIntSOpt, rPotHIntSOpt = [copy.deepcopy(item) for item in [rhoBest, cPotHIntSBest, rPotHIntSBest]]
        if not (rhoOpt > 0 and cPotHIntSOpt is not None and rPotHIntSOpt is not None):
            continue  # go to next GCP
        csDomInt.append(cDomHIntS/scale), rsDomInt.append(rDomHIntS/scale), csPot.append(cPotHIntSOpt/scale), rsPot.append(rPotHIntSOpt/scale)
    csDom, rsDom, csPot, rsPot = [np.asarray(item) for item in [csDomInt, rsDomInt, csPot, rsPot]]
    if len(csDom) == 0:
        cs0, rs0, cs1, rs1 = [np.asarray([]) for _ in range(4)]
        return cs0, rs0, cs1, rs1
    csPat, rsPat = ApplyRSST01240227(RSSTDomPat, csPot, rsPot)  # IMP*
    assert len(csDom) == len(rsDom) == len(csPat) == len(rsPat)
    cs0, rs0, cs1, rs1 = csDom, rsDom, csPat, rsPat
    return cs0, rs0, cs1, rs1
def Integer2Str240227(integer, length=5):  # lm:2024-03-14; lr:2025-02-10
    assert np.isclose(integer, int(np.round(integer)))
    string = str(int(np.round(integer))).zfill(length)
    return string
def InvertRSST240227(RSST01):  # lm:2024-03-18; lr:2025-02-08
    H01 = RSST2H240227(RSST01)
    H10 = np.linalg.inv(H01)
    RSST10 = H2RSST240227(H10)
    return RSST10
def LDesiredNOfConns240227(sizeOfGroup, nOfConns):  # lm:2024-03-04; lr:2025-01-08
    aux0, aux1 = min(sizeOfGroup, nOfConns), nOfConns
    lDesiredNOfConns = sorted(np.linspace(aux0, aux1, 25, dtype=int).tolist())
    assert lDesiredNOfConns[-1] == nOfConns  # check for readability
    return lDesiredNOfConns
def LPWalksAsPairs2LPWalksAsLists240227(lPWalksAsPairs):  # lm:2024-03-03; lr:2025-02-08
    lPWalksAsLists = [PWalkAsPairs2PWalkAsList240227(item) for item in lPWalksAsPairs]
    return lPWalksAsLists
def MConns2Groups240227(mConns, conDegree):  # lm:2024-05-09; lr:2025-02-11
    groups, possRemaining = [], list(range(mConns.shape[0]))
    while len(possRemaining) > 0:  # obtains a group and updates possRemaining in each loop
        mConnsRemaining = mConns[:, possRemaining][possRemaining, :]
        pos0 = possRemaining[np.argmax(np.sum(mConnsRemaining, axis=0))]  # best connected in remaining; np.argmax has repetitivity
        if sum(mConns[pos0, possRemaining]) < conDegree:
            group = [pos0]  # group finished
        else:
            group = [pos0] + [item for item in possRemaining if mConns[pos0, item] > 0]  # initialize group (could be improved)
            while True:  # increase group
                lenOld = len(group)
                for posH in possRemaining:
                    if posH not in group and np.sum(mConns[posH, group]) >= conDegree:
                        group.append(posH)
                assert sorted(set(group)) == sorted(group)  # avoidable check for readability
                if len(group) == lenOld:
                    break
            while True:  # decrease group (deleting bad positions derived from line "group = [pos0] + [item... ")
                lenOld = len(group)
                for posH in reversed(group):  # to be able to remove items
                    if np.sum(mConns[posH, group]) < conDegree:
                        group.remove(posH)  # remove only removes one, if repeated (not the case)
                        break
                if len(group) == lenOld or len(group) == 1:
                    break
        group = sorted(group)
        groups.append(group)
        possRemaining = [item for item in possRemaining if item not in group]
    return groups
def MConns2LSortedPairs240227(mConns):  # lm:2024-03-04; lr:2025-02-11
    lSortedPairs = []
    for pos0, pos1 in itertools.product(range(mConns.shape[0]), range(mConns.shape[0])):
        if pos0 < pos1 and mConns[pos0, pos1] > 0:  # IMP*: pos0 < pos1
            lSortedPairs.append([pos0, pos1])
    assert len(lSortedPairs) == MConns2NOfConns240227(mConns)  # check for readability
    return lSortedPairs
def MConns2NOfConns240227(mConns):  # lm:2024-05-10; lr:2025-02-11
    nOfConns = np.sum(np.reshape(mConns, -1)) / 2
    assert np.isclose(nOfConns, int(np.round(nOfConns)))
    nOfConns = int(np.round(nOfConns))
    return nOfConns
def MConnsAndDPixelsConns240227(names, pathFld05View):  # lm:2024-04-25; lr:2025-02-11
    assert names == sorted(names)
    mConns, dPixelsConns = np.zeros((len(names), len(names)), dtype=int), {}
    for root, _, fns in os.walk(pathFld05View):
        fnsNpz05 = [item for item in fns if item.endswith('.npz')]
        if len(fnsNpz05) == 0:
            continue
        name0 = os.path.split(root)[1]
        posImg0 = names.index(name0)
        for fnNpz05 in fnsNpz05:
            name1 = os.path.splitext(fnNpz05)[0]
            posImg1 = names.index(name1)
            assert posImg0 < posImg1  # check for readability
            mConns[posImg0, posImg1] = 1
            mConns[posImg1, posImg0] = 1
            dInfo = PathNpzOfPairs2DInfo240227(os.path.join(root, fnNpz05), options={})
            str01 = MergeTwoStringsWithTo240227(name0, name1)
            dPixelsConns[str01] = {key: dInfo[key] for key in ['nc0', 'nr0', 'nc1', 'nr1', 'cs0', 'rs0', 'cs1', 'rs1']}
    assert np.allclose(mConns, np.transpose(mConns))
    assert np.isclose(np.sum(np.abs(np.diagonal(mConns))), 0)
    return mConns, dPixelsConns
def MConnsConvex2RandomPWalkAsPairs240227(mConns):  # lm:2024-03-03; lr:2025-02-08
    assert len(MConns2Groups240227(mConns, 1)) == 1  # convex (degree = 1)
    sizeMConns = mConns.shape[0]
    lRemainingPairs = MConns2LSortedPairs240227(mConns)
    lPWalksAsPairs = []
    while len(lRemainingPairs) > 0:  # this is pretty beautiful
        pair = random.sample(lRemainingPairs, 1)[0]
        lPWalksAsPairs.append([pair])
        lPWalksAsPairs = MergeLPWalksAsPairs240227(lPWalksAsPairs)
        lRemainingPairs = UpdateLRemainingPairs240227(lRemainingPairs, lPWalksAsPairs)
    assert len(lPWalksAsPairs) == 1
    randomPWalkAsPairs = lPWalksAsPairs[0]
    assert len(randomPWalkAsPairs) == sizeMConns - 1
    assert set(PWalkAsPairs2PWalkAsList240227(randomPWalkAsPairs)) == set(range(sizeMConns))
    assert all([mConns[item[0], item[1]] == 1 for item in randomPWalkAsPairs])
    return randomPWalkAsPairs
def MConnsLight240227(mConns, conDegree):  # lm:2024-06-06; lr:2025-02-11
    assert len(MConns2Groups240227(mConns, conDegree)) == 1  # convex
    nOfConns = MConns2NOfConns240227(mConns)
    lSortedPairs = MConns2LSortedPairs240227(mConns)  # local and so that [pos0, pos1] with pos0 < pos1
    lDesiredNOfConns = LDesiredNOfConns240227(mConns.shape[0], nOfConns)
    mConnsLightFound = False
    for desiredNOfConns in lDesiredNOfConns:
        mConnsLight = np.zeros_like(mConns, dtype=int)
        pairs = random.sample(lSortedPairs, desiredNOfConns)
        for pair in pairs:
            mConnsLight[pair[0], pair[1]] = 1
            mConnsLight[pair[1], pair[0]] = 1
        if len(MConns2Groups240227(mConnsLight, conDegree)) == 1:  # one group
            mConnsLightFound = True
            break
    if not mConnsLightFound:
        mConnsLight = None
    return mConnsLight
def MFMConns240227(pathFld02View, names):  # lm:2024-06-13; lr:2025-02-06
    assert names == PathFld2Flds240227(pathFld02View)
    mFMConns = np.zeros((len(names), len(names)), dtype=int)
    for pos0 in range(len(names)):
        pathFld02Name0 = os.path.join(pathFld02View, names[pos0])
        assert os.path.exists(pathFld02Name0)  # avoidable; for readability
        for fn in PathFld2Fns240227(pathFld02Name0, ext='.npz'):
            pos1 = names.index(os.path.splitext(fn)[0])
            if np.isclose(pos0, pos1):  # this should not happen
                continue
            mFMConns[pos0, pos1], mFMConns[pos1, pos0] = 1, 1
    return mFMConns
def MFMConnsAddPairs240227(pathFld02View, methodFM, names, ncsOri, nrsOri, kpssSca, dessSca, scaFM, limFM):  # lm:2024-06-19; lr:2025-02-06
    assert names == PathFld2Flds240227(pathFld02View) and limFM > 0
    mFMConns = MFMConns240227(pathFld02View, names)
    vSum = np.sum(mFMConns, axis=0)
    workDone = False
    success, trials = np.zeros(len(names)), np.ones(len(names))
    mTried = 1 * mFMConns
    for c in range(mTried.shape[0]):  # artifact not to try pos0 == pos1
        mTried[c, c] = 1
    for _ in range(min([len(names), limFM]) * len(names)):  # IMP*
        pos0, pos1 = sorted(random.sample(list(range(len(names))), 2))
        if mTried[pos0, pos1] > 0 or kpssSca[pos0] is None or kpssSca[pos1] is None:
            continue
        if vSum[pos0] >= limFM and vSum[pos1] >= limFM:  # IMP*; "and" -> refreshing is allowed thru Deletion
            continue
        if np.clip(min([success[pos0] / trials[pos0], success[pos1] / trials[pos1]]), 0.1, 0.9) > np.random.random():  # IMP*; favour less successful
            continue
        trials[pos0], trials[pos1] = trials[pos0] + 1, trials[pos1] + 1
        mTried[pos0, pos1], mTried[pos1, pos0] = 1, 1  # not to try it again
        kps0Sca, des0Sca, kps1Sca, des1Sca = kpssSca[pos0], dessSca[pos0], kpssSca[pos1], dessSca[pos1]
        nc0Ori, nr0Ori, nc1Ori, nr1Ori = ncsOri[pos0], nrsOri[pos0], ncsOri[pos1], nrsOri[pos1]
        pathNpz02 = os.path.join(pathFld02View, names[pos0], '{:}.npz'.format(names[pos1]))
        doesWrite = WritePathNpz02240227(pathNpz02, methodFM, nc0Ori, nr0Ori, nc1Ori, nr1Ori, kps0Sca, des0Sca, kps1Sca, des1Sca, scaFM)
        if doesWrite:
            success[pos0], success[pos1] = success[pos0] + 1, success[pos1] + 1
            mFMConns[pos0, pos1], mFMConns[pos1, pos0] = 1, 1
            vSum = np.sum(mFMConns, axis=0)
            workDone = True
    return workDone
def MFMConnsDeletePairs240227(pathFld02View, names, limFM):  # lm:2024-06-13; lr:2025-02-06
    assert names == PathFld2Flds240227(pathFld02View) and limFM > 0
    workDone = False
    while True:
        workDoneH = False
        mFMConns = MFMConns240227(pathFld02View, names)
        vSum = np.sum(mFMConns, axis=0)
        poss0 = list(np.where(vSum > limFM)[0])
        if len(poss0) == 0:  # we are done (workDoneH and workDone are False)
            break
        for pos0 in random.sample(poss0, len(poss0)):
            if vSum[pos0] <= limFM:
                continue
            poss1 = list(np.where(mFMConns[:, pos0] > 0)[0])
            poss1 = [poss1[item] for item in np.argsort(vSum[poss1])[::-1]]
            for pos1 in poss1:
                assert not np.isclose(pos0, pos1)  # for readability
                if vSum[pos1] <= limFM:
                    continue
                if pos0 > pos1:
                    pos0, pos1 = pos1, pos0
                assert pos0 < pos1  # for readability
                pathNpz02 = os.path.join(pathFld02View, names[pos0], '{:}.npz'.format(names[pos1]))
                RmFileIfExists(pathNpz02)  # it should always exist
                WriteTimeStampLog(pathFld02View)
                mFMConns[pos0, pos1], mFMConns[pos1, pos0] = 0, 0
                vSum = np.sum(mFMConns, axis=0)
                workDoneH, workDone = True, True
        if not workDoneH:
            break
    return workDone
def MergeLPWalksAsPairs240227(lPWalksAsPairs):  # lm:2024-03-03; lr:2025-02-08
    while True:
        lenOld = len(lPWalksAsPairs)
        lPWalksAsPairs = MergeLPWalksAsPairsOneStep240227(lPWalksAsPairs)
        if len(lPWalksAsPairs) == lenOld:
            break
    return lPWalksAsPairs
def MergeLPWalksAsPairsOneStep240227(lPWalksAsPairs):  # lm:2024-02-08; lr:2025-02-08
    if len(lPWalksAsPairs) <= 1:
        return lPWalksAsPairs
    lPWalksAsPairsOld = lPWalksAsPairs
    lPWalksAsPairsNew = [lPWalksAsPairsOld[0]]  # initialize with the first pWalk
    for pWalkAsPairsOld in lPWalksAsPairsOld[1:]:  # starts 1; tries to merge it some of the paths in new, or add
        isMerged = False
        for posNew, pWalkAsPairsNew in enumerate(lPWalksAsPairsNew):
            if DoTwoPWalksAsPairsTouch240227(pWalkAsPairsOld, pWalkAsPairsNew):
                lPWalksAsPairsNew[posNew] = pWalkAsPairsNew + pWalkAsPairsOld
                isMerged = True
                break  # pWalkAsPairsOld has been merged to a new
        if not isMerged:
            lPWalksAsPairsNew.append(pWalkAsPairsOld)  # IMP*
    lPWalksAsPairs = lPWalksAsPairsNew
    return lPWalksAsPairs
def MergeTwoStringsWithTo240227(str0, str1):  # lm:2024-03-07; lr:2025-02-12
    assert all('_to' not in item for item in [str0, str1])  # WATCH OUT; forbidden names
    assert all('to_' not in item for item in [str0, str1])  # WATCH OUT; forbidden names
    str01 = '{:}_to_{:}'.format(str0, str1)
    return str01
def NForRANSAC(pOutlier, pDesired, nOfPoints): # lm:2022-11-09; lr:2025-02-07
    num = np.log(min(1-1.e-12, max(1.e-12, 1. - pDesired)))
    den = np.log(min(1-1.e-12, max(1.e-12, 1. - (1. - pOutlier) ** nOfPoints)))
    N = int(num / den) + 1
    return N
def NOfBands240227(nBoxes, nc, nr):  # lm:2024-04-16; lr:2025-02-12
    nOfCBands = int(np.ceil(np.sqrt(nBoxes * nc / nr)))
    nOfRBands = int(np.ceil(np.sqrt(nBoxes * nr / nc)))
    return nOfCBands, nOfRBands
def Name2YYYYMMDD240227(name, yyyymmddG=None):  # lm:2025-02-11; lr:2025-02-11
    if name.split('_')[0] in ['00']:  # IMP*; WATCH OUT; the given reference
        assert yyyymmddG is not None
        yyyymmdd = yyyymmddG
    elif name.split('_')[0] in ['L5', 'L7', 'L8', 'L9']:
        yyyymmdd = name.split('_')[3][0:8]  # IMP*; 3
    elif name.split('_')[0] in ['S2']:
        yyyymmdd = name.split('_')[1][0:8]  # IMP*; 1
    else:
        assert False
    return yyyymmdd
def NamesAndMConns2DConns240227(names, mConns):  # lm:2024-12-29; lr:2025-02-11
    dConns = {}
    for pos0, name0 in enumerate(names):  # run through all (symmetry)
        names1 = sorted([names[pos] for pos in range(len(names)) if mConns[pos0, pos] > 0])
        dConns[name0] = names1
    for name0 in dConns:
        names1 = dConns[name0]
        assert name0 not in names1
        for name1 in names1:
            assert name0 in dConns[name1]
    return dConns
def ORBKeypoints(img, options={}): # lm:2021-09-22; lr:2025-02-03
    keys, defaultValues = ['mask', 'nOfFeatures'], [None, 5000]
    options = CompleteADictionary(options, keys, defaultValues)
    try:
        img = PathImgOrImg2Img(img)
        nr, nc = img.shape[0:2]
        orb = cv2.ORB_create(nfeatures=options['nOfFeatures'], scoreType=cv2.ORB_FAST_SCORE)  # HARRIS_SCORE or ORB_FAST_SCORE
        if options['mask'] is not None:
            kps, des = orb.detectAndCompute(img, mask=options['mask'])
        else:
            kps, des = orb.detectAndCompute(img)  # it was (img, None) 2025-02-03
        assert len(kps) == len(des) > 0
        ctrl = True
    except Exception:
        nc, nr, kps, des, ctrl = None, None, None, None, False
    return nc, nr, kps, des, ctrl
def ORBMatches(kps1, des1, kps2, des2, options={}): # lm:2024-02-23; lr:2025-02-03
    keys, defaultValues = ['erMaximum', 'nOfStd', 'resolution1', 'resolution2'], [20., 2., 1., 1.]
    options = CompleteADictionary(options, keys, defaultValues)
    cs1, rs1 = [np.asarray([item.pt[pos] for item in kps1]) for pos in [0, 1]]
    cs2, rs2 = [np.asarray([item.pt[pos] for item in kps2]) for pos in [0, 1]]
    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
    matches12 = sorted(bf.match(des1, des2), key = lambda x:x.distance)
    matches21 = sorted(bf.match(des2, des1), key = lambda x:x.distance)
    poss1 = [match.queryIdx for match in matches12] + [match.trainIdx for match in matches21]
    poss2 = [match.trainIdx for match in matches12] + [match.queryIdx for match in matches21]
    cs1, rs1, cs2, rs2 = cs1[poss1], rs1[poss1], cs2[poss2], rs2[poss2]
    ers = np.asarray([match.distance for match in matches12] + [match.distance for match in matches21])
    cs1, rs1, cs2, rs2, ers = np.unique(np.asarray([cs1, rs1, cs2, rs2, ers]), axis=1)  # IMP* interesting
    if len(cs1) > 0:
        xs1, xs2 = cs1 * options['resolution1'], cs2 * options['resolution2']
        ys1, ys2 = rs1 * options['resolution1'], rs2 * options['resolution2']
        ds = np.sqrt((xs1 - xs2) ** 2 + (ys1 - ys2) ** 2)
        possGood = np.where((ers < options['erMaximum']) & (ds < np.mean(ds) + options['nOfStd'] * np.std(ds) + 1.e-8))[0]
        cs1, rs1, cs2, rs2, ers = [item[possGood] for item in [cs1, rs1, cs2, rs2, ers]]
    return cs1, rs1, cs2, rs2, ers
def PWalkAsPairs2DHsRef240227(pWalkAsPairs, dPixelsConns, strRef, nPairs, rotation):  # lm:2024-04-18; lr:2025-02-11
    strs = PWalkAsPairs2PWalkAsList240227(pWalkAsPairs)
    assert strRef in strs
    strsAlr, dHsRef = [strRef], {MergeTwoStringsWithTo240227(strRef, strRef): np.eye(3)}
    while len(strsAlr) < len(strs):
        for pair in pWalkAsPairs:
            if pair[0] in strsAlr and pair[1] not in strsAlr:
                strAlr, strNew = pair[0], pair[1]
            elif pair[1] in strsAlr and pair[0] not in strsAlr:
                strAlr, strNew = pair[1], pair[0]
            else:
                continue
            csNew, rsNew, csAlr, rsAlr = Str0Str12Cs0Rs0Cs1Rs1240227(strNew, strAlr, dPixelsConns)
            while True:
                possN = random.sample(range(len(csNew)), nPairs)  # WATCH OUT; nPairs
                csNewH, rsNewH, csAlrH, rsAlrH = [item[possN] for item in [csNew, rsNew, csAlr, rsAlr]]
                if rotation == 'free':
                    RSSTNewToAlr = FindForcedRSST01240227(1, csNewH, rsNewH, csAlrH, rsAlrH, nOfIter=10)  # dt = 1
                elif rotation == 'null':
                    RSSTNewToAlr = FindRSST01Translation240227(csNewH, rsNewH, csAlrH, rsAlrH)
                else:
                    assert False
                if RSSTNewToAlr is not None:  # WATCH OUT; to avoid aligned
                    break
            HNewToAlr = RSST2H240227(RSSTNewToAlr)
            strNewRef = MergeTwoStringsWithTo240227(strNew, strRef)
            if strAlr == strRef:  # straightforward
                dHsRef[strNewRef] = HNewToAlr
            else:
                strAlrRef = MergeTwoStringsWithTo240227(strAlr, strRef)
                HNewToRef = H01AndH12ToH02240227(HNewToAlr, dHsRef[strAlrRef])
                dHsRef[strNewRef] = HNewToRef
            strRefNew = MergeTwoStringsWithTo240227(strRef, strNew)
            dHsRef[strRefNew] = np.linalg.inv(dHsRef[strNewRef])
            strsAlr.append(strNew)
    return dHsRef
def PWalkAsPairs2PWalkAsList240227(pWalkAsPairs):  # lm:2024-02-08; lr:2025-02-08
    pWalkAsList = list(set(np.reshape(np.asarray(pWalkAsPairs), -1)))  # not necessarily sorted
    return pWalkAsList
def Par2ViewCodes240227(par):  # lm:2025-01-03; lr:2025-02-12
    viewCodes = ['{:}_{:}'.format(resCode, domainCode) for resCode, domainCode in itertools.product(par['resCodes'], par['domainCode2Limits'])]
    return viewCodes
def ParAndXX2FldsXX240227(par, xx):  # lm:2024-06-18; lr:2025-02-10
    fldsXX = []
    if xx in ['00', '01']:
        fldsXX.append(FldXX240227(xx))
    else:
        if xx == '02':
            auxs = [par[item] for item in ['nsFM']]
            keys = ['nFM']
        elif xx in ['03', '04']:
            auxs = [par[item] for item in ['nsFM', 'errorsC1', 'fractions', 'nsPairs', 'nsBoxes', 'rotations']]
            keys = ['nFM', 'errorC1', 'fraction', 'nPairs', 'nBoxes', 'rotation']
        elif xx in ['05', '15']:
            auxs = [par[item] for item in ['nsFM', 'errorsC1', 'fractions', 'nsPairs', 'nsBoxes', 'rotations', 'errorsC2']]
            keys = ['nFM', 'errorC1', 'fraction', 'nPairs', 'nBoxes', 'rotation', 'errorC2']
        elif xx in ['06', '16']:
            auxs = [par[item] for item in ['nsFM', 'errorsC1', 'fractions', 'nsPairs', 'nsBoxes', 'rotations', 'errorsC2', 'conDegrees']]
            keys = ['nFM', 'errorC1', 'fraction', 'nPairs', 'nBoxes', 'rotation', 'errorC2', 'conDegree']
        else:
            assert False
        for values in itertools.product(*auxs):
            fldsXX.append(FldXX240227(xx, options=dict(zip(keys, values))))
    return fldsXX
def PathFld06View2LSubgroupDInfo240227(pathFld06View):  # lm:2024-05-23; lr:2025-02-10
    maxSizeOfSubgroup = 0
    for root, _, _ in sorted(os.walk(pathFld06View)):
        if not os.path.split(root)[1].startswith('subgroup_'):
            continue
        dInfoH = PathFld06ViewSubgroup2DInfo240227(root)
        if dInfoH['size'] > maxSizeOfSubgroup:
            dInfo = copy.deepcopy(dInfoH)
            maxSizeOfSubgroup = 1 * dInfoH['size']
    if maxSizeOfSubgroup <= 1:  # IMP*
        dInfo = None
    return dInfo
def PathFld06View2Names240227(pathFld06View):  # lm:2025-02-12; lr:2025-02-12
    names = set()
    for fld in PathFld2Flds240227(pathFld06View):
        pathTMP = os.path.join(pathFld06View, fld, 'names.txt')
        namesH = set(np.atleast_1d(np.loadtxt(pathTMP, usecols=0, dtype=str)))
        assert not (names & namesH)
        names.update(namesH)
    names = sorted(names)
    return names
def PathFld06View2Out(pathFld06View, par):  # lm:2024-06-09; lr:2025-02-10
    pathFld06, viewCode = os.path.split(pathFld06View)
    pathFldScratch, fld06 = os.path.split(pathFld06)
    pathFldMain, scratch = os.path.split(pathFldScratch)
    assert scratch == 'scratch'
    conDegree = FldXX2DInfo240227(fld06)['conDegree']
    fld06 = '_'.join(fld06.split('_')[2:])  # to remove '06_affines'
    pathFldOutView = os.path.join(pathFldMain, 'output', fld06, viewCode)
    if ReadTimeStampLog(pathFld06View) > ReadTimeStampLog(pathFldOutView):
        RmTreeIfExists(pathFldOutView)
    if os.path.exists(pathFldOutView):
        return None
    os.makedirs(pathFldOutView, exist_ok=True)
    dInfo06LSubgroup = PathFld06View2LSubgroupDInfo240227(pathFld06View)
    if dInfo06LSubgroup is None:
        return None
    names, nameRef, dRSSTsRef = [dInfo06LSubgroup[item] for item in ['names', 'nameRef', 'dRSSTsRef']]
    if len(names) < conDegree + 1:  # IMP*; scape
        return None
    pathFld00ViewTif = os.path.join(pathFldScratch, '00_zip_tif', viewCode, 'tif')
    pathFld00NameRef = os.path.join(pathFld00ViewTif, nameRef)
    metadataRef = PathTif2Metadata240227(os.path.join(pathFld00NameRef, os.listdir(pathFld00NameRef)[0]))  # WATCH OUT; any tif if OK, already checked
    projectionRef, geotransformRef = [metadataRef[item] for item in ['projection', 'geotransform']]
    nc, nr = metadataRef['width'], metadataRef['height']
    csRef, rsRef = np.asarray([0, nc-1, nc-1, 0]), np.asarray([0, 0, nr-1, nr-1])
    xsRef, ysRef = GeotranformCR2XY240227(geotransformRef, csRef, rsRef)
    pathNameRef = os.path.join(pathFldOutView, '{:}.ref'.format(nameRef))
    os.makedirs(os.path.dirname(pathNameRef), exist_ok=True)
    fileout = open(pathNameRef, 'w')
    for pos in range(len(csRef)):
        fileout.write('{:12.3f} {:12.3f} {:12.3f} {:12.3f}\n'.format(csRef[pos], rsRef[pos], xsRef[pos], ysRef[pos]))
    fileout.close()
    WriteTimeStampLog(pathFldOutView)
    for name in names:
        RSSTNameRefToName = dRSSTsRef[MergeTwoStringsWithTo240227(nameRef, name)]
        pathFldBands = os.path.join(pathFld00ViewTif, name)
        metadata = PathTif2Metadata240227(os.path.join(pathFldBands, os.listdir(pathFldBands)[0]))  # WATCH OUT; any tif if OK, already checked
        projection = metadata['projection']
        csAux, rsAux = ApplyRSST01240227(RSSTNameRefToName, csRef, rsRef)
        xsAux, ysAux = XsYsA2XsYsB240227(projectionRef, projection, xsRef, ysRef)
        geotransformM = CRXY2Geotranform240227(csAux, rsAux, xsAux, ysAux)
        for fn in os.listdir(pathFldBands):
            pathTif = os.path.join(pathFldBands, fn)
            pathTifM = os.path.join(pathFldOutView, 'tif', name, fn)
            os.makedirs(os.path.dirname(pathTifM), exist_ok=True)
            PathTif02PathTif1WithGeotransform1240227(pathTif, geotransformM, pathTifM)
            WriteTimeStampLog(pathFldOutView)
            metadataMH = PathTif2Metadata240227(pathTifM)
            assert np.allclose(np.asarray(metadataMH['geotransform']), np.asarray(geotransformM))
            assert all(metadataMH[item] == metadata[item] for item in ['projection', 'height', 'width'])
    csRef, rsRef = np.meshgrid(np.arange(metadataRef['width']).astype(int), np.arange(metadataRef['height']).astype(int))
    csRef, rsRef = [np.reshape(item, -1) for item in [csRef, rsRef]]
    xsRef, ysRef = GeotranformCR2XY240227(geotransformRef, csRef, rsRef)
    for name in names:
        bandsCodes = par['satCode2RGBBandsCodes'][name.split('_')[0]]
        pathFldBands = os.path.join(pathFldOutView, 'tif', name)
        metadata = PathTif2Metadata240227(os.path.join(pathFldBands, os.listdir(pathFldBands)[0]))  # WATCH OUT; any tif if OK
        projection = metadata['projection']
        nc, nr = metadataRef['width'], metadataRef['height']
        cs, rs = [1 * item for item in [csRef, rsRef]]
        xs, ys = XsYsA2XsYsB240227(projectionRef, projection, xsRef, ysRef)
        imgGeo = PathFldBands2ImgGeo240227(pathFldBands, bandsCodes, nc, nr, cs, rs, xs, ys, margin=2)
        if imgGeo is not None:
            pathPng = os.path.join(pathFldOutView, 'geopng', '{:}.png'.format(name))
            os.makedirs(os.path.dirname(pathPng), exist_ok=True)
            cv2.imwrite(pathPng, imgGeo)
            WriteTimeStampLog(pathFldOutView)
        else:
            continue
    return None
def PathFld06ViewSubgroup2DInfo240227(pathFld06Subgroup):  # lm:2024-06-07; lr:2025-02-10
    assert os.path.exists(os.path.join(pathFld06Subgroup, 'names.txt'))
    assert len([item for item in os.scandir(pathFld06Subgroup) if item.is_file() and item.name.endswith('.ref')]) == 1
    assert len(os.listdir(pathFld06Subgroup)) in [2, 3]
    dInfo = {}
    pathFld06Group, subgroupCode = os.path.split(pathFld06Subgroup)
    assert subgroupCode.split('_')[0] == 'subgroup'
    dInfo['pathFld06Group'], dInfo['subgroupCode'] = pathFld06Group, subgroupCode
    pathFld06View, groupCode = os.path.split(pathFld06Group)
    assert groupCode.split('_')[0] == 'group'
    dInfo['pathFld06View'], dInfo['groupCode'] = pathFld06View, groupCode
    pathFld06, viewCode = os.path.split(pathFld06View)
    dInfo['pathFld06'], dInfo['viewCode'] = pathFld06, viewCode
    pathFldScratch, fld06 = os.path.split(pathFld06)
    dInfo['pathFldScratch'], dInfo['fld06'] = pathFldScratch, fld06
    pathFldMain, scratch = os.path.split(pathFldScratch)
    assert scratch == 'scratch'
    dInfo['pathFldMain'] = pathFldMain
    pathTMP = os.path.join(pathFld06Subgroup, 'names.txt')
    names = np.atleast_1d(np.loadtxt(pathTMP, usecols=0, dtype=str)).tolist()
    dInfo['names'], dInfo['size'] = names, len(names)
    nameRef = [os.path.splitext(item.name)[0] for item in os.scandir(pathFld06Subgroup) if item.is_file() and item.name.endswith('.ref')][0]
    assert nameRef in names  # check for readability
    dInfo['nameRef'] = nameRef
    pathTMP = os.path.join(pathFld06Subgroup, 'RSST.json')
    if os.path.exists(pathTMP):
        with open(pathTMP, 'r') as f:
            dRSSTsRef = json.load(f)
        assert names == DRSSTsRef2Names240227(dRSSTsRef) and nameRef == DRSSTsRef2NameRef240227(dRSSTsRef)
    else:
        dRSSTsRef = None
        assert len(os.listdir(pathFld06Subgroup)) == 2
    dInfo['dRSSTsRef'] = dRSSTsRef
    return dInfo
def PathFld2Flds240227(pathFld):  # non-recursive; lm:2024-11-27; lr:2025-02-10
    if os.path.exists(pathFld):
        flds = sorted([item.name for item in os.scandir(pathFld) if item.is_dir()])
    else:
        flds = []
    return flds
def PathFld2Fns240227(pathFld, ext=None):  # non-recursive; lm:2024-11-27; lr:2025-02-10
    if os.path.exists(pathFld):
        fns = sorted([item.name for item in os.scandir(pathFld) if item.is_file()])
        if ext is not None:
            fns = [item for item in fns if item.endswith(ext)]
    else:
        fns = []
    return fns
def PathFldBands2ImgGeo240227(pathFldBands, bandsCodes, nc, nr, cs, rs, xs, ys, margin=0):  # lm:2024-06-07; lr:2025-02-10
    AreAllMetadataInFldEqual240227(pathFldBands)
    metadata = PathTif2Metadata240227(os.path.join(pathFldBands, os.listdir(pathFldBands)[0]))  # WATCH OUT; all are tif
    geotransform = metadata['geotransform']
    imgRaw = PathFldBands2ImgRGBRaw240227(pathFldBands, bandsCodes)
    nrRaw, ncRaw = imgRaw.shape[0:2]
    if nrRaw == nr and ncRaw == nc:  # could be valid in all cases
        csRaw, rsRaw = GeotranformXY2CR240227(geotransform, xs, ys)
        HRaw = RSST2H240227(FindRSST01240227(csRaw, rsRaw, cs, rs))
        imgGeo = cv2.warpAffine(imgRaw, HRaw[:2, :], (nc, nr))
    else:
        color = [0, 0, 0]
        imgRaw[:margin, :, :], imgRaw[-margin:, :, :], imgRaw[:, :margin, :], imgRaw[:, -margin:, :] = color, color, color, color
        imgGeo = np.zeros((nr, nc, 3))
        csRaw, rsRaw = GeotranformXY2CR240227(geotransform, xs, ys)
        csIARaw, rsIARaw, wsARaw = CR2CRIntegerAroundAndWeights(csRaw, rsRaw)
        csIARaw, rsIARaw = np.clip(csIARaw, 0, ncRaw-1), np.clip(rsIARaw, 0, nrRaw-1)
        wsAAuxRaw = np.asarray([wsARaw, wsARaw, wsARaw])
        for corner in range(4):
            auxTMP = imgRaw[rsIARaw[:, corner], csIARaw[:, corner], :] * np.transpose(wsAAuxRaw[:, :, corner])
            imgGeo[rs, cs, :] = imgGeo[rs, cs, :] + auxTMP
        imgGeo = imgGeo.astype(np.uint8)
    return imgGeo
def PathFldBands2ImgRGBRaw240227(pathFldBands, bandsCodes):  # lm:2023-06-07; lr:2025-02-10
    assert len(bandsCodes) == 3
    auxRGB = []
    for bandCode in bandsCodes:  # for RGB respectively
        fnBand = [item for item in os.listdir(pathFldBands) if item.endswith('.{:}.tif'.format(bandCode))][0]
        datasetH = gdal.Open(os.path.join(pathFldBands, fnBand))
        assert datasetH.RasterCount == 1  # each file is 1-band
        bandH = datasetH.GetRasterBand(1)
        auxRGB.append(np.array([gdn.BandReadAsArray(bandH)]).astype('float32'))
    auxRGB = np.stack(auxRGB, axis=3)[0]
    img = cv2.cvtColor(auxRGB, cv2.COLOR_RGB2BGR)
    img = cv2.normalize(img, dst = None, alpha = 0, beta = 255, norm_type = cv2.NORM_MINMAX, dtype = cv2.CV_8UC3)
    if len(np.where(img == [0, 0, 0])[0]) / (img.shape[0]*img.shape[1]*3) * 100 > 20:  # WATCH OUT
        return None
    img = Img2ImgColorClipped240227(img)
    return img
def PathFldGroup2MaxSizeOfSubgroup240227(pathFldGroup):  # lm:2024-06-06; lr:2025-02-17
    maxSizeOfSubgroup = 0
    for fld in [item for item in PathFld2Flds240227(pathFldGroup) if item.startswith('subgroup_')]:
        pathNamesTxt = os.path.join(pathFldGroup, fld, 'names.txt')
        names = np.atleast_1d(np.loadtxt(pathNamesTxt, usecols=0, dtype=str)).tolist()
        maxSizeOfSubgroup = max(maxSizeOfSubgroup, len(names))
    return maxSizeOfSubgroup
def PathGivenRefTif2PathFldBands240227(pathRefTif, pathFldBands):  # lm:2025-01-03; lr:2025-02-10
    dataset0 = gdal.Open(pathRefTif)
    if dataset0 is None:
        return False
    try:
        os.makedirs(pathFldBands, exist_ok=True)
        nOfBands = dataset0.RasterCount
        assert nOfBands >= 1
        for posBand in range(1, 4):  # the bands start at 1 in gdal
            if nOfBands >= 3:
                band = dataset0.GetRasterBand(posBand)
            else:
                band = dataset0.GetRasterBand(1)
            pathTif1 = os.path.join(pathFldBands, '00_reference.B{:}.tif'.format(str(posBand)))
            driver = gdal.GetDriverByName('GTiff')
            dataset1 = driver.Create(pathTif1, dataset0.RasterXSize, dataset0.RasterYSize, 1, band.DataType)
            dataset1.SetGeoTransform(dataset0.GetGeoTransform())
            dataset1.SetProjection(dataset0.GetProjection())
            band1 = dataset1.GetRasterBand(1)
            band1.WriteArray(band.ReadAsArray())
            band1.SetDescription(band.GetDescription())
            dataset1 = None
        dataset0 = None
    except Exception:
        RmTreeIfExists(pathFldBands)
        return False
    if not AreAllMetadataInFldEqual240227(pathFldBands):
        shutil.rmtree(pathFldBands)
        return False
    return True
def PathImgOrImg2Img(img): # lm:2022-05-14; lr:2025-02-01
    if isinstance(img, str):
        img = cv2.imread(img)
    else:
        pass
    return img
def PathNpzOfPairs2DInfo240227(pathNpz, options={}):  # lm:2024-04-24; lr:2025-02-12
    keys, defaultValues = ['nBoxes', 'par'], None
    options = CompleteADictionary(options, keys, defaultValues)
    dInfo = {}
    pathFldXXViewName0, fnNpz1 = os.path.split(pathNpz)
    dInfo['pathFldXXViewName0'], dInfo['fnNpz1'], dInfo['name1'] = pathFldXXViewName0, fnNpz1, os.path.splitext(fnNpz1)[0]
    pathFldXXView, name0 = os.path.split(pathFldXXViewName0)
    dInfo['pathFldXXView'], dInfo['name0'] = pathFldXXView, name0
    pathFldXX, viewCode = os.path.split(pathFldXXView)
    dInfo['pathFldXX'], dInfo['viewCode'] = pathFldXX, viewCode
    pathFldScratch, fldXX = os.path.split(pathFldXX)
    dInfo['pathFldScratch'], dInfo['fldXX'] = pathFldScratch, fldXX
    pathFldMain, scratch = os.path.split(pathFldScratch)
    assert scratch == 'scratch'
    dInfo['pathFldMain'] = pathFldMain
    dInfo['datetime'] = datetime.datetime.fromtimestamp(os.path.getmtime(pathNpz))
    data = np.load(pathNpz)
    assert set(data) == set(['nc0', 'nr0', 'nc1', 'nr1', 'cs0', 'rs0', 'cs1', 'rs1', 'ers'])  # check for readability
    dInfo.update(data)
    nc0, nr0, cs0, rs0, cs1, rs1, ers = [data[item] for item in ['nc0', 'nr0', 'cs0', 'rs0', 'cs1', 'rs1', 'ers']]
    assert all([len(item) == len(cs0) for item in [cs0, rs0, cs1, rs1, ers]])  # check for readability
    if options['nBoxes'] is not None:
        nOfCBands, nOfRBands = NOfBands240227(options['nBoxes'], nc0, nr0)
        dInfo['nOfCBands'], dInfo['nOfRBands'] = nOfCBands, nOfRBands
        nOfFBoxes = Pixels2NOfFBoxes240227(nc0, nr0, nOfCBands, nOfRBands, cs0, rs0)
        dInfo['nOfFBoxes'] = nOfFBoxes
    if options['par'] is not None:
        res = options['par']['resCode2Res'][viewCode.split('_')[0]]
        dInfo['res'] = res
    return dInfo
def PathTif02PathTif1WithGeotransform1240227(pathTif0, geotransform1, pathTif1):  # lm:2025-01-03; lr:2025-02-11
    dataset0 = gdal.Open(pathTif0, gdal.GA_ReadOnly)
    if dataset0 is None:
        return None
    assert dataset0.RasterCount == 1
    band = dataset0.GetRasterBand(1)
    os.makedirs(os.path.dirname(pathTif1), exist_ok=True)
    driver = gdal.GetDriverByName("GTiff")
    dataset1 = driver.Create(pathTif1, dataset0.RasterXSize, dataset0.RasterYSize, 1, band.DataType)
    dataset1.SetGeoTransform(tuple(geotransform1))
    dataset1.SetProjection(dataset0.GetProjection())
    band1 = dataset1.GetRasterBand(1)
    band1.WriteArray(band.ReadAsArray())
    dataset0, dataset1 = None, None
    return None
def PathTif2Metadata240227(pathTif):  # lm:2025-01-03; lr:2025-02-12
    ds = gdal.Open(pathTif)
    metadata = {'width': ds.RasterXSize, 'height': ds.RasterYSize, 'bands': ds.RasterCount, 'driver': ds.GetDriver().LongName, 'projection': ds.GetProjection(), 'geotransform': ds.GetGeoTransform()}
    return metadata
def PathZip2PathFldBands240227(pathZip, pathFldBands, par):  # lm:2025-02-10; lr:2025-02-10
    fnZip = os.path.split(pathZip)[1]
    satCode = fnZip.split('_')[0]  # IMP*
    if not (fnZip.endswith('.zip') and satCode in par['satCode2ICode']):  # satCode 00 not allowed for zip
        RmTreeIfExists(pathFldBands)
        return False
    try:
        with zipfile.ZipFile(pathZip, 'r') as zip_ref:
            zip_ref.extractall(pathFldBands)
    except Exception:
        RmTreeIfExists(pathFldBands)
        return False
    for fn in os.listdir(pathFldBands):
        if not fn.endswith('.tif') or (satCode in ['L7', 'L8', 'L9'] and '.B8.' in fn):  # IMP*; panchromatic
            os.remove(os.path.join(pathFldBands, fn))
    for fn in os.listdir(pathFldBands):
        os.rename(os.path.join(pathFldBands, fn), os.path.join(pathFldBands, '{:}_{:}'.format(satCode, fn)))
    for bandCode in par['satCode2RGBBandsCodes'][satCode]:  # band by band
        if not any('.{:}.'.format(bandCode) in item for item in os.listdir(pathFldBands)):
            shutil.rmtree(pathFldBands)
            return False
    RGBCodeTMP = par['satCode2RGBBandsCodes'][satCode][0]
    fnTif0 = [item for item in os.listdir(pathFldBands) if '.{:}.'.format(RGBCodeTMP) in item and item.endswith('.tif')][0]
    metadata0 = PathTif2Metadata240227(os.path.join(pathFldBands, fnTif0))
    for fn in os.listdir(pathFldBands):
        try:
            metadata = PathTif2Metadata240227(os.path.join(pathFldBands, fn))
            assert AreTwoMetadataEqual240227(metadata, metadata0)
        except Exception:
            os.remove(os.path.join(pathFldBands, fn))
    assert AreAllMetadataInFldEqual240227(pathFldBands)
    return True
def Pixels2BandPositionsInGrid240227(nc, nr, nOfCBands, nOfRBands, cs, rs):  # lm:2024-01-08; lr:2025-02-08
    if len(cs) == 0:
        bandCs, bandRs, bandGs = [np.asarray([]) for _ in range(3)]
    else:
        cs = np.clip(cs, 0, nc-1.e-11)
        rs = np.clip(rs, 0, nr-1.e-11)
        bandCs = np.floor(cs * nOfCBands / nc).astype(int)  # bandCs: cs=0 -> bandCs=0; cs=nc -> bandCs=nOfCBands (never reached)
        bandRs = np.floor(rs * nOfRBands / nr).astype(int)  # bandRs: rs=0 -> bandRs=0; rs=nr -> bandRs=nOfRBands (never reached)
        bandGs = bandCs * nOfRBands + bandRs  # bandGs: global counter
    return bandCs, bandRs, bandGs
def Pixels2NOfFBoxes240227(nc, nr, nOfCBands, nOfRBands, cs, rs):  # lm:2024-04-10; lr:2025-02-12
    if len(cs) == 0:
        nOfFBoxes = 0
    else:
        bandGs = Pixels2BandPositionsInGrid240227(nc, nr, nOfCBands, nOfRBands, cs, rs)[-1]
        nOfFBoxes = len(set(bandGs))
    return nOfFBoxes
def PosGroup2GroupCode240227(posGroup):  # lm:2024-03-07; lr:2025-02-11
    groupCode = 'group_{:}'.format(Integer2Str240227(posGroup))
    return groupCode
def PosSubgroup2SubgroupCode240227(posSubgroup):  # lm:2024-05-24; lr:2025-02-11
    subgroupCode = 'subgroup_{:}'.format(Integer2Str240227(posSubgroup))
    return subgroupCode
def Poss0AndPoss1InFind2DTransform(n): # lm:2022-07-11; lr:2025-02-08
    poss0 = [2*pos+0 for pos in range(n)]
    poss1 = [2*pos+1 for pos in range(n)]
    return poss0, poss1
def PurgePixelsRANSAC240227(dt01, cs0, rs0, cs1, rs1, ers, nc0, nr0, fraction, errorC, rotation, nOfCBands, nOfRBands, nOfIter=10):  # lm:2024-04-17; lr:2025-02-08
    assert np.max(cs0) - np.min(cs0) >= fraction * nc0 and np.max(rs0) - np.min(rs0) >= fraction * nr0
    assert len(cs0) == len(rs0) == len(cs1) == len(rs1) >= 4  # (comes from "if len(ers) < 4:  # IMP*; BLOCKED")
    parRANSAC = {'pOutlier': 0.5, 'pDesired': 0.9999, 'nOfPoints': 2, 'errorC': errorC}  # IMP*; 0.9999 can be increased to try to get better results
    possGood = FindGoodPossRANSAC240227(dt01, cs0, rs0, cs1, rs1, parRANSAC, nc0, nr0, fraction, rotation, nOfCBands, nOfRBands, nOfIter=nOfIter)
    if len(possGood) == 0:  # IMP*
        cs0, rs0, cs1, rs1, ers = [np.asarray([]) for item in range(5)]
    else:
        cs0, rs0, cs1, rs1, ers = [item[possGood] for item in [cs0, rs0, cs1, rs1, ers]]
        assert np.max(cs0) - np.min(cs0) >= fraction * nc0 and np.max(rs0) - np.min(rs0) >= fraction * nr0  # check for readability (ensured in FindGoodPossRANSAC240227)
    return cs0, rs0, cs1, rs1, ers
def RSST01AndRSST12ToRSST02240227(RSST01, RSST12):  # lm:2024-02-29; lr:2025-02-08
    H01, H12 = [RSST2H240227(item) for item in [RSST01, RSST12]]
    H02 = H01AndH12ToH02240227(H01, H12)
    RSST02 = H2RSST240227(H02)
    return RSST02
def RSST2H240227(RSST):  # lm:2024-02-29; lr:2025-02-08
    H = np.eye(3)
    H[0, 0], H[0, 1], H[0, 2] = +RSST[0], RSST[1], RSST[2]
    H[1, 0], H[1, 1], H[1, 2] = -RSST[1], RSST[0], RSST[3]
    return H
def RSSTNcNr2RMSE2240227(RSST, nc, nr):  # lm:2024-12-29; lr:2025-02-03
    dtcosalpha, dtsinalpha, dc, dr = RSST
    dt = np.sqrt(dtcosalpha ** 2 + dtsinalpha ** 2)  # dt = dilatation
    cosalpha, sinalpha = [item / dt for item in [dtcosalpha, dtsinalpha]]
    ec, es = 2 - 2 * cosalpha, 2 * sinalpha
    kc, kr = -ec * dc - es * dr, es * dc - ec * dr
    rmse2 = 0
    rmse2 = rmse2 + ((nc - 1) * (2 * nc - 1) / 6 + (nr - 1) * (2 * nr - 1) / 6) * ec
    rmse2 = rmse2 + (nc - 1) / 2 * kc
    rmse2 = rmse2 + (nr - 1) / 2 * kr
    rmse2 = rmse2 + dc ** 2 + dr ** 2
    rmse2 = max(0, rmse2)  # to avoid negative epsilons
    return rmse2
def ReadTimeStampLog(pathFldXXView):  # lm:2025-01-03; lr:2025-02-10
    if not os.path.exists(pathFldXXView):
        date = '0' * 17
    else:
        fnsTMP = sorted([item for item in os.listdir(pathFldXXView) if len(item) == 25 and item.startswith('000_') and item.endswith('.log')])
        if len(fnsTMP) == 0:
            date = '0' * 17
        elif len(fnsTMP) == 1:
            date = fnsTMP[0][4:21]
        elif len(fnsTMP) > 1:  # this should not happen
            date = fnsTMP[-1][4:21]
            for fnTMP in fnsTMP[:-1]:
                os.remove(os.path.join(pathFldXXView, fnTMP))
    assert len(date) == 17 and int(date) > -1
    return date
def ReverseTwoStringsWithTo240227(str01):  # lm:2024-03-13; lr:2025-02-12
    str0, str1 = SplitTwoStringsWithTo240227(str01)
    str10 = MergeTwoStringsWithTo240227(str1, str0)
    return str10
def RmEmptyFldsRecursively(pathFld, exceptions):  # lm:2024-11-27; lr:2025-02-01
    for root, flds, _ in os.walk(pathFld, topdown=False):
        for fld in flds:
            pathFldH = os.path.join(root, fld)
            if not os.listdir(pathFldH) and all(item not in pathFldH for item in exceptions):
                os.rmdir(pathFldH)
    if not os.listdir(pathFld) and all(item not in pathFld for item in exceptions):
        os.rmdir(pathFld)
    return None
def RmFileIfExists(pathFile):  # lm:2024-06-18; lr:2025-02-01
    if os.path.exists(pathFile) and os.path.isfile(pathFile):
        os.remove(pathFile)
    return None
def RmTreeIfExists(pathFld):  # lr:2024-06-19; lr:2025-02-03
    if os.path.exists(pathFld) and os.path.isdir(pathFld):
        shutil.rmtree(pathFld)
    return None
def SIFTKeypoints(img, options={}):  # lm:2024-12-09; lr:2025-02-03
    keys, defaultValues = ['mask', 'nOfFeatures'], [None, 5000]
    options = CompleteADictionary(options, keys, defaultValues)
    try:
        img = PathImgOrImg2Img(img)
        nr, nc = img.shape[0:2]
        sift = cv2.SIFT_create(nfeatures=options['nOfFeatures'])
        if options['mask'] is not None:
            kps, des = sift.detectAndCompute(img, None)  # MISSING; WATCH OUT
        else:
            kps, des = sift.detectAndCompute(img, None)
        assert len(kps) == len(des) > 0
        ctrl = True
    except Exception:
        nc, nr, kps, des, ctrl = None, None, None, None, False
    return nc, nr, kps, des, ctrl
def SIFTMatches(kps1, des1, kps2, des2, options={}):  # lm:2025-02-11; lr:2025-02-11
    keys, defaultValues = ['erMaximum', 'rThres', 'nOfStd', 'resolution1', 'resolution2'], [1.e+11, 0.75, 2., 1., 1.]
    options = CompleteADictionary(options, keys, defaultValues)
    cs1, rs1 = [np.asarray([item.pt[pos] for item in kps1]) for pos in [0, 1]]
    cs2, rs2 = [np.asarray([item.pt[pos] for item in kps2]) for pos in [0, 1]]
    bf = cv2.BFMatcher(cv2.NORM_L2)
    matches12 = bf.knnMatch(des1, des2, k=2)
    matches21 = bf.knnMatch(des2, des1, k=2)
    good_matches12 = []
    for m, n in matches12:
        if m.distance < options['rThres'] * n.distance:
            good_matches12.append(m)
    good_matches21 = []
    for m, n in matches21:
        if m.distance < options['rThres'] * n.distance:
            good_matches21.append(m)
    matches12 = good_matches12
    matches21 = good_matches21
    poss1 = [match.queryIdx for match in matches12] + [match.trainIdx for match in matches21]
    poss2 = [match.trainIdx for match in matches12] + [match.queryIdx for match in matches21]
    cs1, rs1, cs2, rs2 = cs1[poss1], rs1[poss1], cs2[poss2], rs2[poss2]
    ers = np.asarray([match.distance for match in matches12] + [match.distance for match in matches21])
    cs1, rs1, cs2, rs2, ers = np.unique(np.asarray([cs1, rs1, cs2, rs2, ers]), axis=1)  # IMP* interesting
    if len(cs1) > 0:
        xs1, xs2 = cs1 * options['resolution1'], cs2 * options['resolution2']
        ys1, ys2 = rs1 * options['resolution1'], rs2 * options['resolution2']
        ds = np.sqrt((xs1 - xs2) ** 2 + (ys1 - ys2) ** 2)
        possGood = np.where((ers < options['erMaximum']) & (ds < np.mean(ds) + options['nOfStd'] * np.std(ds) + 1.e-8))[0]
        possGood = np.where((ers < options['erMaximum']) & (ds < np.mean(ds) + options['nOfStd'] * np.std(ds) + 1.e-8))[0]
        cs1, rs1, cs2, rs2, ers = [item[possGood] for item in [cs1, rs1, cs2, rs2, ers]]
    return cs1, rs1, cs2, rs2, ers
def SelectPossInGrid240227(nc, nr, nOfCBands, nOfRBands, cs, rs, ers, fraction):  # lm:2025-01-07; lr:2025-02-08
    assert np.max(cs) - np.min(cs) >= fraction * nc and np.max(rs) - np.min(rs) >= fraction * nr
    bandCs, bandRs, bandGs = Pixels2BandPositionsInGrid240227(nc, nr, nOfCBands, nOfRBands, cs, rs)
    nOfFBoxes = Pixels2NOfFBoxes240227(nc, nr, nOfCBands, nOfRBands, cs, rs)
    possSelected = []
    for bandGU in set(bandGs):  # unique global positions (characterizing boxes)
        possOfBox = np.where(bandGs == bandGU)[0]  # list of global positions in the box
        if len(possOfBox) == 1:
            possSelected.append(possOfBox[0])
        else:
            assert np.std(bandCs[possOfBox]) < 1.e-11 and np.std(bandRs[possOfBox]) < 1.e-11  # avoidable check
            possSelected.append(possOfBox[np.argmin(ers[possOfBox])])  # IMP*
    csS, rsS = [item[possSelected] for item in [cs, rs]]  # S = Selected
    assert Pixels2NOfFBoxes240227(nc, nr, nOfCBands, nOfRBands, csS, rsS) == nOfFBoxes  # avoidable check
    if np.max(csS) - np.min(csS) >= fraction * nc and np.max(rsS) - np.min(rsS) >= fraction * nr:
        return possSelected
    possSelected0 = []
    for bandGU in set(bandGs):  # unique global positions (characterizing boxes)
        possOfBox = np.where(bandGs == bandGU)[0]  # list of global positions in the box
        if len(possOfBox) == 1:
            possSelected0.append(possOfBox[0])
        else:
            assert np.std(bandCs[possOfBox]) < 1.e-11 and np.std(bandRs[possOfBox]) < 1.e-11  # avoidable check
            csBox, rsBox = [item[possOfBox] for item in [cs, rs]]
            condCMin = np.isclose(np.min(cs), np.min(csBox))  # the point with min cs is in the box
            condCMax = np.isclose(np.max(cs), np.max(csBox))  # the point with max cs is in the box
            condRMin = np.isclose(np.min(rs), np.min(rsBox))  # the point with min rs is in the box
            condRMax = np.isclose(np.max(rs), np.max(rsBox))  # the point with max rs is in the box
            if condCMin or condCMax or condRMin or condRMax:
                if condCMin:
                    possSelected0.append([item for item in possOfBox if np.isclose(cs[item], np.min(csBox))][0])
                if condCMax:  # IMP*; not elif
                    possSelected0.append([item for item in possOfBox if np.isclose(cs[item], np.max(csBox))][0])
                if condRMin:  # IMP*; not elif
                    possSelected0.append([item for item in possOfBox if np.isclose(rs[item], np.min(rsBox))][0])
                if condRMax:  # IMP*; not elif
                    possSelected0.append([item for item in possOfBox if np.isclose(rs[item], np.max(rsBox))][0])
            else:
                possSelected0.append(possOfBox[np.argmin(ers[possOfBox])])  # IMP*
    assert len(possSelected0) >= nOfFBoxes
    csS0, rsS0 = [item[possSelected0] for item in [cs, rs]]
    assert Pixels2NOfFBoxes240227(nc, nr, nOfCBands, nOfRBands, csS0, rsS0) == nOfFBoxes  # avoidable check for readability
    assert np.max(csS0) - np.min(csS0) >= fraction * nc and np.max(rsS0) - np.min(rsS0) >= fraction * nr  # avoidable check for readability
    counterH = 0
    while True:
        counterH = counterH + 1
        possSelected, bandsGUs = [], []
        for posSelected0 in random.sample(possSelected0, len(possSelected0)):  # think OK
            if bandGs[posSelected0] not in bandsGUs:  # think OK
                possSelected.append(posSelected0)
                bandsGUs.append(bandGs[posSelected0])
        csS, rsS = [item[possSelected] for item in [cs, rs]]
        assert Pixels2NOfFBoxes240227(nc, nr, nOfCBands, nOfRBands, csS, rsS) == nOfFBoxes  # avoidable check for readability
        if np.max(csS) - np.min(csS) >= fraction * nc and np.max(rsS) - np.min(rsS) >= fraction * nr:
            break
        if counterH > 50:  # IMP*
            possSelected = []
            break
    return possSelected
def SplitTwoStringsWithTo240227(str01):  # lm:2024-03-07; lr:2025-02-12
    str0, str1 = str01.split('_to_')
    assert MergeTwoStringsWithTo240227(str0, str1) == str01  # useful check for readability
    return str0, str1
def Str0Str12Cs0Rs0Cs1Rs1240227(str0, str1, dPixelsConns):  # lm:2024-04-18; lr:2025-02-11
    str01 = MergeTwoStringsWithTo240227(str0, str1)
    str10 = MergeTwoStringsWithTo240227(str1, str0)
    if str01 in dPixelsConns.keys():
        assert str0 < str1  # check for readability
        cs0, rs0, cs1, rs1 = [dPixelsConns[str01][item] for item in ['cs0', 'rs0', 'cs1', 'rs1']]
    elif str10 in dPixelsConns.keys():
        assert str1 < str0  # check for readability
        cs1, rs1, cs0, rs0 = [dPixelsConns[str10][item] for item in ['cs0', 'rs0', 'cs1', 'rs1']]  # IMP*; from dPixelsConns[str10]
    else:
        assert False
    return cs0, rs0, cs1, rs1
def Str2Float240227(string):  # lm:2024-03-07; lr:2025-02-10
    flt = int(string) / 100  # IMP*: 2 decimals
    return flt
def Str2Integer240227(string):  # lm:2024-03-14; lr:2025-02-10
    integer = int(string)
    return integer
def TwoMatrices2Corr240227(mus, mvs):  # lm:2024-03-18; lr:2025-02-11
    assert mus.shape == mvs.shape
    us, vs = [np.reshape(item, -1) for item in [mus, mvs]]
    assert len(us) == len(vs) > 1
    corr = np.sum((us - np.mean(us)) * (vs - np.mean(vs))) / (len(us) - 1)
    return corr
def UpdateLRemainingPairs240227(lRemainingPairs, lPWalksAsPairs):  # lm:2024-03-03; lr:2025-02-08
    lPWalksAsLists = LPWalksAsPairs2LPWalksAsLists240227(lPWalksAsPairs)
    lRemainingPairsOld = lRemainingPairs
    lRemainingPairsNew = []
    for pairOld in lRemainingPairsOld:
        if any([set(pairOld) <= set(item) for item in lPWalksAsLists]):
            pass  # IMP*; connection is satisfied
        else:
            lRemainingPairsNew.append(pairOld)
    lRemainingPairs = lRemainingPairsNew
    return lRemainingPairs
def ViewCode2N240227(pathFldScratch, viewCode):  # lm:2024-06-05; lr:2025-02-12
    names = ViewCode2Names240227(pathFldScratch, viewCode)
    N = len(names) * (len(names) - 1) / 2
    assert np.isclose(N, int(np.round(N)))
    N = int(np.round(N))
    return N
def ViewCode2Names240227(pathFldScratch, viewCode):  # lm:2024-11-27; lr:2025-02-12
    pathFld01View = os.path.join(pathFldScratch, '01_images', viewCode)
    names = sorted([os.path.splitext(item)[0] for item in PathFld2Fns240227(pathFld01View, ext='.png')])  # IMP*; png
    return names
def ViewCodesInScratch240227(pathFldScratch):  # lm:2025-01-02; lr:2025-02-10
    viewCodes = set()
    for fldXX in PathFld2Flds240227(pathFldScratch):
        pathFldXX = os.path.join(pathFldScratch, fldXX)
        viewCodes.update(PathFld2Flds240227(pathFldXX))
    viewCodes = sorted(viewCodes)
    return viewCodes
def WriteGroups240227(pathFldGroup, names, mConns, group):  # lm:2024-06-06; lr:2025-02-11
    namesGroup = [names[item] for item in group]  # IMP*
    assert namesGroup == sorted(namesGroup)  # check for readability
    pathNamesTxt = os.path.join(pathFldGroup, 'names.txt')
    if os.path.exists(pathNamesTxt):  # if exists, check
        namesGroupA = np.atleast_1d(np.loadtxt(pathNamesTxt, usecols=0, dtype=str)).tolist()
        assert namesGroupA == namesGroup
    else:
        fileout = open(pathNamesTxt, 'w')
        for nameGroup in namesGroup:
            fileout.write('{:}\n'.format(nameGroup))
        fileout.close()
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
    mConnsGroup = mConns[group, :][:, group]  # symmetric
    dConnsFMGroup = NamesAndMConns2DConns240227(namesGroup, mConnsGroup)
    assert set(dConnsFMGroup.keys()) == set(namesGroup)  # check for readability
    pathDConnsFM = os.path.join(pathFldGroup, 'connections_FM.json')
    if os.path.exists(pathDConnsFM):  # if exists, check
        with open(pathDConnsFM, 'r') as f:
            dConnsFMGroupA = json.load(f)
        assert sorted(dConnsFMGroupA.keys()) == sorted(dConnsFMGroup.keys())
        assert all(sorted(dConnsFMGroupA[key]) == sorted(dConnsFMGroup[key]) for key in dConnsFMGroup.keys())
    else:
        fileout = open(pathDConnsFM, 'w')
        json.dump(dConnsFMGroup, fileout, indent=2)
        fileout.close()
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
    nameRefGroup = DConns2NameRef240227(dConnsFMGroup)
    pathNameRef = os.path.join(pathFldGroup, '{:}.ref'.format(nameRefGroup))
    if not os.path.exists(pathNameRef):
        fileout = open(pathNameRef, 'w')
        fileout.close()
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
    assert len(PathFld2Fns240227(pathFldGroup, ext='.ref')) == 1  # check in case there was something
    return namesGroup, mConnsGroup, nameRefGroup
def WritePathNpz02240227(pathNpz02, methodFM, nc0Ori, nr0Ori, nc1Ori, nr1Ori, kps0Sca, des0Sca, kps1Sca, des1Sca, scaFM):  # lm:2024-06-13; lr:2025-01-08
    if methodFM == 'ORB':
        cs0Sca, rs0Sca, cs1Sca, rs1Sca, ers = ORBMatches(kps0Sca, des0Sca, kps1Sca, des1Sca, options={})
    elif methodFM == 'SIFT':
        cs0Sca, rs0Sca, cs1Sca, rs1Sca, ers = SIFTMatches(kps0Sca, des0Sca, kps1Sca, des1Sca, options={})
    else:
        assert False
    if len(ers) < 4:  # IMP*; BLOCKED
        doesWrite = False
    else:
        cs0Ori, rs0Ori, cs1Ori, rs1Ori = [item / scaFM for item in [cs0Sca, rs0Sca, cs1Sca, rs1Sca]]
        os.makedirs(os.path.dirname(pathNpz02), exist_ok=True)
        np.savez(pathNpz02, nc0=nc0Ori, nr0=nr0Ori, nc1=nc1Ori, nr1=nr1Ori, cs0=cs0Ori, rs0=rs0Ori, cs1=cs1Ori, rs1=rs1Ori, ers=ers)
        WriteTimeStampLog(os.path.split(os.path.split(pathNpz02)[0])[0])  # IMP*
        doesWrite = True
    return doesWrite
def WriteSubgroupsFlds240227(pathFldGroup, nc4SatCode, nr4SatCode):  # lm:2024-12-29; lr:2025-02-11
    fld06 = pathFldGroup.split(os.sep)[-3]
    conDegree = FldXX2DInfo240227(fld06)['conDegree']
    pathTMP = os.path.join(pathFldGroup, 'names.txt')
    namesGroup = np.atleast_1d(np.loadtxt(pathTMP, usecols=0, dtype=str)).tolist()
    nameRefGroup = [os.path.splitext(item)[0] for item in os.listdir(pathFldGroup) if item.endswith('.ref')][0]  # it is not 00_reference even if 00_reference in nameGroup
    assert nameRefGroup in namesGroup
    fnsRSSTjson = [item for item in os.listdir(pathFldGroup) if item.startswith('RSST_') and item.endswith('.json')]
    if len(fnsRSSTjson) == 0:
        return None
    assert len(fnsRSSTjson) == 1
    pathRSST = os.path.join(pathFldGroup, fnsRSSTjson[0])
    with open(pathRSST, 'r') as f:
        dRSSTs = json.load(f)
    pathDConnsRecov = pathRSST.replace('RSST_', 'connections_')
    with open(pathDConnsRecov, 'r') as f:
        dConnsRecov = json.load(f)
    mConnsRecov = DConns2NamesAndMConns240227(dConnsRecov)[1]
    subgroups = MConns2Groups240227(mConnsRecov, conDegree)  # local counter in group
    for posSubgroup, subgroup in enumerate(subgroups):
        subgroupCode = PosSubgroup2SubgroupCode240227(posSubgroup)
        pathFldSubgroup = os.path.join(pathFldGroup, subgroupCode)
        os.makedirs(pathFldSubgroup, exist_ok=True)
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
        namesSubgroup = [namesGroup[pos] for pos in subgroup]
        assert namesSubgroup == sorted(namesSubgroup)  # for readability
        pathNamesTxt = os.path.join(pathFldSubgroup, 'names.txt')
        assert not os.path.exists(pathNamesTxt)
        fileout = open(pathNamesTxt, 'w')
        for nameSubgroup in namesSubgroup:
            fileout.write('{:}\n'.format(nameSubgroup))
        fileout.close()
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
        if len(subgroup) == 1:
            nameRefSubgroup = namesSubgroup[0]
        elif len([item for item in namesSubgroup if 'reference' in item]) > 0:  # IMP*
            nameRefSubgroup = [item for item in namesSubgroup if 'reference' in item][0]
        else:  # IMP*
            nameRefSubgroup = CenteredSubgroupName(dRSSTs, nameRefGroup, namesSubgroup, nc4SatCode, nr4SatCode)
        pathNameRef = os.path.join(pathFldGroup, subgroupCode, '{:}.ref'.format(nameRefSubgroup))
        assert not os.path.exists(pathNameRef)  # for readability
        fileout = open(pathNameRef, 'w')
        fileout.close()
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
        if len(subgroup) <= conDegree + 1:  # IMP*; scape
            continue
        dRSSTsRefSubgroup = DRSSTsRefNew240227(dRSSTs, namesSubgroup, nameRefSubgroup)
        pathRSST = os.path.join(pathFldGroup, subgroupCode, 'RSST.json')
        assert not os.path.exists(pathRSST)
        fileout = open(pathRSST, 'w')
        json.dump(dRSSTsRefSubgroup, fileout, indent=2)
        fileout.close()
        WriteTimeStampLog(os.path.split(pathFldGroup)[0])
    return None
def WriteTimeStampLog(pathFldXXView):  # lm:2024-12-29; lr:2025-02-10
    os.makedirs(pathFldXXView, exist_ok=True)
    fnsTMP = [item for item in os.listdir(pathFldXXView) if len(item) == 25 and item.startswith('000_') and item.endswith('.log')]
    for fnTMP in fnsTMP:
        os.remove(os.path.join(pathFldXXView, fnTMP))
    fn = '000_{:}.log'.format(Datetime2Date17(datetime.datetime.now()))
    assert len(fn) == 25  # avoidable check
    pathFile = os.path.join(pathFldXXView, fn)
    fileout = open(pathFile, 'w')
    fileout.close()
    return None
def XsYsA2XsYsB240227(projectionA, projectionB, xsA, ysA):  # lm:2024-06-08; lr:2025-02-11
    if projectionB == projectionA:
        xsB, ysB = xsA, ysA
    else:
        srs_A = osr.SpatialReference()
        srs_A.ImportFromWkt(projectionA)
        srs_B = osr.SpatialReference()
        srs_B.ImportFromWkt(projectionB)
        transform = osr.CoordinateTransformation(srs_A, srs_B)
        xsysB = [transform.TransformPoint(xA, yA) for xA, yA in zip(xsA, ysA)]
        xsB = np.asarray([item[0] for item in xsysB])
        ysB = np.asarray([item[1] for item in xsysB])
    return xsB, ysB
