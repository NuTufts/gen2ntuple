
import os,sys,argparse

import ROOT as rt

from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow

from math import sqrt as sqrt
from math import acos as acos
from math import pi

from array import array
import numpy as np

from event_weighting.event_weight_helper import SumPOT, Weights
from helpers.larflowreco_ana_funcs import *

import torch
from torch import nn
from torch.utils.data import DataLoader
import torchvision.transforms as transforms

sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__)))+'/prongCNN/models')

from models_instanceNorm_reco_2chan_multiTask import ResBlock, ResNet34
from datasets_reco_5ClassHardLabel_multiTask import mean, std

parser = argparse.ArgumentParser("Validate Vertex Reco")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list")
parser.add_argument("-w", "--weightfile", type=str, default="none", help="weights file (pickled python dict)")
parser.add_argument("-m", "--model_path", type=str, required=True, help="path to prong CNN checkpoint file")
parser.add_argument("-d", "--device", type=str, default="cpu", help="gpu/cpu device")
parser.add_argument("-mc", "--isMC", help="running over MC input", action="store_true")
parser.add_argument("-o", "--outfile", type=str, default="prepare_selection_test_output.root", help="output file name")
parser.add_argument("--multiGPU", action="store_true", help="use multiple GPUs")
parser.add_argument("--oldVtxBranch", help="use nufitted_v instead of nuvetoed_v for old reco", action="store_true")
args = parser.parse_args()

if args.isMC and args.weightfile=="none":
  sys.exit("Must supply weight file for MC input. Exiting...")

reco2Tag = "merged_dlana_"
if args.isMC:
  reco2Tag = "merged_dlreco_"
files = getFiles(reco2Tag, args.files, args.truth)


def addClusterCharge(iolcv, cluster, vertexPixels, vertexCharge, threshold):
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()
  clusterPixels = []
  clusterCharge = 0.
  for hit in cluster:
    for p in range(3):
      row = (hit.tick - 2400)//6
      pixel = [ p, hit.tick, hit.targetwire[p] ]
      pixVal = image2Dvec[p].pixel(row, hit.targetwire[p])
      if pixVal < threshold:
        continue
      if pixel not in clusterPixels:
        clusterPixels.append(pixel)
        clusterCharge += pixVal
      if pixel not in vertexPixels:
        vertexPixels.append(pixel)
        vertexCharge += pixVal
  return clusterCharge, vertexPixels, vertexCharge


def makeImage(prong_vv):
  plane0pix_row = np.zeros(prong_vv[0].size(), dtype=int)
  plane0pix_col = np.zeros(prong_vv[0].size(), dtype=int)
  plane0pix_val = np.zeros(prong_vv[0].size(), dtype=float)
  plane1pix_row = np.zeros(prong_vv[1].size(), dtype=int)
  plane1pix_col = np.zeros(prong_vv[1].size(), dtype=int)
  plane1pix_val = np.zeros(prong_vv[1].size(), dtype=float)
  plane2pix_row = np.zeros(prong_vv[2].size(), dtype=int)
  plane2pix_col = np.zeros(prong_vv[2].size(), dtype=int)
  plane2pix_val = np.zeros(prong_vv[2].size(), dtype=float)
  raw_plane0pix_row = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_col = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_val = np.zeros(prong_vv[3].size(), dtype=float)
  raw_plane1pix_row = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_col = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_val = np.zeros(prong_vv[4].size(), dtype=float)
  raw_plane2pix_row = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_col = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_val = np.zeros(prong_vv[5].size(), dtype=float)
  for i, pix in enumerate(prong_vv[0]):
    plane0pix_row[i] = pix.row
    plane0pix_col[i] = pix.col
    plane0pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[1]):
    plane1pix_row[i] = pix.row
    plane1pix_col[i] = pix.col
    plane1pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[2]):
    plane2pix_row[i] = pix.row
    plane2pix_col[i] = pix.col
    plane2pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[3]):
    raw_plane0pix_row[i] = pix.row
    raw_plane0pix_col[i] = pix.col
    raw_plane0pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[4]):
    raw_plane1pix_row[i] = pix.row
    raw_plane1pix_col[i] = pix.col
    raw_plane1pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[5]):
    raw_plane2pix_row[i] = pix.row
    raw_plane2pix_col[i] = pix.col
    raw_plane2pix_val[i] = pix.val
  image = np.zeros((6,512,512))
  image[0, plane0pix_row, plane0pix_col] = plane0pix_val
  image[2, plane1pix_row, plane1pix_col] = plane1pix_val
  image[4, plane2pix_row, plane2pix_col] = plane2pix_val
  image[1, raw_plane0pix_row, raw_plane0pix_col] = raw_plane0pix_val
  image[3, raw_plane1pix_row, raw_plane1pix_col] = raw_plane1pix_val
  image[5, raw_plane2pix_row, raw_plane2pix_col] = raw_plane2pix_val
  image = torch.from_numpy(image).float()
  norm = transforms.Normalize(mean, std)
  image = norm(image).reshape(1,6,512,512)
  return torch.clamp(image, max=4.0)


def getPID(cnnClass):
    if cnnClass == 0:
        return 11
    if cnnClass == 1:
        return 22
    if cnnClass == 2:
        return 13
    if cnnClass == 3:
        return 211
    if cnnClass == 4:
        return 2212
    return 0


#other needed classes
sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()
prongvars = larflow.reco.NuSelProngVars()
wcoverlapvars = larflow.reco.NuSelWCTaggerOverlap()
flowTriples = larflow.prep.FlowTriples()

model = ResNet34(2, ResBlock, outputs=5)
if "cuda" in args.device and args.multiGPU:
  model = nn.DataParallel(model)
checkpoint = torch.load(args.model_path)
try:
  model.load_state_dict(checkpoint['model_state_dict'])
except:
  model.module.load_state_dict(checkpoint['model_state_dict'])
model.to(args.device)
model.eval()

outRootFile = rt.TFile(args.outfile, "RECREATE")

if args.isMC:
  potTree = rt.TTree("potTree","potTree")
  totPOT = array('f', [0.])
  totGoodPOT = array('f', [0.])
  potTree.Branch("totPOT", totPOT, 'totPOT/F')
  potTree.Branch("totGoodPOT", totGoodPOT, 'totGoodPOT/F')

eventTree = rt.TTree("EventTree","EventTree")
maxNTrks = 100
maxNShwrs = 100
fileid = array('i', [0])
run = array('i', [0])
subrun = array('i', [0])
event = array('i', [0])
xsecWeight = array('f', [0.])
trueNuE = array('f', [0.])
trueLepE = array('f', [0.])
trueLepPDG = array('i', [0])
trueMuContained = array('i', [0])
trueNuPDG = array('i', [0])
trueNuCCNC = array('i', [0])
nVertices = array('i', [0])
vtxX = array('f', [0.])
vtxY = array('f', [0.])
vtxZ = array('f', [0.])
vtxIsFiducial = array('i', [0])
vtxDistToTrue = array('f', [0.])
vtxBestComp = array('f', [0.])
vtxScore = array('f', [0.])
vtxFracHitsOnCosmic = array('f', [0.])
nTracks = array('i', [0])
trackNHits = array('i', maxNTrks*[0])
trackHitFrac = array('f', maxNTrks*[0.])
trackCharge = array('f', maxNTrks*[0.])
trackChargeFrac = array('f', maxNTrks*[0.])
trackRecoMuKE = array('f', maxNTrks*[0.])
trackRecoPrKE = array('f', maxNTrks*[0.])
trackClassified = array('i', maxNTrks*[0])
trackRecoPID = array('i', maxNTrks*[0])
trackRecoComp = array('f', maxNTrks*[0.])
nShowers = array('i', [0])
showerNHits = array('i', maxNShwrs*[0])
showerHitFrac = array('f', maxNShwrs*[0.])
showerCharge = array('f', maxNShwrs*[0.])
showerChargeFrac = array('f', maxNShwrs*[0.])
showerRecoPl0E = array('f', maxNShwrs*[0.])
showerRecoPl1E = array('f', maxNShwrs*[0.])
showerRecoPl2E = array('f', maxNShwrs*[0.])
showerClassified = array('i', maxNShwrs*[0])
showerRecoPID = array('i', maxNShwrs*[0])
showerRecoComp = array('f', maxNShwrs*[0])
eventTree.Branch("fileid", fileid, 'fileid/I')
eventTree.Branch("run", run, 'run/I')
eventTree.Branch("subrun", subrun, 'subrun/I')
eventTree.Branch("event", event, 'event/I')
eventTree.Branch("xsecWeight", xsecWeight, 'xsecWeight/F')
eventTree.Branch("trueNuE", trueNuE, 'trueNuE/F')
eventTree.Branch("trueLepE", trueLepE, 'trueLepE/F')
eventTree.Branch("trueLepPDG", trueLepPDG, 'trueLepPDG/I')
eventTree.Branch("trueMuContained", trueMuContained, 'trueMuContained/I')
eventTree.Branch("trueNuPDG", trueNuPDG, 'trueNuPDG/I')
eventTree.Branch("trueNuCCNC", trueNuCCNC, 'trueNuCCNC/I')
eventTree.Branch("nVertices", nVertices, 'nVertices/I')
eventTree.Branch("vtxX", vtxX, 'vtxX/F')
eventTree.Branch("vtxY", vtxY, 'vtxY/F')
eventTree.Branch("vtxZ", vtxZ, 'vtxZ/F')
eventTree.Branch("vtxIsFiducial", vtxIsFiducial, 'vtxIsFiducial/I')
eventTree.Branch("vtxDistToTrue", vtxDistToTrue, 'vtxDistToTrue/F')
eventTree.Branch("vtxBestComp", vtxBestComp, 'vtxBestComp/F')
eventTree.Branch("vtxScore", vtxScore, 'vtxScore/F')
eventTree.Branch("vtxFracHitsOnCosmic", vtxFracHitsOnCosmic, 'vtxFracHitsOnCosmic/F')
eventTree.Branch("nTracks", nTracks, 'nTracks/I')
eventTree.Branch("trackNHits", trackNHits, 'trackNHits[nTracks]/I')
eventTree.Branch("trackHitFrac", trackHitFrac, 'trackHitFrac[nTracks]/F')
eventTree.Branch("trackCharge", trackCharge, 'trackCharge[nTracks]/F')
eventTree.Branch("trackChargeFrac", trackChargeFrac, 'trackChargeFrac[nTracks]/F')
eventTree.Branch("trackRecoMuKE", trackRecoMuKE, 'trackRecoMuKE[nTracks]/F')
eventTree.Branch("trackRecoPrKE", trackRecoPrKE, 'trackRecoPrKE[nTracks]/F')
eventTree.Branch("trackClassified", trackClassified, 'trackClassified[nTracks]/I')
eventTree.Branch("trackRecoPID", trackRecoPID, 'trackRecoPID[nTracks]/I')
eventTree.Branch("trackRecoComp", trackRecoComp, 'trackRecoComp[nTracks]/F')
eventTree.Branch("nShowers", nShowers, 'nShowers/I')
eventTree.Branch("showerNHits", showerNHits, 'showerNHits[nShowers]/I')
eventTree.Branch("showerHitFrac", showerHitFrac, 'showerHitFrac[nShowers]/F')
eventTree.Branch("showerCharge", showerCharge, 'showerCharge[nShowers]/F')
eventTree.Branch("showerChargeFrac", showerChargeFrac, 'showerChargeFrac[nShowers]/F')
eventTree.Branch("showerRecoPl0E", showerRecoPl0E, 'showerRecoPl0E[nShowers]/F')
eventTree.Branch("showerRecoPl1E", showerRecoPl1E, 'showerRecoPl1E[nShowers]/F')
eventTree.Branch("showerRecoPl2E", showerRecoPl2E, 'showerRecoPl2E[nShowers]/F')
eventTree.Branch("showerClassified", showerClassified, 'showerClassified[nShowers]/I')
eventTree.Branch("showerRecoPID", showerRecoPID, 'showerRecoPID[nShowers]/I')
eventTree.Branch("showerRecoComp", showerRecoComp, 'showerRecoComp[nShowers]/F')


if args.isMC:
  totPOT_ = 0.
  totGoodPOT_ = 0.
  weights = Weights(args.weightfile)


#-------- begin file loop -----------------------------------------------------#
for filepair in files:

  ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
  ioll.add_in_filename(filepair[1])
  ioll.open()

  iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
  iolcv.add_in_file(filepair[1])
  iolcv.reverse_all_products()
  iolcv.initialize()

  kpsfile = rt.TFile(filepair[0])
  kpst = kpsfile.Get("KPSRecoManagerTree")

  try:
    nKPSTEntries = kpst.GetEntries()
  except:
    print("%s is empty. skipping..."%(filepair[0]))
    ioll.close()
    iolcv.finalize()
    kpsfile.Close()
    continue

  if args.isMC:
    potInFile, goodPotInFile = SumPOT(filepair[1])
    totPOT_ = totPOT_ + potInFile
    totGoodPOT_ = totGoodPOT_ + goodPotInFile

  #++++++ begin entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=
  for ientry in range(ioll.get_entries()):

    print("reached entry:", ientry)

    ioll.go_to(ientry)
    iolcv.read_entry(ientry)
    kpst.GetEntry(ientry)
  
    if kpst.run != ioll.run_id() or kpst.subrun != ioll.subrun_id() or kpst.event != ioll.event_id():
      print("EVENTS DON'T MATCH!!!")
      print("truth run/subrun/event: %i/%i/%i"%(ioll.run_id(),ioll.subrun_id(),ioll.event_id()))
      print("reco run/subrun/event: %i/%i/%i"%(kpst.run,kpst.subrun,kpst.event))
      continue


    if args.isMC:  

      mctruth = ioll.get_data(larlite.data.kMCTruth, "generator")
      nuInt = mctruth.at(0).GetNeutrino()
      lep = nuInt.Lepton()
      mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
      trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

      if not isFiducial(trueVtxPos):
        continue

      try:
        xsecWeight[0] = weights.get(kpst.run, kpst.subrun, kpst.event)
      except:
        print("Couldn't find weight for run %i, subrun %i, event %i in %s!!!"%(kpst.run, kpst.subrun, kpst.event, args.weightfile))
        continue

      if nuInt.CCNC() == 0:

        if lep.PdgCode() not in [11,13]:
          continue

        lepPDG = lep.PdgCode()

        if lepPDG == 13:
          mcleptons = ioll.get_data(larlite.data.kMCTrack, "mcreco")
        if lepPDG == 11:
          mcleptons = ioll.get_data(larlite.data.kMCShower, "mcreco")
        for mclepton in mcleptons:
          if mclepton.PdgCode() == lepPDG and mclepton.Process() == 'primary':
            mcLeptonUnCorr = mclepton
            break
  
        if not MCLeptonOkay(lep, mcLeptonUnCorr):
          print("Couldn't find MC lepton match!!!")
          continue

        trueMuContained[0] = -1
        if lepPDG == 13:
          if isInDetector(mcLeptonUnCorr.End()):
            trueMuContained[0] = 1
          else:
            trueMuContained[0] = 0

        trueLepPDG[0] = lepPDG
        trueLepE[0] = lep.Momentum().E()

        totLepPixI, lepTickLists, lepPixelDictList = getLeptonPixels(lepPDG, ioll, iolcv)

      else: #from "if nuInt.CCNC"
        trueLepPDG[0] = 0
        trueLepE[0] = -9.

      trueNuPDG[0] = nuInt.Nu().PdgCode()
      trueNuCCNC[0] = nuInt.CCNC()
      trueNuE[0] = nuInt.Nu().Momentum().E()
      

    else: #from "if args.isMC"
      trueLepPDG[0] = 0
      trueMuContained[0] = -1
      trueLepE[0] = -9.
      trueNuPDG[0] = 0
      trueNuCCNC[0] = -1
      trueNuE[0] = -9.

    fileid[0] = int(filepair[0][filepair[0].find("fileid")+6:filepair[0].find("fileid")+10])
    run[0] = kpst.run
    subrun[0] = kpst.subrun
    event[0] = kpst.event

    vertices = kpst.nuvetoed_v
    if args.oldVtxBranch:
      vertices = kpst.nufitted_v
  
    nVertices[0] = 0
    vtxScore[0] = -1.
    foundVertex = False
    for vtx in vertices:
      if vtx.keypoint_type != 0:
        continue
      nVertices[0] += 1
      foundVertex = True
      if vtx.netNuScore > vtxScore[0]:
        vtxScore[0] = vtx.netNuScore
        vertex = vtx

    if not foundVertex:
      vtxX[0] = -9999.
      vtxY[0] = -9999.
      vtxZ[0] = -9999.
      vtxIsFiducial[0] = -1
      vtxDistToTrue[0] = -99.
      vtxBestComp[0] = -1.
      vtxFracHitsOnCosmic[0] = -1.
      nTracks[0] = 0
      nShowers[0] = 0
      eventTree.Fill()
      continue

    if args.isMC:
      vtxDistToTrue[0] = getVertexDistance(trueVtxPos, vertex)
      if nuInt.CCNC() == 0:
        vtxBestComp[0] = getBestCompleteness(iolcv, vertex, lepPDG, totLepPixI, lepTickLists, lepPixelDictList)
      else:
        vtxBestComp[0] = -1.
    else:
      vtxDistToTrue[0] = -99.
      vtxBestComp[0] = -1.

    vtxX[0] = vertex.pos[0]
    vtxY[0] = vertex.pos[1]
    vtxZ[0] = vertex.pos[2]
    vtxTVec3 = rt.TVector3(vertex.pos[0], vertex.pos[1], vertex.pos[2])
    vtxIsFiducial[0] = int(isFiducial(vtxTVec3))

    nusel = larflow.reco.NuSelectionVariables()
    wcoverlapvars.analyze(vertex, nusel, iolcv)
    vtxFracHitsOnCosmic[0] = nusel.frac_allhits_on_cosmic

    nTracks[0] = vertex.track_v.size()
    nShowers[0] = vertex.shower_v.size()

    evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
    csmImage2D = iolcv.get_data(larcv.kProductImage2D, "thrumu")
    adc_v = evtImage2D.Image2DArray()
    thrumu_v = csmImage2D.Image2DArray()

    vertexPixels = []
    vertexCharge = 0.
    vertexNHits = 0


    for iTrk, trackCls in enumerate(vertex.track_hitcluster_v):

      vertexNHits += trackCls.size()
      trackNHits[iTrk] = trackCls.size()
      trackCharge[iTrk], vertexPixels, vertexCharge = addClusterCharge(iolcv,trackCls,vertexPixels,vertexCharge,10.)
      trackRecoMuKE[iTrk] = vertex.track_kemu_v[iTrk]
      trackRecoPrKE[iTrk] = vertex.track_keproton_v[iTrk]

      cropPt = vertex.track_v[iTrk].End()
      prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,trackCls,cropPt,10.,512,512)
      skip = False
      for p in range(3):
        if prong_vv[p].size() < 10:
          skip = True
          break
      if skip:
        trackClassified[iTrk] = 0
        trackRecoPID[iTrk] = 0
        trackRecoComp[iTrk] = -1.
        continue

      prongImage = makeImage(prong_vv).to(args.device)
      prongCNN_out = model(prongImage)
      trackClassified[iTrk] = 1
      trackRecoPID[iTrk] = getPID(prongCNN_out[0].argmax(1).item())
      trackRecoComp[iTrk] = prongCNN_out[1].item()


    for iShw, shower in enumerate(vertex.shower_v):

      vertexNHits += shower.size()
      showerNHits[iShw] = shower.size()
      showerCharge[iShw], vertexPixels, vertexCharge = addClusterCharge(iolcv,shower,vertexPixels,vertexCharge, 10.)
      showerRecoPl0E[iShw] = vertex.shower_plane_mom_vv[iShw][0].E()
      showerRecoPl1E[iShw] = vertex.shower_plane_mom_vv[iShw][1].E()
      showerRecoPl2E[iShw] = vertex.shower_plane_mom_vv[iShw][2].E()

      cropPt = vertex.shower_trunk_v[iShw].Vertex()
      prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,shower,cropPt,10.,512,512)
      skip = False
      for p in range(3):
        if prong_vv[p].size() < 10:
          skip = True
          break
      if skip:
        showerClassified[iShw] = 0
        showerRecoPID[iShw] = 0
        showerRecoComp[iShw] = -1.
        continue

      prongImage = makeImage(prong_vv).to(args.device)
      prongCNN_out = model(prongImage)
      showerClassified[iShw] = 1
      showerRecoPID[iShw] = getPID(prongCNN_out[0].argmax(1).item())
      showerRecoComp[iShw] = prongCNN_out[1].item()

    for i in range(nTracks[0]):
      trackHitFrac[i] = trackNHits[i] / (1.0*vertexNHits)
      trackChargeFrac[i] = trackCharge[i] / vertexCharge

    for i in range(nShowers[0]):
      showerHitFrac[i] = showerNHits[i] / (1.0*vertexNHits)
      showerChargeFrac[i] = showerCharge[i] / vertexCharge

    eventTree.Fill()

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  iolcv.finalize()
  kpsfile.Close()

#-------- end file loop -----------------------------------------------------#


if args.isMC:
  totPOT[0] = totPOT_
  totGoodPOT[0] = totGoodPOT_
  potTree.Fill()

outRootFile.cd()
eventTree.Write("",rt.TObject.kOverwrite)
if args.isMC:
  potTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()

