
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow
from math import sqrt as sqrt
from math import acos as acos
from math import pi
import ctypes
import os

def get_kpsana_file_samhash( filename ):
  # looks like: larflowreco_fileid0000_d9679e9b-3be3-4411-bc25-6e2cea860827_kpsrecomanagerana.root
  bname = os.path.basename( filename.strip() )
  samhash = bname.split("_")[-2]
  return samhash

def parse_merged_dl_filelist( mdlfiles ):
  # we return a dict from hash to file path
  # looks like: mcc9_v29e_dl_run3_G1_extbnb_dlana/merged_dlana_d9679e9b-3be3-4411-bc25-6e2cea860827.root
  filedict = {}
  with open(mdlfiles,'r') as f:
    ll = f.readlines()
    for l in ll:
      l = l.strip()
      fbase = os.path.basename(l.strip())
      tag = fbase.split("_")[-1].replace(".root","")
      filedict[tag] = l
  print("parse dlmerged input file list. number of files in samtag -> filepath dict: ",len(filedict))
  return filedict

def getFiles(mdlTag, kpsfiles, mdlfiles):
  files = []

  tagdict = parse_merged_dl_filelist( mdlfiles )
  
  for kpsfile in kpsfiles:
    # dlrecofilelist = open(mdlfiles, "r")
    # for line in dlrecofilelist:
    #   dlrecofile = line.replace("\n","")
    #   samtag = dlrecofile[dlrecofile.find(mdlTag):].replace(mdlTag,"").replace(".root","")
    #   if samtag in kpsfile:
    #     files.append([kpsfile, dlrecofile])
    #     break
    # dlrecofilelist.close()
    recotag = get_kpsana_file_samhash( kpsfile )
    if recotag in tagdict:
      files.append( [ kpsfile.strip(), tagdict[recotag].strip() ] )
    else:
      print("DID NOT FIND KPS Reco File sam-tag (",recotag,") in source file tag dict (len=",len(tagdict))
  return files

detCrds = [[0., 256.35], [-116.5, 116.5], [0, 1036.8]]
fidCrds = [ [detCrds[0][0] + 20. , detCrds[0][1] - 20.] ]
fidCrds.append( [detCrds[1][0] + 20. , detCrds[1][1] - 20.] )
fidCrds.append( [detCrds[2][0] + 20. , detCrds[2][1] - 60.] )
fidCrdsBig = [ [detCrds[0][0] + 10. , detCrds[0][1] - 10.] ]
fidCrdsBig.append( [detCrds[1][0] + 10. , detCrds[1][1] - 10.] )
fidCrdsBig.append( [detCrds[2][0] + 10. , detCrds[2][1] - 30.] )

def inRange(x, bnd):
  return (x >= bnd[0] and x <= bnd[1])

def isFiducial(p):
  return (inRange(p.X(),fidCrds[0]) and inRange(p.Y(),fidCrds[1]) and inRange(p.Z(),fidCrds[2]))

def isFiducialBig(p):
  return (inRange(p.X(),fidCrdsBig[0]) and inRange(p.Y(),fidCrdsBig[1]) and inRange(p.Z(),fidCrdsBig[2]))

def isInDetector(p):
  return (inRange(p.X(),detCrds[0]) and inRange(p.Y(),detCrds[1]) and inRange(p.Z(),detCrds[2]))

class WCFiducial(ctypes.Structure):
  pass
libpath = os.path.dirname(os.path.realpath(__file__))
libwc = ctypes.cdll.LoadLibrary("%s/lib_wirecell_fiducial_volume.so"%libpath)
libwc.WCFiducial_new.argtypes = ()
libwc.WCFiducial_new.restype = ctypes.POINTER(WCFiducial)
libwc.WCFiducial_insideFV.argtypes = ctypes.POINTER(WCFiducial), ctypes.c_double, ctypes.c_double, ctypes.c_double
libwc.WCFiducial_insideFV.restype = ctypes.c_bool
WCFiducialClass = libwc.WCFiducial_new()

def isFiducialWCSCE(p):
  return libwc.WCFiducial_insideFV(WCFiducialClass, p.X(), p.Y(), p.Z())

def getVertexDistance(pos3v, recoVtx):
  xdiffSq = (pos3v.X() - recoVtx.pos[0])**2
  ydiffSq = (pos3v.Y() - recoVtx.pos[1])**2
  zdiffSq = (pos3v.Z() - recoVtx.pos[2])**2
  return sqrt( xdiffSq + ydiffSq + zdiffSq )

def getTrackLength(track):
  tracklen = 0
  for i in range(track.NumberTrajectoryPoints()-1):
    step = track.LocationAtPoint(i)
    nextStep = track.LocationAtPoint(i+1)
    tracklen += (nextStep-step).Mag();
  return tracklen;

def getDistance(a, b):
  return sqrt( (a.X() - b.X())**2 + (a.Y() - b.Y())**2 + (a.Z() - b.Z())**2)

def getDirection(a, b):
  dirX = b.X() - a.X()
  dirY = b.Y() - a.Y()
  dirZ = b.Z() - a.Z()
  mag = sqrt(dirX**2 + dirY**2 + dirZ**2)
  return (dirX/mag, dirY/mag, dirZ/mag)

def getCosTVecAngle(a, b):
  aMag = sqrt(a.X()**2 + a.Y()**2 + a.Z()**2)
  bMag = sqrt(b.X()**2 + b.Y()**2 + b.Z()**2)
  return (a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z())/(aMag*bMag)

def getTVecAngle(a, b):
  return acos(getCosTVecAngle(a,b))

def getCosThetaBeamTrack(track):
  beamDir = rt.TVector3(0, 0, 1)
  recoStart = track.Vertex()
  for i in range(track.NumberTrajectoryPoints()):
    anglePoint = track.LocationAtPoint(i)
    if getDistance(recoStart, track.LocationAtPoint(i)) > 5.:
      break
  trackDir = rt.TVector3(anglePoint.X() - recoStart.X(),
                         anglePoint.Y() - recoStart.Y(),
                         anglePoint.Z() - recoStart.Z())
  return getCosTVecAngle(trackDir, beamDir)

def getCosThetaGravTrack(track):
  gravDir = rt.TVector3(0, -1, 0)
  recoStart = track.Vertex()
  for i in range(track.NumberTrajectoryPoints()):
    anglePoint = track.LocationAtPoint(i)
    if getDistance(recoStart, track.LocationAtPoint(i)) > 5.:
      break
  trackDir = rt.TVector3(anglePoint.X() - recoStart.X(),
                         anglePoint.Y() - recoStart.Y(),
                         anglePoint.Z() - recoStart.Z())
  return getCosTVecAngle(trackDir, gravDir)

def getCosThetaBeamShower(showerTrunk):
  beamDir = rt.TVector3(0, 0, 1)
  showerDir = rt.TVector3(showerTrunk.VertexDirection().X(),
                          showerTrunk.VertexDirection().Y(),
                          showerTrunk.VertexDirection().Z())
  return getCosTVecAngle(showerDir, beamDir)

def getCosThetaGravShower(showerTrunk):
  gravDir = rt.TVector3(0, -1, 0)
  showerDir = rt.TVector3(showerTrunk.VertexDirection().X(),
                          showerTrunk.VertexDirection().Y(),
                          showerTrunk.VertexDirection().Z())
  return getCosTVecAngle(showerDir, gravDir)

def getCosThetaBeamVector(x, y, z): 
  gravDir = rt.TVector3(0, 0, 1)
  vecDir = rt.TVector3(x, y, z)
  return getCosTVecAngle(vecDir, gravDir)

def getCosThetaGravVector(x, y, z): 
  gravDir = rt.TVector3(0, -1, 0)
  vecDir = rt.TVector3(x, y, z)
  return getCosTVecAngle(vecDir, gravDir)

def getSCECorrectedPos(point, sce):
  offset = sce.GetPosOffsets(point.X(), point.Y(), point.Z())
  corrected = rt.TVector3(point.X() - offset[0] + 0.7, point.Y() + offset[1], point.Z() + offset[2])
  return corrected

