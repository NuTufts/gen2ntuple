
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow
from math import sqrt as sqrt
from math import acos as acos
from math import pi


def getFiles(mdlTag, kpsfiles, mdlfiles):
  files = []
  for kpsfile in kpsfiles:
    dlrecofilelist = open(mdlfiles, "r")
    for line in dlrecofilelist:
      dlrecofile = line.replace("\n","")
      samtag = dlrecofile[dlrecofile.find(mdlTag):].replace(mdlTag,"").replace(".root","")
      if samtag in kpsfile:
        files.append([kpsfile, dlrecofile])
        break
    dlrecofilelist.close()
  return files

def MCLeptonOkay(nuIntLep, mcObjLep):
  if mcObjLep.PdgCode() != nuIntLep.PdgCode():
    return False
  nuIntStart = nuIntLep.Position()
  mcObjStart = mcObjLep.Start()
  xdiffSq = (nuIntStart.X() - mcObjStart.X())**2
  ydiffSq = (nuIntStart.Y() - mcObjStart.Y())**2
  zdiffSq = (nuIntStart.Z() - mcObjStart.Z())**2
  if sqrt(xdiffSq + ydiffSq + zdiffSq) > 1.:
    return False
  return True

detCrds = [[0., 256.35], [-116.5, 116.5], [0, 1036.8]]
fidCrds = [ [detCrds[0][0] + 20. , detCrds[0][1] - 20.] ]
fidCrds.append( [detCrds[1][0] + 20. , detCrds[1][1] - 20.] )
fidCrds.append( [detCrds[2][0] + 20. , detCrds[2][1] - 60.] )

def inRange(x, bnd):
  return (x >= bnd[0] and x <= bnd[1])

def isFiducial(p):
  return (inRange(p.X(),fidCrds[0]) and inRange(p.Y(),fidCrds[1]) and inRange(p.Z(),fidCrds[2]))

def isInDetector(p):
  return (inRange(p.X(),detCrds[0]) and inRange(p.Y(),detCrds[1]) and inRange(p.Z(),detCrds[2]))

def getVertexDistance(pos3v, recoVtx):
  xdiffSq = (pos3v.X() - recoVtx.pos[0])**2
  ydiffSq = (pos3v.Y() - recoVtx.pos[1])**2
  zdiffSq = (pos3v.Z() - recoVtx.pos[2])**2
  return sqrt( xdiffSq + ydiffSq + zdiffSq )

def getLeptonPixels(pdg, ioll, iolcv):
  # get lepton node
  mcpg = ublarcvapp.mctools.MCPixelPGraph()
  mcpg.set_adc_treename("wire")
  mcpg.buildgraph(iolcv, ioll)
  foundLepNode = False
  for node in mcpg.node_v:
    if node.pid == pdg and node.tid == node.mtid:
      lepNode = node
      foundLepNode = True
      break
  if not foundLepNode:
    sys.exit("Couldn't find lepton node in MCPixelPGraph!!")
  # get pixel value array
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()
  # calculate pixel info return variables
  totLepPixI = 0.
  lepTickLists = [ [], [], [] ]
  lepPixelDictList = [ {}, {}, {} ]
  for p in range(3):
    lepPix = lepNode.pix_vv[p]
    for iP in range(lepPix.size()//2):
      simTick = lepPix[2*iP]
      simWire = lepPix[2*iP+1]
      row = (simTick - 2400)//6
      totLepPixI = totLepPixI + image2Dvec[p].pixel(row,simWire)
      if simTick not in lepTickLists[p]:
        lepTickLists[p].append(simTick)
        lepPixelDictList[p][simTick] = [simWire]
      else:
        lepPixelDictList[p][simTick].append(simWire)
  return totLepPixI, lepTickLists, lepPixelDictList

def getBestCompleteness(iolcv, vertex, lepPDG, totLepPixI, lepTickLists, lepPixelDictList):
  bestComp = -1.
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()
  if lepPDG == 13:
    clusterVec = vertex.track_hitcluster_v
  if lepPDG == 11:
    clusterVec = vertex.shower_v
  for cluster in clusterVec:
    matchedPixels = []
    matchedUniquePixI = 0.
    totPixI = 0.
    for hit in cluster:
      for p in range(3):
        row = (hit.tick - 2400)//6
        pixVal = image2Dvec[p].pixel(row, hit.targetwire[p])
        totPixI = totPixI + pixVal
        if hit.tick in lepTickLists[p]:
          if hit.targetwire[p] in lepPixelDictList[p][hit.tick]:
            matchedPixel = [ p, hit.tick, hit.targetwire[p] ]
            if matchedPixel not in matchedPixels:
              matchedUniquePixI = matchedUniquePixI + pixVal
              matchedPixels.append(matchedPixel)
    comp = matchedUniquePixI/totLepPixI
    if comp > bestComp:
      bestComp = comp
  return bestComp


def getDistance(a, b):
  return sqrt( (a.X() - b.X())**2 + (a.Y() - b.Y())**2 + (a.Z() - b.Z())**2)

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

def getCosThetaBeamShower(showerTrunk):
  beamDir = rt.TVector3(0, 0, 1)
  showerDir = rt.TVector3(showerTrunk.VertexDirection().X(),
                          showerTrunk.VertexDirection().Y(),
                          showerTrunk.VertexDirection().Z())
  return getCosTVecAngle(showerDir, beamDir)



