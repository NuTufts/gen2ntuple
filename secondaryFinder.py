
import sys, argparse
import numpy as np
import ROOT as rt

from helpers.larflowreco_ana_funcs import getCosThetaGravVector


parser = argparse.ArgumentParser("Make energy histograms from a bnb nu overlay ntuple file")
parser.add_argument("-i", "--infile", type=str, required=True, help="input ntuple file")
parser.add_argument("-o", "--outfile", type=str, default="example_ntuple_analysis_script_output.root", help="output root file name")
parser.add_argument("-fc", "--fullyContained", action="store_true", help="only consider fully contained events")
parser.add_argument("-ncc", "--noCosmicCuts", action="store_true", help="don't apply cosmic rejection cuts")
args = parser.parse_args()

#needed for proper scaling of error bars:
rt.TH1.SetDefaultSumw2(rt.kTRUE)

#open input file and get event and POT trees
ntuple_file = rt.TFile(args.infile)
eventTree = ntuple_file.Get("EventTree")
potTree = ntuple_file.Get("potTree")

#we will scale histograms to expected event counts from POT in runs 1-3: 6.67e+20
targetPOT = 6.67e+20
targetPOTstring = "6.67e+20" #for plot axis titles

#calculate POT represented by full ntuple file after applying cross section weights
ntuplePOTsum = 0.
for i in range(potTree.GetEntries()):
  potTree.GetEntry(i)
  ntuplePOTsum = ntuplePOTsum + potTree.totGoodPOT

#define histograms to fill
#we will write histograms to output file for:
neutralPhotonHist = rt.TH1F("neutralPhotonHist", "Energy of NCC events with secondary photons",30,0,6)



#set histogram axis titles and increase line width
def configureHist(h):
  h.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h.GetXaxis().SetTitle("energy (GeV)")
  h.SetLineWidth(2)
  return h
neutralPhotonHist = configureHist(neutralPhotonHist)
#h_trueMuE = configureHist(h_trueMuE)
#h_recoMuE = configureHist(h_recoMuE)

counter = 0
#begin loop over events in ntuple file

for i in range(eventTree.GetEntries()):

  eventTree.GetEntry(i)

#Determine whether there are photons as secondary particles
  photonInSecondary = False
  primList = []
  for x in range(len(eventTree.trueSimPartPDG)):
    for x in range(len(eventTree.trueSimPartTID)):
      if eventTree.trueSimPartTID[x] == eventTree.trueSimPartMID[x]:
        primList.append(eventTree.trueSimPartTID[x])
    if eventTree.trueSimPartPDG[x] == 22:
      if eventTree.trueSimPartMID[x] in primList:
        photonInSecondary = True
  if photonInSecondary == False:
    continue


  #filter for neutral current
  if eventTree.trueNuCCNC != 1:
    continue

  
  #Filter out cosmic cuts, because that seems like the right thing to do
 # if not args.noCosmicCuts:
  #skip events where all hits overlap with tagged cosmic rays
  #  if eventTree.vtxFracHitsOnCosmic >= 1.:
   #   continue



  #fill histograms, weighting each event by cross-section weight to allow for proper POT scaling

  #fill our histogram using energy
  neutralPhotonHist.Fill(eventTree.trueNuE, eventTree.xsecWeight)

#----- end of event loop ---------------------------------------------#

print("there are", counter, "secondary photons")
#scale histograms to target POT
neutralPhotonHist.Scale(targetPOT/ntuplePOTsum)

#create output root file and write histograms to file
outFile = rt.TFile(args.outfile, "RECREATE")
neutralPhotonHist.Write()
