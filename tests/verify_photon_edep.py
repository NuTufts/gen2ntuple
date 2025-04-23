#!/bin/env python3
from __future__ import print_function
import os,sys,argparse
from math import fabs

parser = argparse.ArgumentParser("Test MCPixelPGraph")
parser.add_argument("-idl", "--input-dlmerged",required=True,type=str,help="Input DL Merged (larcv+larlite) input file")
parser.add_argument("-int", "--input-ntuple",required=True,type=str,help="Input Gen2 NTuple")
parser.add_argument("-adc", "--adc",type=str,default="wire",help="Name of tree with Wire ADC values [default: wire]")
#parser.add_argument("-tb",  "--tick-backward",action='store_true',default=False,help="Input LArCV data is tick-backward [default: false]")
parser.add_argument("-d",   "--debug", action='store_true', default=False, help="Run in debug mode")
#parser.add_argument('-n',   "--nentries", required=False, type=int, default=-1, help="number of entries to run")
#parser.add_argument('-e',   "--entry", required=False, type=int, default=-1, help="start at given entry number")
parser.add_argument('-v',   "--vis", required=False, default=False, action='store_true', help="if flag provided, will visualize event")
args = parser.parse_args()

import ROOT as rt
rt.gROOT.SetBatch(True)

from larcv import larcv
from larlite import larlite
from ublarcvapp import ublarcvapp
from larflow import larflow

"""
test script that demos the MCPixelPGraph class.
"""

rt.gStyle.SetOptStat(0)

# ================================================
# SETUP INPUTS
# ================================================
ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.set_data_to_read( "mctrack",  "mcreco" ) # larlite.data.kMCTrack
ioll.set_data_to_read( "mcshower", "mcreco" ) # larlite.data.kMCShower
ioll.set_data_to_read( "mctruth", "generator" ) # larlite.data.kMCTruth
ioll.add_in_filename(  args.input_dlmerged )
ioll.open()

# MCC9 DL Gen Input LArCV Files are Tick-Backward
iolcv = larcv.IOManager( larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward )
iolcv.add_in_file( args.input_dlmerged )
iolcv.reverse_all_products()
iolcv.initialize()

tot_nentries = ioll.get_entries()
start_entry = 0
print("Number of entries: ",tot_nentries)
#if args.entry > 0:
#    start_entry = args.entry
#if args.nentries>0:
#    nentries = args.nentries
#else:
nentries = tot_nentries
end_entry = start_entry + nentries
if end_entry>tot_nentries:
    end_entry = tot_nentries

input_ntuple = rt.TFile(args.input_ntuple,"READ")
tree = input_ntuple.Get("EventTree")
nentries_ntuple = tree.GetEntries()

if nentries_ntuple!=tot_nentries:
    print("Number of entries do not match. dlmerged=",tot_nentries," vs. ntuple=",nentries_ntuple,". Stopping now.")
    sys.exit(0)

print("Start loop.")

# ==================================
# Classes we use to do analysis
# ==================================
mcpg = ublarcvapp.mctools.MCPixelPGraph()
mcpg.set_adc_treename( args.adc )
if args.debug:
    mcpg.set_verbosity( "debug" )
else:
    mcpg.set_verbosity( "info" )

photonTruthMetrics = larflow.reco.ShowerTruthMetricsMaker()    

tmp = rt.TFile("output_verify_photon_edep.root","recreate")
c = rt.TCanvas("c","c",2400,600)
c.Divide(3,1)

for ientry in range( start_entry, end_entry ):

    print() 
    print("==========================")
    print("===[ EVENT ",ientry," ]===")
    ioll.go_to(ientry)
    iolcv.read_entry(ientry)
    tree.GetEntry(ientry)

    # use mctrack, mcshower, (genie) generator to build truth particle graph
    mcpg.clear()    
    mcpg.set_cluster_neutrino_particles(True)    
    mcpg.buildgraphonly( ioll )

    # we are interested in checking gamma info, so we make sure we have some
    ## isolate events
    has_photons = False
    for inode in range(mcpg.node_v.size()):
        node = mcpg.node_v.at(inode)
        if node.pid in [22]:
            has_photons = True
        if has_photons:
            break
        
    if not has_photons:
        continue

    mcpg.clear()

    print("images needed for pixel matching to particles")
    ev_adc = iolcv.get_data( "image2d", args.adc )
    ev_instance = iolcv.get_data( "image2d", "instance" )
    ev_ancestor = iolcv.get_data( "image2d", "ancestor" )
    ev_segment  = iolcv.get_data( "image2d", "segment" )
    ev_larflow  = iolcv.get_data( "image2d", "larflow" )
    print("number of adc images: ",ev_adc.Image2DArray().size())
    print("number of larflow images: ",ev_larflow.Image2DArray().size())
    adc_v = ev_adc.Image2DArray()
    for p in range(adc_v.size()):
        print(" image[",p,"] ",adc_v[p].meta().dump())

    # rebuild the graph, but now using images with true information about the simulation
    # this information can tell us
    # (1) a better first energy deposit location for photons
    # (2) a better measure of how much energy has been deposited in the TPC
    mcpg.buildgraph( iolcv, ioll )
    mcpg.printGraph(0,False)    

    # look for photons and extract info from simulation truth
    # using MCPixelPGraph and same helper class larflow.reco.ShowerTruth
    mcpg_photon_starts = {}
    mcpg_pixelsums = {}
    mcpg_photon_imgpos = {}
    for i in range(mcpg.node_v.size()):
        node = mcpg.node_v.at(i)
        if node.pid!=22:
            continue;
        print("running fixingPhotonStartPoints(...) for photons")
        gamma_start = mcpg.getParticleEDepPos( node.tid )
        mcpg_photon_starts[node.tid] = gamma_start
        mcpg_pixelsums[node.tid] = mcpg.getTruePhotonTrunkPlanePixelSums( node.tid )
        if mcpg_pixelsums[node.tid].size()!=3:
            mcpg_pixelsums[node.tid] = (0,0,0)
        mcpg_photon_imgpos[node.tid] = [ node.imgpos4_edep[v] for v in range(node.imgpos4_edep.size()) ]

    # look through our ntuple entry
    ntuple_photon_starts = {}
    ntuple_photon_ends = {}
    ntuple_photon_pixsum = {}
    for ids in range( tree.nTrueSimParts ):
        if tree.trueSimPartPDG[ids]!=22:
            continue

        tid = tree.trueSimPartTID[ids]        
        if tid not in mcpg_photon_starts:
            continue

        match_edep   = True
        ntuple_photon_ends[tid] = [0.0,0.0,0.0]
        for v,pos in enumerate([tree.trueSimPartEDepX[ids], tree.trueSimPartEDepY[ids], tree.trueSimPartEDepZ[ids] ]):
            edep_diff = fabs( pos-mcpg_photon_starts[ tid ][v] )
            if edep_diff>1.0e-5:
                match_edep = False
        for v,pos in enumerate([tree.trueSimPartEndX[ids],tree.trueSimPartEndY[ids], tree.trueSimPartEndZ[ids]]):
            ntuple_photon_ends[tid][v] = pos
        match_pixsum = True
        for v,pixsum in enumerate([tree.trueSimPartPixelSumUplane[ids],tree.trueSimPartPixelSumVplane[ids],tree.trueSimPartPixelSumYplane[ids]]):
            pixsum_diff = fabs( pixsum-mcpg_pixelsums[ tree.trueSimPartTID[ids] ][v] )
            if pixsum_diff>1.0e-5:
                match_pixsum = False

        print("Photon[ trackid=",tree.trueSimPartTID[ids]," ] match_edep=",match_edep," match_pixsum=",match_pixsum)
              
            
    
    if args.vis:
        # We're going to make a canvas that we'll save, visualizing the photons
        
        #print("====================================================")
        #print("CONSTRUCTED PARTICLE GRAPH [WITH NU VERTEX GROUPING]")
        #mcpg.printGraph(0,False)
        #print("====================================================")
        # make histogram
        marker_v = []
        #hist_v = mcpg_nu.makeTH2D( "hentry%d"%(ientry) )
        hist_v = larcv.rootutils.as_th2d_v( ev_adc.as_vector(), "hentry%d"%(ientry) )
        segment_v = larcv.rootutils.as_th2d_v( ev_segment.as_vector(), "hsegment%d"%(ientry))
        instance_v = larcv.rootutils.as_th2d_v( ev_instance.as_vector(), "hinstance%d"%(ientry))
        ancestor_v = larcv.rootutils.as_th2d_v( ev_ancestor.as_vector(), "hancestor%d"%(ientry))
        
        tmp.cd()
        c.cd()
        c.Clear()
        c.Divide(3,1)
        c.Draw()
        c.Update()
        # we zoom around the photon start position
        iphoton = 0
        for tid in mcpg_photon_starts:
            text_v = []
            marker_v = []
            
            # end position
            if tid in ntuple_photon_ends:
                trunk_end = ntuple_photon_ends[tid]
                end_imgpos = ublarcvapp.mctools.MCPos2ImageUtils.Get().to_imagepos( trunk_end[0], trunk_end[1], trunk_end[2], 0.0 )
            else:
                print("WEIRD: tid=",tid," not in ntuple_photon_ends.")
                print("ntuple_photon_ends keys: ")
                for k,v in ntuple_photon_ends.items():
                    print("  ",k," : ",v)
                end_imgpos = None
            
            for p in range(3):
                c.cd(p+1)
                h = hist_v[p]
                wire = mcpg_photon_imgpos[tid][p]
                tick = mcpg_photon_imgpos[tid][3]
                h.Draw("colz")
                t = rt.TText( wire, tick+50*6, "trackid=%d pixsum=%.2f"%(tid,mcpg_pixelsums[tid][p]))
                t.SetNDC(False)                
                t.Draw()
                text_v.append(t)

                #x = node.imgpos4_edep[ih]
                #y = node.imgpos4_edep[3]
                #print("nodeidx=",node.nodeidx," pid=",node.pid," (x,y)=",(x,y))
                
                m = rt.TMarker( wire, tick, 4 )
                m.SetNDC(False)
                m.SetMarkerColor(rt.kMagenta)
                m.Draw()
                marker_v.append(m)

                if end_imgpos is not None:
                    m2 = rt.TMarker( end_imgpos[p], end_imgpos[3], 4 )
                    m2.SetNDC(False)
                    m2.SetMarkerColor(rt.kCyan)
                    m2.Draw()
                    marker_v.append(m2)
                
                h.GetZaxis().SetRangeUser(0,300)
                h.GetXaxis().SetRangeUser( wire-100, wire+101 )
                h.GetYaxis().SetRangeUser( tick-100*6, tick+101*6 )
                
                c.Update()
            tmp.cd()
            c.Write("c_event%d_tid%d_plane%d"%(ientry,tid,p))
            c.SaveAs("c_event%d_tid%d_plane%d.png"%(ientry,tid,p))
            print("save canvas to root file: c_event%d_tid%d_plane%d"%(ientry,tid,p))
    # end of if visualize photons

print("Save output")
tmp.Write()
tmp.Close()
print("=== FIN ==")
