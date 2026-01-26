import os,sys
import ROOT as rt
from larcv import larcv
from larflow import larflow
#import hdf5 as h5

# larpid code
# ./../../photon_analysis/prongCNN/models/

from filelists import get_filelist
from larpid_interface import run_larpid

def get_selected_event_info_fromrootfile( selection_rootfile ):
    """
******************************************************************************
*Tree    :analysis_tree: Processed Events                                       *
******************************************************************************
*Br    0 :eventselect_run : eventselect_run/I                                *
*............................................................................*
*Br    1 :eventselect_subrun : eventselect_subrun/I                          *
*............................................................................*
*Br    2 :eventselect_event : eventselect_event/I                            *
*............................................................................*
*Br    3 :eventselect_fileid : eventselect_fileid/I                          *
*............................................................................*
*Br    4 :eventselect_vtxIdx : eventselect_vtxIdx/I                          *
*............................................................................*
*Br    5 :eventselect_vertexScore : eventselect_vertexScore/F                *
*............................................................................*
*Br    6 :eventselect_recoNuE : eventselect_recoNuE/F                        *
*............................................................................*
*Br    7 :eventselect_selectevent_passes : eventselect_selectevent_passes/I  *
*............................................................................*
*Br    8 :eventselect_vertex : eventselect_vertex[3]/F                       *
*............................................................................*
"""
    rinput = rt.TFile( selection_rootfile )
    tree = rinput.Get("analysis_tree")

    event_list = []
    for ientry in range(tree.GetEntries()):
        tree.GetEntry(ientry)
        event_info = {'run':tree.eventselect_run,
                      'subrun':tree.eventselect_subrun,
                      'event':tree.eventselect_event,
                      'fileid':tree.eventselect_fileid,
                      'vtxIdx':tree.eventselect_vtxIdx}
        event_list.append( event_info )
        
    print(f"Number of events in the selection root file: {len(event_list)}")
    return event_list

def extract_spacepoints( nuvtx ):
    ntracks  = nuvtx.track_v.size()
    iprong = 0
    prong_spacepoints = []
    for iTrk in range(ntracks):
        trackCls = nuvtx.track_hitcluster_v.at(iTrk)
        nhits = trackCls.size()
        hitarray = np.zeros( (nhits,5) )
        for ihit in range(trackCls.size()):
            hit = trackCls.at(ihit)
            hitarray[ihit,3] = iprong
            hitarray[ihit,4] = 0
            for i in range(3):
                hitarray[ihit,i] = hit[i]
        prong_spacepoints.append(hitarray)
        iprong += 1
        
    for iShw in range(nshowers):
        shower = nuvtx.shower_v.at(iShw)
        nhits  = shower.size()
        hitarray = np.zeros( (nhits,5) )
        for ihit in range(shower.size()):
            hit = shower.at(ihit)
            hitarray[ihit,3] = iprong
            hitarray[ihit,4] = 1
            for i in range(3):
                hitarray[ihit,i] = hit[i]
        prong_spacepoints.append(hitarray)
        iprong += 1

    spacepoints = np.concatenate( prong_spacepoints, axis=0 )
    return spacepoints
    

def extract_event_info( event_info, iolcv_path, iolarlite_path, kpsreco_path ):
    """
    We extract info for visualization/study.
    1. image crop around vertex: pixel values
    2. image crop around vertex: ssnet
    3. 
    """
    print("=== EXTRACT EVENT ====")
    print(f"iolcv_path: {iolcv_path}")
    print(f"iolarlite_path: {iolarlite_path}")
    print(f"kps reco path: {kpsreco_path}")
    print("event info: ",event_info)

    run = event_info['run']
    subrun = event_info['subrun']
    event = event_info['event']
    vtxIdx = event_info['vtxIdx']

    model = None

    iolcv = larcv.IOManager(larcv.IOManager.kREAD,'larcv',larcv.IOManager.kTickBackward)
    iolcv.add_in_file( iolcv_path )
    iolcv.reverse_all_products()
    iolcv.initialize()

    # Get the nu reco object.
    rfile_reco = rt.TFile( kpsreco_path, 'read' )
    recotree = rfile_reco.Get("KPSRecoManagerTree")
    nentries = recotree.GetEntries()
    found  = False
    larpid_output = None
    for ientry in range(nentries):
        recotree.GetEntry(ientry)
        if run!=recotree.run or subrun!=recotree.subrun or event!=recotree.event:
            continue
        iolcv.read_entry(ientry)
        
        nuvertices = recotree.nuvetoed_v.size()
        if nuvertices<=vtxIdx:
            raise ValueError(f"Number of vertices in this event ({nvertices}) is less than target vertex index ({vtxIdx})")
        nuvtx = recotree.nuvetoed_v.at(vtxIdx)
        found = True
        print(f"Found vertex index: {vtxIdx}")

        vtx_imgcol  = [ nuvtx.col_v.at(i) for i in range(3) ]
        vtx_imgrow  = nuvtx.row
        vtx_imgtick = nuvtx.tick

        larpid_output = run_larpid( nuvtx, iolcv, model )
        for k,vdict in larpid_output.items():
            if 'larpid_img' in vdict:
                print(k,": ",vdict['larpid_img'].shape)

        # extract 3D points
        prong_spacepoints = extract_spacepoints( nuvtx )

        # save as root histograms
        rootfile = f"selected_event_{fileid}_{run}_{subrun}_{event}_{vtxIdx}.root"
        outroot = rt.TFile(rootfile,'recreate')
        for k,vdict in larpid_output.items():
            prongtype,prongindex = k
            if 'larpid_img' in vdict:
                prongimages = vdict['larpid_img']
                for ii in range(prongimages.shape[0]):
                    #print(k,": ",vdict['larpid_img'].shape)
                    h2d = rt.TH2D(f"h{prongtype}_{prongindex}_{ii}","",512,0,512,512,0,512)
                    for ix in range(512):
                        for iy in range(512):
                            h2d.SetBinContent(ix+1,iy+1,prongimages[ii,iy,ix])
                    h2d.Write()
        outroot.Close()
        
        
        if True:
            break
        
    rfile_reco.Close()
    return found

if __name__=="__main__":
    
    event_info_list = get_selected_event_info_fromrootfile( "mcc9_v28_wctagger_run3_bnb1e19_eventselection_20260124_134343.root" )

    fileid_dict = get_filelist( "run3_opendata_1e19" )

    for event_info in event_info_list:
        fileid = event_info['fileid']
        fileid_files = fileid_dict[fileid]
        if "reco" not in fileid_files:
            continue

        dlmerged = fileid_files['dlmerged']
        kpsreco_path = fileid_files['reco']
        iolcv_path = dlmerged
        iolarlite_path = dlmerged
        
        found = extract_event_info( event_info, iolcv_path, iolarlite_path, kpsreco_path )

        
        if found:
            break

        
    

