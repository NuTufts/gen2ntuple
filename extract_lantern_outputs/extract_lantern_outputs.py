import os,sys
import ROOT as rt

from filelists import get_filelist

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

    # Get the nu reco object.
    rfile_reco = rt.TFile( kpsreco_path, 'read' )
    #rfile_reco.Get()

    recotree = rfile_reco.Get("KPSRecoManagerTree")
    nentries = recotree.GetEntries()
    found  = False
    for ientry in range(nentries):
        recotree.GetEntry(ientry)
        if run!=recotree.run or subrun!=recotree.subrun or event!=recotree.event:
            continue
        nuvertices = recotree.nuvetoed_v.size()
        if nuvertices<=vtxIdx:
            raise ValueError(f"Number of vertices in this event ({nvertices}) is less than target vertex index ({vtxIdx})")
        nuvtx = recotree.nuvetoed_v.at(vtxIdx)
        found = True
        print(f"Found vertex index: {vtxIdx}")
        
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

        
    

