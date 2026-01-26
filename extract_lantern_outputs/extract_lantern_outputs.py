import os,sys
import numpy as np
import h5py
import ROOT as rt
from larcv import larcv
from larflow import larflow

# larpid code
# ./../../photon_analysis/prongCNN/models/

from filelists import get_filelist,NTUPLE_FILES
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

def get_ntuple_rsemap( dataset_name ):
    ntuple_path = NTUPLE_FILES[dataset_name]
    print("Making RSE Map for: ",ntuple_path)
    rfile = rt.TFile( ntuple_path, 'read' )
    ntuple = rfile.Get("EventTree")
    ntuple.SetBranchStatus("*",0)
    ntuple.SetBranchStatus("run",1)
    ntuple.SetBranchStatus("subrun",1)
    ntuple.SetBranchStatus("event",1)
    rsemap = {}
    runmax = -1
    runmin = 10000000
    for ientry in range(ntuple.GetEntries()):
        ntuple.GetEntry(ientry)
        if ntuple.run>runmax:
            runmax = ntuple.run
        if ntuple.run<runmin:
            runmin = ntuple.run
        rsemap[ (int(ntuple.run),int(ntuple.subrun),int(ntuple.event)) ]  = ientry
    rfile.Close()
    print(f" RSE Map made. Run Min={runmin} Max={runmax}")
    return rsemap

def extract_spacepoints( nuvtx ):
    ntracks  = nuvtx.track_v.size()
    nshowers = nuvtx.shower_v.size()
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

    if len(prong_spacepoints) > 0:
        spacepoints = np.concatenate( prong_spacepoints, axis=0 )
    else:
        spacepoints = np.zeros( (0,5) )
    return spacepoints, ntracks, nshowers


def write_event_to_hdf5( h5file, event_info, spacepoints, larpid_output,
                         trackinfo, showerinfo, ntracks, nshowers ):
    """
    Write event data to an HDF5 file.

    Structure:
        /event_<idx>/
            event_info/  (attributes: run, subrun, event, fileid, vertexid)
            spacepoints  (N, 5) array: [x, y, z, prong_id, prong_type]
            ntracks      scalar
            nshowers     scalar
            track_images/
                track_0  (6, 512, 512) array
                track_1  ...
            shower_images/
                shower_0 (6, 512, 512) array
                shower_1 ...
            track_info/
                track_0  (5,) array
                track_1  ...
            shower_info/
                shower_0 (5,) array
                shower_1 ...
    """
    # Create event group
    run = event_info['run']
    subrun = event_info['subrun']
    event = event_info['event']
    fileid = event_info['fileid']
    vtxIdx = event_info['vtxIdx']

    event_group_name = f"event_{fileid}_{run}_{subrun}_{event}_{vtxIdx}"
    event_grp = h5file.create_group(event_group_name)

    # Store event indexing info as attributes
    event_grp.attrs['run'] = run
    event_grp.attrs['subrun'] = subrun
    event_grp.attrs['event'] = event
    event_grp.attrs['fileid'] = fileid
    event_grp.attrs['vertexid'] = vtxIdx
    event_grp.attrs['ntracks'] = ntracks
    event_grp.attrs['nshowers'] = nshowers

    # Store spacepoints
    event_grp.create_dataset('spacepoints', data=spacepoints, compression='gzip')

    # Create groups for track and shower images
    track_grp = event_grp.create_group('track_images')
    shower_grp = event_grp.create_group('shower_images')

    # Create groups for track and shower images
    trackinfo_grp = event_grp.create_group('track_info')
    showerinfo_grp = event_grp.create_group('shower_info')
    
    # Store images by prong type
    for k, vdict in larpid_output.items():
        prongtype, prongindex = k
        if 'larpid_img' in vdict:
            prongimages = vdict['larpid_img']
            if prongtype == 'track':
                track_grp.create_dataset(f'track_{prongindex}',
                                         data=prongimages,
                                         compression='gzip')
                trackinfo_grp.create_dataset(f'track_{prongindex}',
                                             data=trackinfo[trackinfo[:,1]==prongindex,:],
                                             compression='gzip')
            elif prongtype == 'shower':
                shower_grp.create_dataset(f'shower_{prongindex}',
                                          data=prongimages,
                                          compression='gzip')
                showerinfo_grp.create_dataset(f'shower_{prongindex}',
                                              data=showerinfo[showerinfo[:,1]==prongindex,:],
                                              compression='gzip')
    

    return event_group_name


def extract_event_info( event_info, iolcv_path, iolarlite_path, kpsreco_path,
                        ntuple_path, ntuple_rsemap,
                        h5file=None ):
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

    # Load LArCV data file
    iolcv = larcv.IOManager(larcv.IOManager.kREAD,'larcv',larcv.IOManager.kTickBackward)
    iolcv.add_in_file( iolcv_path )
    iolcv.reverse_all_products()
    iolcv.initialize()

    # Load ntuple
    ntuple_rfile = rt.TFile( ntuple_path, 'read' )
    ntuple_tree = ntuple_rfile.Get("EventTree")

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
        prong_spacepoints, ntracks, nshowers = extract_spacepoints( nuvtx )

        # Get Prong Data from the LANTERN ntuple
        rse = (run,subrun,event)
        ntuple_entry = ntuple_rsemap[rse]
        ntuple_bytes = ntuple_tree.GetEntry(ntuple_entry)
        if ntuple_bytes==0:
            raise ValueError(f"Could not find ntuple entry for RSE={rse}")
        # We make arrays for the tracks and showers
        # we store flags for each particle:
        #  [0]: 0 (track) or 1 (shower)
        #  [1]: index in track or shower array        
        #  [2]: 0 (primary) or 1 (secondary)
        #  [3]: pdgcode from lantern
        trackinfo_array  = np.zeros( (ntuple_tree.nTracks,4), dtype=np.int64 )
        showerinfo_array = np.zeros( (ntuple_tree.nShowers,4), dtype=np.int64 )        
        for i in range(ntuple_tree.nTracks):
            trackinfo_array[i,0] = 0
            trackinfo_array[i,1] = i
            trackinfo_array[i,2] = ntuple_tree.trackIsSecondary[i]
            trackinfo_array[i,3] = ntuple_tree.trackPID[i]
        for i in range(ntuple_tree.nShowers):
            showerinfo_array[i,0] = 1
            showerinfo_array[i,1] = i
            showerinfo_array[i,2] = ntuple_tree.showerIsSecondary[i]
            showerinfo_array[i,3] = ntuple_tree.showerPID[i]

        # Write to HDF5 if file handle provided
        if h5file is not None:
            event_group_name = write_event_to_hdf5(
                h5file, event_info, prong_spacepoints, larpid_output,
                trackinfo_array, showerinfo_array, ntracks, nshowers
            )
            print(f"Wrote event data to HDF5 group: {event_group_name}")
        
        
        if True:
            break
        
    rfile_reco.Close()
    return found

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Extract LArTPC event data to HDF5')
    parser.add_argument('-s', '--selection-file', type=str,
                        default="mcc9_v28_wctagger_run3_bnb1e19_eventselection_20260126_175956.root",
                        help='ROOT file with selected event info')
    parser.add_argument('-l', '--filelist', type=str,
                        default="run3_opendata_1e19",
                        help='Name of filelist to use')
    parser.add_argument('-o', '--output', type=str,
                        default="lantern_output.h5",
                        help='Output HDF5 filename')
    parser.add_argument('--max-events', type=int, default=-1,
                        help='Maximum number of events to process (-1 for all)')
    args = parser.parse_args()

    event_info_list = get_selected_event_info_fromrootfile( args.selection_file )

    fileid_dict = get_filelist( args.filelist )

    ntuple_path = NTUPLE_FILES[args.filelist]
    ntuple_rsemap = get_ntuple_rsemap( args.filelist )
    print("RSEMap Entries: ",len(ntuple_rsemap))

    # Create HDF5 output file
    with h5py.File(args.output, 'w') as h5file:
        # Store metadata
        h5file.attrs['selection_file'] = args.selection_file
        h5file.attrs['filelist'] = args.filelist

        nevents_processed = 0
        for event_info in event_info_list:
            fileid = event_info['fileid']
            fileid_files = fileid_dict[fileid]
            if "reco" not in fileid_files:
                continue

            dlmerged = fileid_files['dlmerged']
            kpsreco_path = fileid_files['reco']
            iolcv_path = dlmerged
            iolarlite_path = dlmerged

            found = extract_event_info( event_info, iolcv_path, iolarlite_path,
                                        kpsreco_path, ntuple_path, ntuple_rsemap, h5file=h5file )

            if found:
                nevents_processed += 1
                if args.max_events > 0 and nevents_processed >= args.max_events:
                    print(f"Reached max events limit: {args.max_events}")
                    break

        h5file.attrs['nevents'] = nevents_processed
        print(f"Processed {nevents_processed} events. Output written to: {args.output}")
