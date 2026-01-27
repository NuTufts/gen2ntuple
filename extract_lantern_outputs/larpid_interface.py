import os
import gc
from larcv import larcv
from larflow import larflow
import torch

from larpid_makeimage import makeImage

def run_larpid( nuvtx, iolcv, model ):
    """
    For a vertex, we run larpid on the track and shower prongs.
    We save not only the scores, but the crops used.
    """
    
    flowTriples = larflow.prep.FlowTriples()
    
    flowTriples.make_trackshower_images_from_sparse_uresnet_output( iolcv )
    mod_thrumu_v = flowTriples.make_thrumu_image_with_restored_ssnet_shower_pixels( iolcv, "ubspurn_plane", "thrumu" )    
    evtImage2D = iolcv.get_data("image2d", "wire")
    csmImage2D = iolcv.get_data("image2d", "thrumu")
    adc_v = evtImage2D.Image2DArray()
    thrumu_v = csmImage2D.Image2DArray()
    #cosmictagged_pixels_v = iolcv.get_data( larcv.kProductImage2D, "thrumu" )
    
    ntracks  = nuvtx.track_v.size()
    nshowers = nuvtx.shower_v.size()
    larpid_output = {}
    for iTrk in range(ntracks):
        trackCls = nuvtx.track_hitcluster_v.at(iTrk)
        nhits = trackCls.size()
        nTrajPoints = nuvtx.track_v[iTrk].NumberTrajectoryPoints()
        #trackLength = getDistance(vertex.track_v[iTrk].Vertex(),vertex.track_v[iTrk].End()) if (nTrajPoints > 1) else -9.
        trackLength = 1.0
        goodTrack = nTrajPoints > 1 and trackLength > 1e-6

        if goodTrack:
            cropPt = nuvtx.track_v[iTrk].End()
            print(" track loop[",iTrk,"] calling make_cropped_initial_sparse_prong_image_reco(...)",flush=True)
            prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,trackCls,cropPt,10.,512,512)
            with torch.no_grad():
                print("make prong image: ",prong_vv.size(),flush=True)
                #prongImage, prongImage_np = makeImage(prong_vv).to(args.device)
                prongImage_np = makeImage(prong_vv)
                print("run prongCNN on track image",flush=True)
                #prongCNN_out = model(prongImage)
                #trackClassified[iTrk] = 1
                #trackPID[iTrk] = getPID(prongCNN_out[0].argmax(1).item())
                #trackElScore[iTrk] = prongCNN_out[0][0][0].item()
                #trackPhScore[iTrk] = prongCNN_out[0][0][1].item()
                #trackMuScore[iTrk] = prongCNN_out[0][0][2].item()
                #trackPiScore[iTrk] = prongCNN_out[0][0][3].item()
                larpid_output[('track',iTrk)] = {'larpid_img':prongImage_np}
        else:
            #larpid_output[('track',iTrk)] = {}
            pass

    for iShw in range(nshowers):
        shower = nuvtx.shower_v.at(iShw)
        num_shower_hits = shower.size()
        cropPt = nuvtx.shower_trunk_v[iShw].Vertex()
        print(" shower loop[",iShw,"] calling make_cropped_initial_sparse_prong_image_reco(...)",flush=True)
        prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,mod_thrumu_v,shower,cropPt,10.,512,512)
        with torch.no_grad():
            print("   make prong image: ",prong_vv.size(),flush=True)
            #prongImage, prongImage_np = makeImage(prong_vv).to(args.device)
            prongImage_np = makeImage(prong_vv)
            print("   run prongCNN on shower image",flush=True)
            larpid_output[('shower',iShw)] = {'larpid_img':prongImage_np}

    # Clean up large intermediate objects
    del flowTriples
    del mod_thrumu_v
    del adc_v
    del thrumu_v
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    gc.collect()

    return larpid_output
