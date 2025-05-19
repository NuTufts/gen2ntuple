import os,sys

def highest_nu_keypoint_score(nuvetoed_v=None):
    """
    Select based on neutrino keypoints only. Keep the highest score.
    parameters
    ----------

    inputs
    ------
    nuvetoed_v [vector< larflow::reco::VertexNuCandidate >] Neutrino vertex candidates in the KPSRecoManagerTree TTree

    returns
    -------
    foundVertex [int] will be set to 1 if a vertex is found, 0 if not.
    vtxScore [float] set to -1. if no vertex is found, otherwise returns the maximum nu keypoint score. using 'netNuScore' value.
    vtxindex [int] index of max vertex in the nuvetoed_v container. set to -1, if no qualifying vertex found.
    """
    if nuvetoed_v is None:
        raise ValueError("No list of vertices passed to nu selection function 'highest_nu_keypoint_score'")
    
    foundVertex = 0
    vtxScore = -1.
    vtxindex = -1
    nvertices = nuvetoed_v.size()
    for ivtx in range(nvertices):
        vtx = kpst.nuvetoed_v.at(ivtx)
        if vtx.keypoint_type != 0:
            continue
        foundVertex = 1
        if vtx.netNuScore > vtxScore:
            vtxScore = vtx.netNuScore
            vtxindex = ivtx

    return foundVertex, vtxScore, vtxindex

def highest_kprank_with_visenergy(nuvetoed_v=None, nuselvar_v=None, min_num_showers=0, min_num_tracks=0):
    """
    Select based on keypoint type and visible energy.
    We favor neutrino keypoints (type==0).
    If multiple nu keypoints found, pick the one with the highest visible energy.
    If no nu keypoints found, pick the one with the highest visible energy.

    inputs
    ------
    nuvetoed_v [vector< larflow::reco::VertexNuCandidate >] Neutrino vertex candidates in the KPSRecoManagerTree TTree

    returns
    -------
    foundVertex [int] will be set to 1 if a vertex is found, 0 if not.
    vtxScore [float] set to -1. if no vertex is found, otherwise returns the maximum nu keypoint score. using 'netNuScore' value.
    vtxindex [int] index of max vertex in the nuvetoed_v container. set to -1, if no qualifying vertex found.
    """
    if nuvetoed_v is None:
        raise ValueError("No list of vertices passed to nu selection function 'highest_nu_keypoint_score'")

    # first test if there is a neutrino vertex
    foundNuVertex = 0
    nuvtxScore = -1.
    nuvtxindex = -1
    nvertices = nuvetoed_v.size()
    for ivtx in range(nvertices):
        vtx = nuvetoed_v.at(ivtx)
        selvar = nuselvar_v.at(ivtx)
        if vtx.keypoint_type != 0:
            continue
        if vtx.shower_v.size()<min_num_showers:
            continue
        if vtx.track_v.size()<min_num_tracks:
            continue
        foundNuVertex = 1
        # use nu score to rank
        #if vtx.netNuScore > vtxScore:
        #    nuvtxScore = vtx.netNuScore
        #    nuvtxindex = ivtx
        # use visible energy to rank
        if selvar.approx_vis_energy_MeV>nuvtxScore:
            nuvtxScore = selvar.approx_vis_energy_MeV
            nuvtxindex = ivtx
    if nuvtxindex>=0:
        # found neutrino vertex, return that.
        return foundNuVertex, nuvtxScore, nuvtxindex

    # if no neutrino vertex, find the vertex with the highest visible energy
    foundVertex = 0
    vtxScore = -1.
    vtxindex = -1
    for ivtx in range(nvertices):
        vtx = nuvetoed_v.at(ivtx)
        selvar = nuselvar_v.at(ivtx)        
        if vtx.keypoint_type == 0:
            # this time skip neutrino vertices
            continue
        if vtx.shower_v.size()<min_num_showers:
            continue
        if vtx.track_v.size()<min_num_tracks:
            continue        
        foundVertex = 1
        if selvar.approx_vis_energy_MeV > vtxScore:
            vtxScore = selvar.approx_vis_energy_MeV
            vtxindex = ivtx

    return foundVertex, vtxScore, vtxindex

def highest_intime_reco_frac( nuvetoed_v, nuselvar_v, prioritize_by_keypoint=True ):
    nvertices = nuvetoed_v.size()
    kp_intime_reco_frac = {0:10.0,
                           1:10.0,
                           2:10.0,
                           3:10.0}
    kp_index = {0:-1,
                1:-1,
                2:-1,
                3:-1}


    max_intime_reco_frac = 10.0
    max_kpindex = -1
    
    for ivtx in range(nvertices):
        nuvtx = nuvetoed_v.at(ivtx)
        if nuvtx.track_v.size()==0:
            continue
        if nuvtx.shower_v.size()==0:
            continue
        
        nusel = nuselvar_v.at(ivtx)
        kptype = int(nuvtx.keypoint_type)
        if kptype>3:
            continue
        frac_intime_reco = [ nusel.unreco_fraction_v[i] for i in range(3) ]
        frac_intime_reco.sort()
        if frac_intime_reco[1]<kp_intime_reco_frac[kptype]:
            kp_intime_reco_frac[kptype]=frac_intime_reco[1]
            kp_index[kptype] = ivtx
        if frac_intime_reco[1]<max_intime_reco_frac:
            max_kpindex = ivtx

    if prioritize_by_keypoint:
        # we go in kptype prioritization    
        for ikp in [0,3,1,2]:
            if kp_index[ikp]!=-1:
                max_idx = kp_index[ikp]
                kpscore = nuvetoed_v.at(max_idx).netScore
                print("select_nu_vertex.py: highest_intime_reco_frac : kptype=",ikp," idx=",max_idx," intime_reco_frac=",kp_intime_reco_frac[ikp])
                return 1, kpscore, kp_index[ikp]
            
    if max_kpindex>=0:
        kpscore = nuvetoed_v.at(max_kpindex).netScore
        return 1, kpscore, max_kpindex
    
    return 0, 0.0, -1

def select_nu_vertex(selector='highest_nu_keypoint_score', kwargs={}):
    if selector=='highest_nu_keypoint_score':
        return highest_nu_keypoint_score(**kwargs)
    elif selector=="highest_kprank_with_visenergy":
        return highest_kprank_with_visenergy(**kwargs)
    elif selector=="highest_intime_reco_frac":
        return highest_intime_reco_frac(**kwargs)
