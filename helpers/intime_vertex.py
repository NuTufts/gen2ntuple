import os,sys

def get_nucandidate_intime_charge( nuvertex, wireplane_images_v, cosmictagged_pixels_v ):
    import ROOT
    from larflow import larflow

    masker = larflow.reco.ClusterImageMask()
    masker.storePixelValue( True )

    ntracks = nuvertex.track_hitcluster_v.size()
    for i in range(ntracks):
        trackcluster = nuvertex.track_hitcluster_v.at(i)
        masker.maskClusterAndStore( trackcluster, wireplane_images_v.as_vector(), 10.0, 2, False )
    print("[get_nucandidate_intime_charge] ntracks=",ntracks," npix-masked=",masker._npix)
    if masker._npix>0:
        for p in range(wireplane_images_v.as_vector().size()):
            print("  plane[",p,"]: ",masker.getPlaneMaskSum(p))
        
    # make an inverse mask using the cosmic pixels
    intime_charge = 0.0
    if masker._npix>0:
        print(cosmictagged_pixels_v.as_vector())
        masker.maskStoredImage( cosmictagged_pixels_v.as_vector(), 10.0, True )
        plane_charge = []
        print("[get_nucandidate_intime_charge] post-cosmictagger mask")
        for p in range( wireplane_images_v.as_vector().size() ):
            planesum = masker.getPlaneMaskSum(p)
            print("  plane[",p,"]: ",planesum)
            plane_charge.append( planesum )
        plane_charge.sort()
        intime_charge = plane_charge[-2]
    # take the median charge
    return intime_charge
    

