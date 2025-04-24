import os,sys

from larlite import larlite
from larflow import larflow


input_larlite = sys.argv[1]

ioll = larlite.storage_manager( larlite.storage_manager.kREAD )
ioll.add_in_filename( input_larlite )
ioll.open()

outll = larlite.storage_manager( larlite.storage_manager.kWRITE )
outll.set_out_filename( "compressed_cosmic_components.root" )
outll.set_data_to_write( "track", "boundarycosmicreduced" )
outll.set_data_to_write( "track", "containedcosmicreduced" )
outll.set_data_to_write( "track", "cosmicprotonreduced" )
outll.open()

nentries = ioll.get_entries()

compressor  = larflow.recoutils.CompressRecoTrack()
max_saggita_cm = 1.0
max_step_size = 10.0


producers = ["boundarycosmicnoshift","containedcosmic","cosmicproton"]
for ientry in range(nentries):
    print("Process ENTRY[",ientry,"]")
    ioll.go_to(ientry)
    for producer in producers:
        ev_input_track = ioll.get_data( "track", producer )
        outname = producer+"reduced"
        if producer=="boundarycosmicnoshift":
            outname = "boundarycosmicreduced"
        ev_out_track = outll.get_data( "track", outname )
        print(outname)
        ninput_tracks = 0
        try:
            ninput_tracks = ev_input_track.size()
        except:
            ninput_tracks = 0
        for itrack in range( ninput_tracks ):
            track = ev_input_track.at( itrack )
            reduced_track = compressor.compress( track, max_saggita_cm, max_step_size )
            norig = track.NumberTrajectoryPoints()
            nreduced = reduced_track.NumberTrajectoryPoints()
            print(f"  {outname}[{itrack}]: {norig} --> {nreduced}")
            ev_out_track.push_back( reduced_track )
    outll.set_id( ioll.run_id(), ioll.subrun_id(), ioll.event_id() )
    outll.next_event()
    
ioll.close()
outll.close()

