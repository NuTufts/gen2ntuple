#!/bin/bash

# Make Photon Selection Study Tree
# Usage: make_photon_selection_study_tree [OPTIONS]

# Photon Selection Variables for Study
# Creates possible variables for selecting photon events.

# Required Arguments:
#   --input-dlmerged FILE          DL Merged file holding wire plane larcv images
#   --input-reco FILE              LANTERN reco file holding neutrino candidates
#   --siren-model FILE             Siren Model weights (from make_siren_trace.py)
#   --output FILE                  Output ROOT file with matched data

# Optional Arguments:
#   --max-events N            Maximum number of events to process
#   --start-event N           Starting event number (default: 0)
#   --verbosity N             Verbosity larcv level 0-2 (default: 2. Most verbose=0)
#   --help                    Display this help message

# Examples:
#   make_photon_selection_study_tree --input-dlmerged dlmerged.root --input-reco output_lanternreco_kpsanafile.root --siren-model siren_extbnb.pt --output output_photonsel_variables.root


# lanternreco_mcc9_v40a_dl_run3b_NC_pi0_overlay_CV_000f0010-c197-47cf-bf02-7470b8e7cde5_kpsrecomanagerana.root
# lanternreco_mcc9_v40a_dl_run3b_NC_pi0_overlay_CV_000f0010-c197-47cf-bf02-7470b8e7cde5_larcv.root
# lanternreco_mcc9_v40a_dl_run3b_NC_pi0_overlay_CV_000f0010-c197-47cf-bf02-7470b8e7cde5_larlite.root
# larmatchme_easywave_mcc9_v40a_dl_run3b_NC_pi0_overlay_CV_000f0010-c197-47cf-bf02-7470b8e7cde5.root
# merged_dlana_000f0010-c197-47cf-bf02-7470b8e7cde5.root
# merged_dlreco_mcc9_v40a_dl_run3b_NC_pi0_overlay_CV_aa444faa-530a-4fd7-b43f-b501bc221880.root
# test.root

mergedfile="ncpi0/merged_dlana_000f0010-c197-47cf-bf02-7470b8e7cde5.root"
recobase="ncpi0/lanternreco_mcc9_v40a_dl_run3b_NC_pi0_overlay_CV_000f0010-c197-47cf-bf02-7470b8e7cde5_kpsrecomanagerana.root"
sirenweights="siren_model_corsika_lemon_snowflake_ckpt382k.pt"
outfile="test.root"

CMD="gdb --args ./build/installed/bin/make_photon_selection_study_tree --input-dlmerged ${mergedfile} --input-reco ${recobase} --siren-model ${sirenweights} --output ${outfile}"
CMD="${CMD} --verbosity 1"
#CMD="${CMD} --max-events 1"
echo ${CMD}
${CMD}

