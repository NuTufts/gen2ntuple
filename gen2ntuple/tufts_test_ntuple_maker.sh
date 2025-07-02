#!/bin/bash

dlmerged_input=/cluster/tufts/wongjiradlab/larbys/data/mcc9/mcc9_v28_wctagger_nueintrinsics/p00/out/v08_00_00_29/27939982_98/merged_dlreco_3989b1ff-c89a-41d2-8fe7-8b9bbd6af111.root
reco_input=/cluster/tufts/wongjiradlabnu/nutufts/data//v3dev_reco_retune/mcc9_v28_wctagger_nueintrinsics/larflowreco/ana/000/larflowreco_fileid0000_3989b1ff-c89a-41d2-8fe7-8b9bbd6af111_kpsrecomanagerana.root
larpid_model=/cluster/tufts/wongjiradlabnu/nutufts/larpid_weights/Scripted_LArPID_default_network_weights.pt
weight_file=/cluster/tufts/wongjiradlabnu/mrosen25/gen2ntuple/event_weighting/weights_forCV_v48_Sep24_intrinsic_nue_run1.pkl

  # # Process MC data with CNN
  # make_gen2_ntuples --mc \
  #   -f kpsreco_file1.root kpsreco_file2.root \
  #   -t merged_dlreco_list.txt \
  #   -m /path/to/model.pt \
  #   -w weights.pkl \
  #   -o output.root


make_gen2_ntuples --mc -f $reco_input -t $dlmerged_input -m $larpid_model -w $weight_file -o test.root -x photon_vertex_selection

