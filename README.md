# gen2ntuple
scripts for making and analyzing DL Gen2 ntuples

## software setup

### dependencies

* ROOT is needed for all scripts
* pytorch is needed for the LArPID network
* The ubdl software environment must be setup to run the ntuple maker (`make_dlgen2_flat_ntuples.py`). See the quick start section of the readme in the [ubdl repository](https://github.com/LArbys/ubdl) for instructions

ROOT, pytorch, and all of the dependencies needed to compile the ubdl repository are available in a singularity container Taritree has prepared, which you can find on the Tufts cluster at `/cluster/tufts/wongjiradlabnu/larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif`. To load the container at Tufts, do:

```
module load singularity/3.5.3
singularity shell --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab /cluster/tufts/wongjiradlabnu/larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
```

Then enter `bash` when prompted to select a bash shell for the container. 

### first time setup

Before running the ntuple maker, you need to compile a C++ script in the `helpers` directory (this creates a shared C++ library used for the fiducial volume calculation). To do this, simply run:

```
cd helpers/
source compile_wirecell_fiducial_volume.sh
```

This only needs to be done once after first downloading the repository.

## Running the ntuple maker

The `make_dlgen2_flat_ntuples.py` script produces ntuple files from input larflowreco files and merged_dlreco files. It allows or requires the following options:

* __-f__ / __--files__: The input larflowreco files that will be processed (__required__). You can either provide a text file where each line in the file is the full path to a larflowreco file, the path to a single larflowreco file, or a space separated list of larflowreco files.
* __-t__ / __--truth__: The merged_dlreco files matching the input larflowreco files (__required__). This must be a text file where each line is the full path to a merged_dlreco file. This text file can include a full list of all merged_dlreco files in some sample, it does not need to be restricted to the merged_dlreco files that match the larflowreco files supplied with --files. The script will go through the input merged_dlreco files and find the ones with the same hash identifier as the input larflowreco files.
* __-w__ / __--weightfile__: The path to the file (a pickled python dictionary) containing the cross section weights (__required for MC inputs__). Check the `event_weighting` directory to see if there is a weights file there that matches your sample. If not, there are instructions in that directory on how to create one. This should not be supplied for data input but must be supplied for MC input.
* __-m__ / __--model_path__: The path to the pytorch checkpoint file containing the weights to run the LArPID network. The path to the most recent weights file on the Tufts cluster is `/cluster/tufts/wongjiradlabnu/mrosen25/prongCNN/models/checkpoints/run3bOverlays_quadTask_plAll_2inChan_5ClassHard_minHit10_b64_oneCycleLR_v2me05_noPCTrainOrVal/ResNet34_recoProng_5class_epoch20.pt`
* __-d__ / __--device__: The device to run the LArPID network on ("cpu" or "gpu"). Default: "cpu"
* __-mc__ / __--isMC__: A flag that __must be supplied when running over MC input__
* __-ana__ / __--dlana_input__: A flag that __must be supplied when running over dlana input__ (truth files that begin with "merged_dlana" instead of "merged_dlreco"
* __-o__ / __--outfile__: The name of the output ntuple file (default: "dlgen2_flat_ntuple.root")
* __-nkp__ / __--noKeypoints__: A flag, when supplied the keypoint variables will not be saved in the output ntuple. __Required when running over larflowreco files with version < v2_me_06__. (Older larflowreco files don't have the keypoint info.)
* __--ignoreWeights__: A flag, if supplied don't look up xsec weights for MC events. Set weights to 1 and process all events. (Default: attempt to lookup xsec weights and exit with error if not found)
* __--skipNoWeightEvts__: A flag, if supplied: print a warning, skip the event, and continue processing when we can't find an xsec weight for an MC event. (Default: attempt to lookup xsec weights and exit with error if not found)
* __--multiGPU__: A flag, supply to run LArPID network on multiple gpus. (This really isn't necessary, the network runs pretty quickly on just a cpu.)

There are two scripts included in this repository to provide an example of how to submit ntuple maker jobs on the Tufts cluster using slurm:
* tufts_submit_ntuple_job_example.sh
* tufts_run_ntuple_maker.sh

We produce ntuple files as part of the standard processing, but if you want to add a variable or two and reprocess the ntuples from existing larflowreco files, then you can take a look at these scripts for an example on how that can be done. If you'd like to use them, then __first copy__ both scripts and then edit the following lines:
* __tufts_submit_ntuple_job_example.sh__
    * line 18: supply the full file path to a text file containing all the larflowreco files you'd like to process. Each line should be the full path to one larflowreco file
    * line 19: supply the full file path to a text file containing all the merged_dlreco or merged_dlana files that match the larflowreco files. Each line should be the full path to one merged_dlreco or merged_dlana file. The file on a particular line of this list need not match the larflowreco file on the same line of the larflowreco list. The scripts select larflowreco files to process, then searches through the supplied merged_dlreco/dlana list to find the matching files.
    * line 13: change this NFILES variable to the number of files you want to process per job. I usually do 100
    * lines 3-11: edit the slurm sbatch options as necessary. You'll want to at least change the `--array` option, which dictates how many jobs are submitted. This option works with the NFILES variable on line 13 to determine which larflowreco files from the line 18 list are processed in which job. For example, using `--array=0-3,6` and `NFILES=100` will submit five jobs (with arrayIDs 0, 1, 2, 3, and 6), each of which will process 100 files given by those specified on lines arrayID\*NFILES + 1 to arrayID\*NFILES + NFILES.
    * line 17: provide the full file path to the `tufts_run_ntuple_maker.sh` script
    * line 20: the filename (not path!) of the weight file in the event_weighting directory that will get passed to `make_dlgen2_flat_ntuples.py` with the `-w` option (see above).
    * line 21: supply an output tag that will be included in output file names
    * line 26: edit, as necessary, the end of this line where you see `-mc`. In this example, `-mc` (and any other options you might add to the end of the source command after it) is passed to the `make_dlgen2_flat_ntuples.py` script as a flag. So, for example, if you were running over data, you would remove `-mc`. If you were running over data with merged_dlana input and you didn't want to save keypoint information, you would change `-mc` to `-ana -nkp`. Note that `tufts_run_ntuple_maker.sh` is currently only configured to accept up to two flags supplied in this way. You'll have to edit lines 59-68 of `tufts_run_ntuple_maker.sh` to get it to accept more.
* __tufts_run_ntuple_maker.sh__
    * line 5: supply the path to your copy of this gen2ntuple repository
    * line 22: supply the path to your compiled and installed copy of the ubdl repository
    * line 41: supply the directory you would like output files copied to
    * line 42: supply the directory you would like log files copied to

To submit the jobs after the above editing, just do:
```
sbatch modified_tufts_submit_ntuple_job_example.sh
```

All of the files processed in a single job will, at the end, get merged together before being copied over to the specified output directory. So each job will produce one output file named make_dlgen2_flat_ntuples_OUTTAG_output_ARRAYID.root and one log file (in addition to the usual slurm .out file) named make_dlgen2_flat_ntuples_OUTTAG_ARRAYID.log, where OUTTAG is the output tag (specified on line 21 of `tufts_submit_ntuple_job_example.sh`) and ARRAYID is the slurm array ID (all specified on line 8 of `tufts_submit_ntuple_job_example.sh` will be run) of the job. You can use the root `hadd` command to merge the output root files

## using ntuple files

See the following section for full documentation of all ntuple variables.

The ntuples are flat ROOT files containing information about reconstructed neutrino interactions and (for Monte Carlo samples) truth information about the simulated neutrino interactions. If you are new to ROOT, this [ROOT Guide For Beginners](https://root.cern.ch/root/htmldoc/guides/primer/ROOTPrimerLetter.pdf) provides a nice, detailed introduction.

For Monte Carlo samples, the ntuples also contain information about the number of protons on target (POT) that would produce the full sample contained in the ntuple.

I have provided an example script, `example_ntuple_analysis_script.py`, that you can take a look at to see how to use the POT data and how to pull some basic information out of the ntuple files using ROOT in python (PyROOT). This script takes as input an ntuple file made from an MC bnb nu overlay sample, applies an inclusive CCnumu selection, makes histograms of true and reconstructed neutrino and primary muon energies, scales the output histograms to the number of events expected for 6.67e+20 POT, and writes these histograms to an output root file. You can supply the flags `--fullyContained` to filter out events that are not fully contained inside the fiducial volume and `--noCosmicCuts` to skip the cosmic rejection cuts in the CCnumu selection.

To run:
```
python example_ntuple_analysis_script.py -i input_ntuple_file.root -o output_histogram_file.root
```

You can then view the histograms in `output_histogram_file.root` in e.g. a ROOT TBrowser:
```
root -l output_histogram_file.root
root [0] new TBrowser
```

## ntuple variable documentation

Please update this readme if you add new variables!

### Acronyms / terms

* __MC__: Monte Carlo simulation
* __CC__: charged current
* __NC__: neutral current
* __QE__: quasi-elastic
* __Res__: resonance production
* __DIS__: deep inelastic scattering
* __Coh__: coherent 
* __MEC__: meson exchange current
* __POT__: protons on target
* __PDG code__: particle data group code, an integer unique to every type of particle
* __SCE__: space charge effect
* __True__: refers to simulated quantities from MC
* __Reco / reconstructed__: refers to quantities calculated from detector measurements (or simulated detector measurements in the case of MC events)
* __Detsim tracked__: simulated particle tracked by the detector simulation
* __trackID__: a unique integer associated with a simulated particle. Used to identify a particular simulated particle in an event
* __Prong__: a reconstructed track or shower
* __Truth-matched particle__: the simulated particle that deposits the most energy in a reconstructed prong
* __Prong purity__: the fraction of a reco prong produced by it’s truth-matched particle. (Calculated by summing pixel charges for pixels associated with the prong’s reco 3D hits.)
* __Prong completeness__: the fraction of a reco prong’s truth-matched particle that is included in the reco prong. (Calculated by summing pixel charges for pixels associated with the prong’s reco 3D hits.)
* __LArPID__: a convolutional neural network used for particle identification, particle production process classification, and purity and completeness predictions for reconstructed prongs
* __Prong charge__: the sum of all pixel values for pixels associated with prong’s reco 3D hits
* __Hit__: a reconstructed 3D point created from measured/simulated wire plane signals
* __Fiducial Volume__: the region of the detector where we select vertices / events for analyses
* __Wire Cell Fiducial Volume__: the fiducial volume defined in the wire-cell analysis (> 3cm from SCE corrected detector edge)

### Event filters:

* An entry will appear in the flat ntuple files for every input data event
* For MC simulation, events that meet the following criteria are excluded and will not appear in these file:
    * SCE corrected neutrino interaction vertex is outside the wire cell fiducial volume (criteria removed 04/2025)
    * The event’s cross-section weight is unknown or infinite (should never happen once weights files are updated)

### potTree variables:

The potTree (present only in ntuples produced from MC samples) provide information about the number of protons on target (POT) that would produce the full sample contained in the ntuple. There is one entry for each merged_dlreco file used to create the ntuple file. To get the full POT count, add up all the entries in the tree.

Variables:

* __totPOT__: total POT for the file represented by this entry
* __totGoodPOT__: total good POT for the file represented by this entry

I am not entirely sure what the distinction is between the "total POT" and "total good POT" reported in merged_dlreco files, but the two are always the same in every sample I've looked at. I would recomment using `totGoodPOT` and ignoring the `totPOT` entries.

### EventTree file/event identification variables:

* __fileid__: File ID number for larflowreco file that this event was from. Can be used to locate the larflowreco and merged_dlreco files for a particular event
* __run__: run number for event
* __subrun__: subrun number for event
* __event__: event number for event

### EventTree simulation variables for MC events:

These variables are not included for data.

* __xsecWeight__: cross section weight for MC events, use to weight events for POT scaling.
* __trueNuE__: true energy in GeV of simulated neutrino.
* __trueNuPDG__: PDG code of simulated neutrino.
* __trueNuCCNC__: integer indicating if simulated neutrino interaction is CC or NC. 0 for CC, 1 for NC.
* __trueNuMode__: integer indicating the simulated neutrino interaction mode: 0 for QE, 1 for Res, 2 for DIS, 3 for Coh, 10 for MEC.
* __trueNuIntrxnType__: integer indicating the simulated neutrino interaction type. A list of categories corresponding to values can be found [here](https://github.com/NuTufts/larlite/blob/7aadee7f8a0c72429774578539842e8d4936d457/larlite/Base/MCConstants.h). 
* __trueVtx{X,Y,Z}__: SCE corrected X,Y,Z coordinates (in cm) of simulated neutrino interaction vertex.
* __trueVtxIsFiducial__: 1 if SCE corrected simulated neutrino interaction vertex is inside Wire-Cell fiducial volume, 0 otherwise.
* __trueLepE__: true energy in GeV of primary lepton.
* __trueLepPDG__: PDG code of simulated primary lepton.
* __nTruePrimParts__: number of true primary particles (number of stable final state particles from the neutrino interaction that are passed to the detector simulation)
* __truePrimPartPDG__: PDG code of true primary particle. Array of length nTruePrimParts
* __truePrimPart{X,Y,Z}__: SCE corrected X,Y,Z coordinates (in cm) of true primary particle start position. Array of length nTruePrimParts
* __truePrimPart{Px,Py,Pz,E}__: initial 4-momentum of true primary particle in GeV/c for 3-momentum coordinates and GeV for total energy. Array of length nTruePrimParts
* __truePrimPartContained__: integer indicating if true primary particle was contained inside the detector. 1 if SCE corrected true end position is inside wire cell fiducial volume, 0 otherwise
* __nTrueSimParts__: number of simulated detsim-tracked particles (particles tracked by the detector simulation)
* __trueSimPartPDG__: PDG code of detsim-tracked particle. Array of length nTrueSimParts
* __trueSimPartTID__: trackID of detsim-tracked particle. Array of length nTrueSimParts
* __trueSimPartMID__: trackID of detsim-tracked particle’s “mother” (ancestor particle that e.g. decayed, producing this particle). Array of length nTrueSimParts
* __trueSimPartProcess__: Integer indicating the process by which a detsim-tracked particle was created. 0 for primary particles from neutrino interaction, 1 for particles produced in a decay, 2 for all other processes
* __trueSimPart{X,Y,Z}__: SCE corrected X,Y,Z coordinates (in cm) of true detsim-tracked particle start position. Array of length nTrueSimParts
* __trueSimPartEDep{X,Y,Z}__: SCE corrected X,Y,Z coordinates (in cm) of point at which detsim-tracked particle begins depositing energy. Point of first photo-conversion for photons, identical to trueSimPart{X,Y,Z} for all other particles. Array of length nTrueSimParts
* __trueSimPart{Px,Py,Pz,E}__: initial 4-momentum of true primary particle in MeV/c for 3-momentum coordinates and MeV for total energy. Array of length nTrueSimParts
* __trueSimPartContained__: integer indicating if detsim-tracked particle was contained inside the detector. 1 if SCE corrected true end position is inside wire cell fiducial volume, 0 otherwise

### EventTree reconstruction variables:

#### Vertex/event variables

* __foundVertex__: 1 if a neutrino interaction vertex was reconstructed, 0 otherwise. (Technical details: a neutrino interaction vertex was considered to be reconstructed if there was at least one LArMatch keypoint with a neutrino keypoint type. If there are multiple neutrino keypoints, the one with the highest keypoint score is selected as the reconstructed neutrino vertex.)
* __vtx{X,Y,Z}__: X,Y,Z coordinates (in cm) of reco neutrino vertex. Defaults to -999 if foundVertex = 0
* __vtxIsFiducial__: 1 if the reco neutrino vertex is inside the wire cell fiducial volume, 0 otherwise. Defaults to -1 if foundVertex = 0
* __vtxContainment__: 2 if all hits in all prongs associated with the reco neutrino vertex are inside the wire cell fiducial volume. 1 if at least one hit in a neutrino-vertex-associated prong is outside the wire cell fiducial volume, but the reco neutrino vertex is inside. 0 if the reco neutrino vertex is outside the wire cell fiducial volume. Defaults to -1 if foundVertex = 0
* __vtxScore__: keypoint score of reco neutrino vertex. Defaults to -1 if foundVertex = 0
* __vtxFracHitsOnCosmic__: fraction of all hits associated with the reco neutrino vertex that match cosmic-tagged pixels
* __vtxDistToTrue__: distance (in cm) between reco and SCE-corrected true neutrino interaction vertex. Not present for data
* __recoNuE__: reconstructed neutrino energy in MeV. Calculated by summing the reconstructed energies of all tracks and showers attached to the reco neutrino vertex
* __eventPCAxis{0,1,2}__: 3-vectors indicating the direction of the first, second, and third principal component axes in a PCA analysis using the coordinates of all hits in all prongs associated with the reconstructed neutrino vertex
* __eventPCEigenVals__: array of three numbers: the eigenvalues assocated with the PCA analysis described above (e.g. eventPCEigenVals[1] gives the eigenvalue associated with eventPCAxis1)
* __eventPCAxis0TSlope__: diagnostic variable verifying that the first PCA axis (eventPCAxis0) is pointing in direction of increasing hit times (1 if yes, -1 if not)
* __eventPCProjMaxGap__: array of five numbers indicating the largest gap (in cm) of hit projections onto the first PCA axis (eventPCAxis0). All hits in all prongs associated with the reco nuetrino vertex are used for the 0th element (eventPCProjMaxGap[0]). In the 1st, 2nd, 3rd, and 4th elements, the 10%, 20%, 30%, and 40%, respectively, of hits furthest from the vertex are excluded. All elements default to -9 if foundVertex = 0
* __eventPCProjMaxDist__: array of five numbers indicating the largest distance between hit projections onto the first PCA axis (eventPCAxis0) without a gap of distance greater than X, where X = 2cm, 4cm, 6cm, 8cm, or 10cm for the 0th (eventPCProjMaxDist[0]), 1st, 2nd, 3rd, and 4th elements, respectively. All hits in all prongs associated with the reco neutrino vertex are used. All elements default to -9 if foundVertex = 0
* __nKeypoints__: number of LArMatch keypoint clusters for this event
* __kpClusterType__: keypoint cluster type (0: kNuVertex, 1: kTrackStart, 2: kTrackEnd, 3: kShowerStart, 4: kShowerMichel, 5: kShowerDelta, 6: kVertexActivity, 7: kStopMuVtx, 8: kNumKeyPoints). Array of length nKeypoints
* __kpFilterType__: 0 if neutrino activity keypoint cluster, 1 if cosmic keypoint cluster. Array of length nKeypoints
* __kpMaxScore__: Maximum score of any keypoint in cluster. Array of length nKeypoints
* __kpMaxPos{X,Y,Z}__: X,Y,Z coordinates (in cm) of the keypoint in cluster with the highest score. Array of length nKeypoints

#### Track variables

* __nTracks__: number of reco tracks attached to reco neutrino vertex. Defaults to 0 if foundVertex = 0
* __trackIsSecondary__: 1 if track was attached to reco neutrino vertex as a secondary particle (a particle produced downstream from the interaction vertex as a result of e.g. the decay of a primary particle produced in the neutrino interaction itself). Array of length nTracks
* __trackNHits__: number of reco 3D hits in track. Array of length nTracks
* __trackHitFrac__: trackNHits divided by total number of hits associated with reco neutrino vertex. Array of length nTracks
* __trackCharge__: Sum of pixel values for all above-noise-threshold pixels associated with track hits (noise threshold: 10). Array of length nTracks
* __trackChargeFrac__: trackCharge divided by sum of all above-noise-threshold pixels associated with reco neutrino vertex. Array of length nTracks
* __trackCosTheta__: cosine of angle between track start direction and beam (z-axis). Track start direction calculated as unit vector pointing from track start position to point on track 5cm from start position (or end of track if length < 5cm). Defaults to -9 if track has fewer than two trajectory points or a length of < 1e-6 cm. Array of length nTracks
* __trackCosThetaY__: cosine of angle between track start direction and gravity (y-axis). Track start direction calculated as unit vector pointing from track start position to point on track 5cm from start position (or end of track if length < 5cm). Defaults to -9 if track has fewer than two trajectory points or a length of < 1e-6 cm. Array of length nTracks
* __trackDistToVtx__: distance (in cm) between track start point and reco neutrino vertex. Defaults to -9 if track has fewer than two trajectory points or a length of < 1e-6 cm. Array of length nTracks
* __trackStartPos{X,Y,Z}__: X,Y,Z coordinates (in cm) of track start position. Defaults to -9 if track has fewer than two trajectory points or a length of < 1e-6 cm. Array of length nTracks
* __trackStartDir{X,Y,Z}__: Coordinates of unit vector pointing in initial track direction, calculated as the vector pointing from the beginning of the track to the first trajectory point at least 5cm from track start (or end of track if track length is less than 5cm). Defaults to -9 if track has fewer than two trajectory points or a length of < 1e-6 cm. Array of length nTracks
* __trackEndPos{X,Y,Z}__: X,Y,Z coordinates (in cm) of track end position. Defaults to -9 if track has fewer than two trajectory points or a length of < 1e-6 cm. Array of length nTracks
* __trackClassified__: 1 if track was classified, 0 otherwise. Track is “classified” if it was processed by LArPID, which will happen if the track has at least 2 trajectory points, a length of >1e-6 cm, and at least 10 above-noise-threshold pixels in all three wire planes (noise threshold: 10). Array of length nTracks
* __trackPID__: Predicted PDG code of track, from LArPID. Defaults to 0 if trackClassified = 0. Array of length nTracks
* __trackElScore__: Track’s electron score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackPhScore__: Track’s photon score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackMuScore__: Track’s muon score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackPiScore__: Track’s pion score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackPrScore__: Track’s proton score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackComp__: Track’s predicted completeness from LArPID. Defaults to -1 if trackClassified = 0. Array of length nTracks
* __trackPurity__: Track’s predicted purity from LArPID. Defaults to -1 if trackClassified = 0. Array of length nTracks
* __trackProcess__: Track's predicted particle production process from LArPID (0: primary particle from neutrino interaction, 1: secondary particle with a neutral parent, 2: secondary particle with a charged parent). Defaults to -1 if trackClassified = 0. Array of length nTracks
* __trackPrimaryScore__: Track's primary particle score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackFromNeutralScore__: Track's secondary particle with neutral parent score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackFromChargedScore__: Track's secondary particle with charged parent score from LArPID. Defaults to -99 if trackClassified = 0. Array of length nTracks
* __trackRecoE__: reconstructed energy of track in MeV, calculated from range assuming it is a muon, proton, or pion: whichever has the highest LArPID score. Defaults to -1 if trackClassified = 0. Array of length nTracks
* __trackTruePID__: The PDG code of the track’s truth-matched particle. Defaults to 0 if no simulated particle deposited energy in track's pixels. Not filled for data events. Array of length nTracks
* __trackTrueTID__: The trackID of the track’s truth-matched particle. Defaults to -1 if no simulated particle deposited energy in track's pixels. Not filled for data events. Array of length nTracks
* __trackTrueE__: The true energy (in MeV) of the track’s truth-matched particle. Defaults to -1 if no simulated particle deposited energy in track's pixels. Not filled for data events. Array of length nTracks
* __trackTruePurity__: The true purity of the reco track, calculated from its truth-matched particle. Defaults to 0 if no simulated particle deposited energy in track's pixels. Not filled for data events. Array of length nTracks
* __trackTrueComp__: The true completeness of the reco track, calculated from its truth-matched particle. Defaults to 0 if no simulated particle deposited energy in track's pixels. Not filled for data events. Array of length nTracks
* __trackTrueElPurity__: The fraction of the charge in the reco track produced by simulated electrons. Not filled for data events. Array of length nTracks
* __trackTruePhPurity__: The fraction of the charge in the reco track produced by simulated photons. Not filled for data events. Array of length nTracks
* __trackTrueMuPurity__: The fraction of the charge in the reco track produced by simulated muons. Not filled for data events. Array of length nTracks
* __trackTruePiPurity__: The fraction of the charge in the reco track produced by simulated pions. Not filled for data events. Array of length nTracks
* __trackTruePrPurity__: The fraction of the charge in the reco track produced by simulated protons. Not filled for data events. Array of length nTracks

#### Shower variables

* __nShowers__: number of reco showers attached to reco neutrino vertex. Defaults to 0 if foundVertex = 0
* __showerIsSecondary__: 1 if shower was attached to reco neutrino vertex as a secondary particle (a particle produced downstream from the interaction vertex as a result of e.g. the decay of a primary particle produced in the neutrino interaction itself). Array of length nShowers
* __showerNHits__: number of reco 3D hits in shower. Array of length nShowers
* __showerHitFrac__: showerNHits divided by total number of hits associated with reco neutrino vertex. Array of length nShowers
* __showerCharge__: Sum of pixel values for all above-noise-threshold pixels associated with shower hits (noise threshold: 10). Array of length nShowers
* __showerChargeFrac__: showerCharge divided by sum of all above-noise-threshold pixels associated with reco neutrino vertex. Array of length nShowers
* __showerCosTheta__: cosine of angle between shower start direction and beam (z-axis). Shower start direction calculated as unit vector pointing from shower trunk start position to shower trunk end position. Array of length nShowers
* __showerCosThetaY__: cosine of angle between shower start direction and gravity (y-axis). Shower start direction calculated as unit vector pointing from shower trunk start position to shower trunk end position. Array of length nShowers
* __showerDistToVtx__: distance (in cm) between shower start point and reco neutrino vertex. Array of length nShowers
* __showerStartPos{X,Y,Z}__: X,Y,Z coordinates (in cm) of shower start position. Array of length nShowers
* __showerStartDir{X,Y,Z}__: Coordinates of unit vector pointing in initial shower direction, calculated as the vector pointing from the beginning to the end of the shower trunk. Array of length nShowers
* __showerClassified__: 1 if shower was classified, 0 otherwise. Shower is “classified” if it has at least 10 above-noise-threshold pixels in all three wire planes (noise threshold: 10). Array of length nShowers
* __showerPID__: Predicted PDG code of shower, from LArPID. Defaults to 0 if showerClassified = 0. Array of length nShowers
* __showerElScore__: Shower’s electron score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerPhScore__: Shower’s photon score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerMuScore__: Shower’s muon score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerPiScore__: Shower’s pion score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerPrScore__: Shower’s proton score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerComp__: Shower’s predicted completeness from LArPID. Defaults to -1 if showerClassified = 0. Array of length nShowers 
* __showerPurity__: Shower’s predicted purity from LArPID. Defaults to -1 if showerClassified = 0. Array of length nShowers
* __showerProcess__: Shower's predicted particle production process from LArPID (0: primary particle from neutrino interaction, 1: secondary particle with a neutral parent, 2: secondary particle with a charged parent). Defaults to -1 if showerClassified = 0. Array of length nShowers
* __showerPrimaryScore__: Shower's primary particle score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerFromNeutralScore__: Shower's secondary particle with neutral parent score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerFromChargedScore__: Shower's secondary particle with charged parent score from LArPID. Defaults to -99 if showerClassified = 0. Array of length nShowers
* __showerRecoE__: reconstructed energy of shower in MeV, calculated from calibrated collection plane charge. Defaults to -1 if showerClassified = 0. Array of length nShowers
* __showerTruePID__: The PDG code of the shower’s truth-matched particle. Defaults to 0 if no simulated particle deposited energy in shower's pixels. Not filled for data events. Array of length nShowers
* __showerTrueTID__: The trackID of the shower’s truth-matched particle. Defaults to -1 if no simulated particle deposited energy in shower's pixels. Not filled for data events. Array of length nShowers
* __showerTrueE__: The true energy (in MeV) of the shower’s truth-matched particle. Defaults to -1 if no simulated particle deposited energy in shower's pixels. Not filled for data events. Array of length nShowers
* __showerTruePurity__: The true purity of the reco shower, calculated from its truth-matched particle. Defaults to 0 if no simulated particle deposited energy in shower's pixels. Not filled for data events. Array of length nShowers
* __showerTrueComp__: The true completeness of the reco shower, calculated from its truth-matched particle. Defaults to 0 if no simulated particle deposited energy in shower's pixels. Not filled for data events. Array of length nShowers
* __showerTrueElPurity__: The fraction of the charge in the reco shower produced by simulated electrons. Not filled for data events. Array of length nShowers
* __showerTruePhPurity__: The fraction of the charge in the reco shower produced by simulated photons. Not filled for data events. Array of length nShowers
* __showerTrueMuPurity__: The fraction of the charge in the reco shower produced by simulated muons. Not filled for data events. Array of length nShowers
* __showerTruePiPurity__: The fraction of the charge in the reco shower produced by simulated pions. Not filled for data events. Array of length nShowers
* __showerTruePrPurity__: The fraction of the charge in the reco shower produced by simulated protons. Not filled for data events. Array of length nShowers


