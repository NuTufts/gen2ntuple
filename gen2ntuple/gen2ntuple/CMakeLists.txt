# Gen2Ntuple Core Library and Executable

# Find Python
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Include directories for UBDL dependencies
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Python3_INCLUDE_DIRS}  
  ${ROOT_INCLUDE_DIRS}
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}
  ${UBLARCVAPP_INCLUDE_DIR}
  ${LARPID_INCLUDE_DIR}
  ${LARFLOW_INCLUDE_DIR}  
  ${TORCH_INCLUDE_DIRS}
)

# Core library sources
set(GEN2NTUPLE_SOURCES
  # Data structures
  EventData.cxx
  RecoData.cxx
  BranchManager.cxx
  
  # Processing modules
  EventProcessor.cxx
  FileManager.cxx
  MCTruthProcessor.cxx
  VertexSelector.cxx
  TrackProcessor.cxx
  ShowerProcessor.cxx
  VertexSummaryProcessor.cxx
  PCAAnalyzer.cxx
  GeometryUtils.cxx
  ChargeCalculator.cxx

  # VertexSelection modules
  VertexSelectionBase.cxx
  SelectNuVertex.cxx
  HighestKPRankWithVisEnergy.cxx
  GetNuCandidateIntimeCharge.cxx
  FlashPEVertexSelection.cxx
  
  # Utilities
  WCFiducial.cxx
  Configuration.cxx
  Logger.cxx
  pionRange2T.cxx

)

set(GEN2NTUPLE_HEADERS
  # Data structures
  EventData.h
  RecoData.h
  BranchManager.h
  
  # Processing modules
  EventProcessor.h
  FileManager.h
  MCTruthProcessor.h
  VertexSelector.h
  TrackProcessor.h
  ShowerProcessor.h
  VertexSummaryProcessor.h
  PCAAnalyzer.h
  GeometryUtils.h
  ChargeCalculator.h

  # VertexSelection modules
  VertexSelectionBase.h
  SelectNuVertex.h
  HighestKPRankWithVisEnergy.h
  GetNuCandidateIntimeCharge.h
  FlashPEVertexSelection.h
  
  # Utilities
  WCFiducial.h
  Configuration.h
  Logger.h
  pionRange2T.h

)

# Create the library
add_library(gen2ntuple SHARED ${GEN2NTUPLE_SOURCES})

# Link libraries
target_link_libraries(gen2ntuple
  ${ROOT_LIBRARIES}
  -L$ENV{LARCV_LIBDIR} -lLArCVCoreBase -lLArCVCoreDataFormat -lLArCVCoreProcessor
  -L$ENV{LARLITE_LIBDIR} -llarlite_Base -llarlite_DataFormat -llarlite_LArUtil
  -L$ENV{LARPID_LIBDIR} -lLArPID_Data -lLArPID_Model -lLArPID_Interface
  -L$ENV{LARFLOW_LIBDIR} -lLArFlow_Reco -lLArFlow_SCBoundary -lLArFlow_RecoUtils LArFlow_ProngCNN
  -L$ENV{UBLARCVAPP_LIBDIR} -lLArCVApp_ubdllee -lLArCVApp_UBPhotonLib -lLArCVApp_ParticleToPixelUtils -lLArCVApp_MCTools -lLArCVApp_dbscan
  ${TORCH_LIBRARIES}
  ${Python3_LIBRARIES}
)

# Create the main executable
add_executable(make_gen2_ntuples main.cxx)
target_link_libraries(make_gen2_ntuples gen2ntuple)

# Create test executable for weight loading
add_executable(test_weight_loading ../test_weight_loading.cxx)
target_link_libraries(test_weight_loading gen2ntuple)

# Install
install(TARGETS gen2ntuple make_gen2_ntuples test_weight_loading
        EXPORT gen2ntuple
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

install(FILES ${GEN2NTUPLE_HEADERS}
        DESTINATION include/gen2ntuple)
