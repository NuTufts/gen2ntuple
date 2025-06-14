# Source code for flash prediction executables

# Manual include directories using environment variables
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
  $ENV{LARCV_INCDIR}
  $ENV{LARLITE_COREDIR}/..
  $ENV{LARLITE_COREDIR}
  $ENV{LARFLOW_INCDIR}
  $ENV{UBLARCVAPP_INCDIR}
)

# Add executables
add_executable(calculate_flash_predictions calculate_flash_predictions.cxx)
add_executable(calculate_flash_predictions_vectorized calculate_flash_predictions_vectorized.cxx)
#add_executable(run_flashprediction run_flashprediction.cxx)

# Link libraries for calculate_flash_predictions
target_link_libraries(calculate_flash_predictions
  ${ROOT_LIBRARIES}
  -L$ENV{LARCV_LIBDIR} -lLArCVCoreBase -lLArCVCoreDataFormat -lLArCVCoreProcessor
  -L$ENV{LARLITE_LIBDIR} -llarlite_Base -llarlite_DataFormat -llarlite_LArUtil
  -L$ENV{LARFLOW_LIBDIR} -lLArFlow_Reco -lLArFlow_SCBoundary -lLArFlow_RecoUtils
  -L$ENV{UBLARCVAPP_LIBDIR} -lLArCVApp_ubdllee -lLArCVApp_UBPhotonLib -lLArCVApp_ParticleToPixelUtils
)

# Link libraries for calculate_flash_predictions_vectorized
target_link_libraries(calculate_flash_predictions_vectorized
  ${ROOT_LIBRARIES}
  -L$ENV{LARCV_LIBDIR} -lLArCVCoreBase -lLArCVCoreDataFormat -lLArCVCoreProcessor
  -L$ENV{LARLITE_LIBDIR} -llarlite_Base -llarlite_DataFormat -llarlite_LArUtil
  -L$ENV{LARFLOW_LIBDIR} -lLArFlow_Reco -lLArFlow_SCBoundary -lLArFlow_RecoUtils
  -L$ENV{UBLARCVAPP_LIBDIR} -lLArCVApp_ubdllee -lLArCVApp_UBPhotonLib -lLArCVApp_ParticleToPixelUtils
)

## Link libraries for run_flashprediction  
#target_link_libraries(run_flashprediction
#  ${ROOT_LIBRARIES}
#)

# Install executables
install(TARGETS calculate_flash_predictions calculate_flash_predictions_vectorized
        DESTINATION bin
        EXPORT flashprediction)

# Export targets
export(TARGETS calculate_flash_predictions calculate_flash_predictions_vectorized
       FILE "${CMAKE_BINARY_DIR}/flashpredictionTargets.cmake")
