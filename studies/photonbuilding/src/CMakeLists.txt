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
  ${FLASHMATCH_DATAPREP_INCLUDE_DIR}
)

# Add executables
add_executable(photonbuilding photonbuilding.cxx)


# Link libraries for photonbuilding
target_link_libraries(photonbuilding
  ${ROOT_LIBRARIES}
  -L$ENV{LARCV_LIBDIR} -lLArCVCoreBase -lLArCVCoreDataFormat -lLArCVCoreProcessor
  -L$ENV{LARLITE_LIBDIR} -llarlite_Base -llarlite_DataFormat -llarlite_LArUtil
  -L$ENV{LARFLOW_LIBDIR} -lLArFlow_Reco -lLArFlow_SCBoundary -lLArFlow_RecoUtils
  -L$ENV{UBLARCVAPP_LIBDIR} -lLArCVApp_ubdllee -lLArCVApp_UBPhotonLib -lLArCVApp_ParticleToPixelUtils -lLArCVApp_MCTools
  flashmatch::flashmatch_dataprep
)

## Link libraries for run_photonbuilding  
#target_link_libraries(run_photonbuilding
#  ${ROOT_LIBRARIES}
#)

# Install executables
install(TARGETS photonbuilding
        DESTINATION bin
        EXPORT photonbuilding)

# Export targets
export(TARGETS photonbuilding
       FILE "${CMAKE_BINARY_DIR}/photonbuildingTargets.cmake")
