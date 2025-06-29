# Create the main executable

# Include directories for UBDL dependencies
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${ROOT_INCLUDE_DIRS}
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}
  ${UBLARCVAPP_INCLUDE_DIR}
  ${LARPID_INCLUDE_DIR}
  ${LARFLOW_INCLUDE_DIR}  
  ${TORCH_INCLUDE_DIRS}
)

add_executable(make_photon_selection_study_tree make_photon_selection_study_tree.cxx )

# Link libraries
target_link_libraries( make_photon_selection_study_tree PUBLIC
  ${ROOT_LIBRARIES}
  LArCVCoreBase LArCVCoreDataFormat LArCVCoreProcessor LArCVApp_UBImageMod
  larlite_Base  larlite_DataFormat  larlite_LArUtil
  LArPID_Data LArPID_Model LArPID_Interface
  LArFlow_Reco LArFlow_SCBoundary LArFlow_RecoUtils
  LArCVApp_ubdllee LArCVApp_UBPhotonLib LArCVApp_ParticleToPixelUtils LArCVApp_MCTools LArCVApp_dbscan
  ${TORCH_LIBRARIES}
#  ${Python3_LIBRARIES}
)


# Install
install(TARGETS make_photon_selection_study_tree
        EXPORT photon_selection_study
        RUNTIME DESTINATION bin )

