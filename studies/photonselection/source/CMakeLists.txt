# Create the main executable

add_executable(make_photon_selection_study_tree make_photon_selection_study_tree.cxx )


# Install
install(TARGETS make_photon_selection_study_tree
        EXPORT photon_selection_study
        RUNTIME DESTINATION bin )

