#include "MCTruthProcessor.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <Python.h>

// LArLite includes
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/mctruth.h"
#include "larlite/DataFormat/mcpart.h"
#include "larlite/DataFormat/mcnu.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctrack.h"

// LArCV includes
#include "larcv/core/DataFormat/IOManager.h"

#include "ublarcvapp/MCTools/NeutrinoVertex.h"

#include "WCFiducial.h"

namespace gen2ntuple {

MCTruthProcessor::MCTruthProcessor()
: pSCE(nullptr) 
{
    initializePython();
    pSCE = new larutil::SpaceChargeMicroBooNE;
}

MCTruthProcessor::~MCTruthProcessor() {
    finalizePython();
    delete pSCE;
    pSCE = nullptr;
}

bool MCTruthProcessor::processEvent(larlite::storage_manager* larlite_io, 
                                    larcv::IOManager* larcv_io,
                                    EventData* event_data) {
    
    if (!event_data) {
        std::cerr << "MCTruthProcessor: EventData pointer is null" << std::endl;
        return false;
    }
    
    // Set cross-section weight
    event_data->xsecWeight = getEventWeight(event_data->run, event_data->subrun, event_data->event);
    
    // Extract neutrino truth information
    if (!extractNeutrinoTruth(larlite_io, event_data)) {
        return false;
    }
    
    // Extract true vertex
    if (!extractTrueVertex(larlite_io, event_data)) {
        return false;
    }
    
    // Extract primary particles from generator
    if (!extractPrimaryParticles(larlite_io, event_data)) {
        return false;
    }
    
    // Extract simulated particles from Geant4
    if (!extractSimulatedParticles(larlite_io, event_data)) {
        return false;
    }
    
    return true;
}

bool MCTruthProcessor::extractNeutrinoTruth(larlite::storage_manager* larlite_io, EventData* event_data) {
    
    // Get MC truth from generator
    auto ev_mctruth = larlite_io->get_data<larlite::event_mctruth>("generator");
    if (!ev_mctruth || ev_mctruth->size() == 0) {
        // Set default values
        event_data->trueNuE = -1.0f;
        event_data->trueNuPDG = 0;
        event_data->trueNuCCNC = -1;
        event_data->trueNuMode = -1;
        event_data->trueNuIntrxnType = -1;
        event_data->trueLepE = -1.0f;
        event_data->trueLepPDG = 0;
        return true;
    }
    
    const auto& mctruth = ev_mctruth->at(0);
    
    if (!isNeutrinoInteraction(mctruth)) {
        // Not a neutrino interaction
        event_data->trueNuE = -1.0f;
        event_data->trueNuPDG = 0;
        event_data->trueNuCCNC = -1;
        event_data->trueNuMode = -1;
        event_data->trueNuIntrxnType = -1;
        event_data->trueLepE = -1.0f;
        event_data->trueLepPDG = 0;
        return true;
    }
    
    // Extract neutrino information
    const auto& mcnu = mctruth.GetNeutrino();
    
    event_data->trueNuE = mcnu.Nu().Momentum().E(); // Keep GeV for now // * 1000.0f; // Convert to MeV
    event_data->trueNuPDG = mcnu.Nu().PdgCode();
    event_data->trueNuCCNC = mcnu.CCNC(); // CC=0, NC=1
    event_data->trueNuMode = mcnu.Mode();
    event_data->trueNuIntrxnType = mcnu.InteractionType();
    
    // Extract lepton information
    if ( mcnu.CCNC()==0 ) {
        event_data->trueLepE = mcnu.Lepton().Momentum().E(); // Keep GeV for now // * 1000.0f; // Convert to MeV
        event_data->trueLepPDG = mcnu.Lepton().PdgCode();
    }
    else {
        event_data->trueLepPDG = 0;
        event_data->trueLepE   = -9.0;
    }
    
    return true;
}

bool MCTruthProcessor::extractTrueVertex(larlite::storage_manager* larlite_io, EventData* event_data) {
    
    auto ev_mctruth = larlite_io->get_data<larlite::event_mctruth>("generator");
    if (!ev_mctruth || ev_mctruth->size() == 0) {
        event_data->trueVtxX = event_data->trueVtxY = event_data->trueVtxZ = -999.0f;
        return true;
    }
    
    const auto& mctruth = ev_mctruth->at(0);
    
    if (!isNeutrinoInteraction(mctruth)) {
        event_data->trueVtxX = event_data->trueVtxY = event_data->trueVtxZ = -999.0f;
        return true;
    }
    
    const auto& mcnu = mctruth.GetNeutrino();
    
    // Note: These are the raw vertex positions - SCE correction would be applied elsewhere
    ublarcvapp::mctools::NeutrinoVertex mcNuVertexer;
    auto mcNuVertex = mcNuVertexer.getPos3DwSCE( *larlite_io, pSCE );

    event_data->trueVtxX = mcNuVertex[0];
    event_data->trueVtxY = mcNuVertex[1];
    event_data->trueVtxZ = mcNuVertex[2];
    
    return true;
}

bool MCTruthProcessor::extractPrimaryParticles(larlite::storage_manager* larlite_io, EventData* event_data) {
    
    auto ev_mctruth = larlite_io->get_data<larlite::event_mctruth>("generator");
    if (!ev_mctruth || ev_mctruth->size() == 0) {
        event_data->nTruePrimParts = 0;
        return true;
    }
    
    const auto& mctruth = ev_mctruth->at(0);
    
    int n_primary = 0;
    
    for (int i = 0; i < mctruth.NParticles(); i++) {
        const auto& mcpart = mctruth.GetParticle(i);
        
        if (isPrimaryParticle(mcpart)) {

            if ( n_primary< EventData::MAX_TRUE_PRIMARY_PARTS ) {

                event_data->truePrimPartPDG[n_primary] = mcpart.PdgCode();

                event_data->truePrimPartE[n_primary]   = mcpart.Momentum().E();  // * 1000.0f; // Convert to MeV
                event_data->truePrimPartPx[n_primary]  = mcpart.Momentum().Px(); // * 1000.0f;
                event_data->truePrimPartPy[n_primary]  = mcpart.Momentum().Py(); // * 1000.0f;
                event_data->truePrimPartPz[n_primary]  = mcpart.Momentum().Pz(); // * 1000.0f;

                // Get Start Point: Applying Space Charge Effect
                auto startpos = mcpart.Position(0);
                // dont i need to check if inside?
                auto corrected = getSCECorrectedPos(startpos.X(),startpos.Y(),startpos.Z());
                event_data->truePrimPartX[n_primary] = corrected.X();
                event_data->truePrimPartY[n_primary] = corrected.Y();
                event_data->truePrimPartZ[n_primary] = corrected.Z();

                // Determine Endpoint Containment
                auto endpos = mcpart.Position(mcpart.Trajectory().size()-1);
                TVector3 corrected_end
                    = getSCECorrectedPos( endpos.X(), endpos.Y(), endpos.Z() );
                bool isContained = WCFiducial::getME()->insideFV( corrected_end.X(), 
                                                                  corrected_end.Y(),
                                                                  corrected_end.Z() );
                event_data->truePrimPartContained[n_primary] = (isContained) ? 1 : 0;

                n_primary++;
            }
        }
    }
    
    event_data->nTruePrimParts = n_primary;
    
    return true;
}

bool MCTruthProcessor::extractSimulatedParticles(larlite::storage_manager* larlite_io, EventData* event_data) {
    
    // Get MCParticle information from Geant4
    auto mctracks  = larlite_io->get_data<larlite::event_mctrack>("mcreco");
    auto mcshowers = larlite_io->get_data<larlite::event_mcshower>("mcreco");
    auto mcshower_detectable = larlite_io->get_data<larlite::event_mcshower>("mcdetectableshower");

    if ( !mctracks && !mcshowers ) {
        event_data->nTrueSimParts = 0;
        return true;
    }
    
    int n_sim = 0;
    
    // track loop first
    for (size_t i = 0; i < mctracks->size() && n_sim < EventData::MAX_TRUE_SIM_PARTS; i++) {
        const auto& mcpart = mctracks->at(i);
        
        // Store all particles above energy threshold
        if (mcpart.Start().E() > 0.001) { // 1 MeV threshold

            event_data->trueSimPartPDG[n_sim] = mcpart.PdgCode();
            event_data->trueSimPartTID[n_sim] = mcpart.TrackID();
            event_data->trueSimPartMID[n_sim] = mcpart.MotherTrackID();

            if ( mcpart.Process()=="primary") {
                event_data->trueSimPartProcess[n_sim] = 0;
            }
            else if ( mcpart.Process()=="Decay" ) {
                event_data->trueSimPartProcess[n_sim] = 1;
            }
            else {
                event_data->trueSimPartProcess[n_sim] = 2;
            }
            auto startpos = mcpart.Start();
            auto endpos   = mcpart.End();
            TVector3 startsce = getSCECorrectedPos( startpos.X(), startpos.Y(), startpos.Z() );
            TVector3 endsce   = getSCECorrectedPos( endpos.X(),   endpos.Y(),   endpos.Z() );

            // Start position and first Edep position the same
            event_data->trueSimPartX[n_sim] = startsce.X();
            event_data->trueSimPartY[n_sim] = startsce.Y();
            event_data->trueSimPartZ[n_sim] = startsce.Z();
            event_data->trueSimPartEDepX[n_sim] = startsce.X();
            event_data->trueSimPartEDepY[n_sim] = startsce.Y();
            event_data->trueSimPartEDepZ[n_sim] = startsce.Z();
            
            event_data->trueSimPartE[n_sim]   = mcpart.Start().E();  // * 1000.0f; // Convert to MeV
            event_data->trueSimPartPx[n_sim]  = mcpart.Start().Px(); // * 1000.0f;
            event_data->trueSimPartPy[n_sim]  = mcpart.Start().Py(); // * 1000.0f;
            event_data->trueSimPartPz[n_sim]  = mcpart.Start().Pz(); // * 1000.0f;
            
            // End position (trajectory endpoint) 
            event_data->trueSimPartEndX[n_sim] = endsce.X();
            event_data->trueSimPartEndY[n_sim] = endsce.Y();
            event_data->trueSimPartEndZ[n_sim] = endsce.Z();
            event_data->trueSimPartEndE[n_sim]   = mcpart.End().E();  // * 1000.0f; // Convert to MeV
            event_data->trueSimPartEndPx[n_sim]  = mcpart.End().Px(); // * 1000.0f;
            event_data->trueSimPartEndPy[n_sim]  = mcpart.End().Py(); // * 1000.0f;
            event_data->trueSimPartEndPz[n_sim]  = mcpart.End().Pz(); // * 1000.0f;

            // End position containment test
            bool isContained 
              = WCFiducial::getME()->insideFV( endsce.X(), endsce.Y(), endsce.Z() );

            event_data->trueSimPartContained[n_sim] = (isContained) ? 1 : 0;
            
            n_sim++;
        }
    }

    // shower loop
    for (size_t i = 0; i < mcshowers->size() && n_sim < EventData::MAX_TRUE_SIM_PARTS; i++) {
        const auto& mcpart = mcshowers->at(i);
        
        // Store all particles above energy threshold
        if (mcpart.Start().E() > 0.001) { // 1 MeV threshold

            event_data->trueSimPartPDG[n_sim] = mcpart.PdgCode();
            event_data->trueSimPartTID[n_sim] = mcpart.TrackID();
            event_data->trueSimPartMID[n_sim] = mcpart.MotherTrackID();

            if ( mcpart.Process()=="primary") {
                event_data->trueSimPartProcess[n_sim] = 0;
            }
            else if ( mcpart.Process()=="Decay" ) {
                event_data->trueSimPartProcess[n_sim] = 1;
            }
            else {
                event_data->trueSimPartProcess[n_sim] = 2;
            }
            auto startpos = mcpart.Start();
            auto endpos   = mcpart.End();
            TVector3 startsce = getSCECorrectedPos( startpos.X(), startpos.Y(), startpos.Z() );
            TVector3 endsce   = getSCECorrectedPos( endpos.X(),   endpos.Y(),   endpos.Z() );

            // Start position and first Edep position the same
            event_data->trueSimPartX[n_sim] = startsce.X();
            event_data->trueSimPartY[n_sim] = startsce.Y();
            event_data->trueSimPartZ[n_sim] = startsce.Z();
            event_data->trueSimPartEDepX[n_sim] = startsce.X();
            event_data->trueSimPartEDepY[n_sim] = startsce.Y();
            event_data->trueSimPartEDepZ[n_sim] = startsce.Z();
            
            event_data->trueSimPartE[n_sim]   = mcpart.Start().E();  // * 1000.0f; // Convert to MeV
            event_data->trueSimPartPx[n_sim]  = mcpart.Start().Px(); // * 1000.0f;
            event_data->trueSimPartPy[n_sim]  = mcpart.Start().Py(); // * 1000.0f;
            event_data->trueSimPartPz[n_sim]  = mcpart.Start().Pz(); // * 1000.0f;
            
            // End position (trajectory endpoint) 
            event_data->trueSimPartEndX[n_sim] = endsce.X();
            event_data->trueSimPartEndY[n_sim] = endsce.Y();
            event_data->trueSimPartEndZ[n_sim] = endsce.Z();
            event_data->trueSimPartEndE[n_sim]   = mcpart.End().E();  // * 1000.0f; // Convert to MeV
            event_data->trueSimPartEndPx[n_sim]  = mcpart.End().Px(); // * 1000.0f;
            event_data->trueSimPartEndPy[n_sim]  = mcpart.End().Py(); // * 1000.0f;
            event_data->trueSimPartEndPz[n_sim]  = mcpart.End().Pz(); // * 1000.0f;

            // End position containment test
            bool isContained 
              = WCFiducial::getME()->insideFV( endsce.X(), endsce.Y(), endsce.Z() );

            event_data->trueSimPartContained[n_sim] = (isContained) ? 1 : 0;
            
            n_sim++;
        }
    }
    
    event_data->nTrueSimParts = n_sim;
    
    return true;
}

bool MCTruthProcessor::isNeutrinoInteraction(const larlite::mctruth& mctruth) const {
    return mctruth.Origin() == 1; // kBeamNeutrino
}

bool MCTruthProcessor::isPrimaryParticle(const larlite::mcpart& mcpart) const {
    // Primary particles have status code 1 in GENIE
    return mcpart.StatusCode() == 1;
}

int MCTruthProcessor::getParticleProcess(const larlite::mcpart& mcpart) const {
    // Simplified process identification
    // In practice, would use mcpart.Process() and map to specific codes
    return 0; // Placeholder
}

bool MCTruthProcessor::loadWeights() {
    if (weight_file_ == "none" || weight_file_.empty()) {
        return true;
    }
    
    if (!python_initialized_) {
        std::cerr << "MCTruthProcessor: Python not initialized, cannot load weights" << std::endl;
        return false;
    }
    
    // Clear existing weights
    event_weights_.clear();
    
    // Python code to load the pickle file
    PyObject* pModule = nullptr;
    PyObject* pPickle = nullptr;
    PyObject* pFile = nullptr;
    PyObject* pDict = nullptr;
    
    try {
        // Import pickle module
        pPickle = PyImport_ImportModule("pickle");
        if (!pPickle) {
            PyErr_Print();
            throw std::runtime_error("Failed to import pickle module");
        }
        
        // Open the file in binary read mode
        PyObject* pBuiltins = PyEval_GetBuiltins();
        PyObject* pOpen = PyDict_GetItemString(pBuiltins, "open");
        PyObject* pArgs = PyTuple_Pack(2, PyUnicode_FromString(weight_file_.c_str()), 
                                          PyUnicode_FromString("rb"));
        pFile = PyObject_CallObject(pOpen, pArgs);
        Py_DECREF(pArgs);
        
        if (!pFile) {
            PyErr_Print();
            throw std::runtime_error("Failed to open weight file: " + weight_file_);
        }
        
        // Load the pickle data
        PyObject* pLoad = PyObject_GetAttrString(pPickle, "load");
        pArgs = PyTuple_Pack(1, pFile);
        pDict = PyObject_CallObject(pLoad, pArgs);
        Py_DECREF(pArgs);
        Py_DECREF(pLoad);
        
        if (!pDict) {
            PyErr_Print();
            throw std::runtime_error("Failed to load pickle data");
        }
        
        // Close the file
        PyObject* pClose = PyObject_GetAttrString(pFile, "close");
        PyObject_CallObject(pClose, nullptr);
        Py_DECREF(pClose);
        
        // Parse the nested dictionary structure: dict[run][subrun][event] = weight
        PyObject *key1, *value1;
        Py_ssize_t pos1 = 0;
        
        while (PyDict_Next(pDict, &pos1, &key1, &value1)) {
            int run = PyLong_AsLong(key1);
            
            if (PyDict_Check(value1)) {
                PyObject *key2, *value2;
                Py_ssize_t pos2 = 0;
                
                while (PyDict_Next(value1, &pos2, &key2, &value2)) {
                    int subrun = PyLong_AsLong(key2);
                    
                    if (PyDict_Check(value2)) {
                        PyObject *key3, *value3;
                        Py_ssize_t pos3 = 0;
                        
                        while (PyDict_Next(value2, &pos3, &key3, &value3)) {
                            int event = PyLong_AsLong(key3);
                            float weight = (float)PyFloat_AsDouble(value3);
                            
                            auto key = std::make_tuple(run, subrun, event);
                            event_weights_[key] = weight;
                        }
                    }
                }
            }
        }
        
        std::cout << "MCTruthProcessor: Loaded " << event_weights_.size() 
                  << " event weights from " << weight_file_ << std::endl;
        
        // Cleanup
        Py_XDECREF(pDict);
        Py_XDECREF(pFile);
        Py_XDECREF(pPickle);
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "MCTruthProcessor: Error loading weights: " << e.what() << std::endl;
        
        // Cleanup on error
        Py_XDECREF(pDict);
        Py_XDECREF(pFile);
        Py_XDECREF(pPickle);
        
        return false;
    }
}

float MCTruthProcessor::getEventWeight(int run, int subrun, int event) const {
    auto key = std::make_tuple(run, subrun, event);
    
    auto it = event_weights_.find(key);
    if (it != event_weights_.end()) {
        return it->second;
    }
    
    // Default weight
    return 1.0f;
}

bool MCTruthProcessor::initializePython() {
    if (python_initialized_) {
        return true;
    }
    
    // Check if Python is already initialized (e.g., by another component)
    if (!Py_IsInitialized()) {
        Py_Initialize();
        
        // Add current directory to Python path to ensure we can import modules
        PyObject* sysPath = PySys_GetObject("path");
        PyList_Append(sysPath, PyUnicode_FromString("."));
    }
    
    python_initialized_ = true;
    return true;
}

void MCTruthProcessor::finalizePython() {
    // Note: We don't call Py_Finalize() here because Python might be
    // used by other components. In a real application, Python finalization
    // should be handled at the application level, not component level.
    python_initialized_ = false;
}

TVector3 MCTruthProcessor::getSCECorrectedPos( double x, double y, double z ) {
    std::vector<double> offset
        = pSCE->GetPosOffsets(x, y, z);
    TVector3 corrected( x-offset[0]+0.7, 
                        y+offset[1], 
                        z+offset[2]);
    return corrected;
}

} // namespace gen2ntuple