#include "WCFiducial.h"
#include <cmath>

namespace gen2ntuple {

  WCFiducial* WCFiducial::_gWCFiducial = nullptr;

  WCFiducial* WCFiducial::getME( double boundary_dis_cut_cm ) {

    if ( _gWCFiducial==nullptr ) {
      if ( boundary_dis_cut_cm<0 )
        _gWCFiducial = new WCFiducial;
      else
        _gWCFiducial = new WCFiducial( boundary_dis_cut_cm );
    }
    else {
      if ( boundary_dis_cut_cm>=0.0 ) {
        //  check boundary distance and refresh if necessary
    
        double dis_diff = std::fabs( boundary_dis_cut_cm-_gWCFiducial->getCurrentBoundary() );
        if ( dis_diff>1.e-4 ) {
          delete _gWCFiducial;
          _gWCFiducial = new WCFiducial( boundary_dis_cut_cm );
        }
      }
    }

    return _gWCFiducial;

  }

}


