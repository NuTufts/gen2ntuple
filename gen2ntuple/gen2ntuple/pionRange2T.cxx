#include "pionRange2T.h"
#include "TGraph.h"

#include <array>

namespace gen2ntuple {

pionRange2T* pionRange2T::_gPionRange2T = nullptr;

pionRange2T::pionRange2T() 
{
    
    cutoff = 59.5;
    slope = 2.0699248120300755;
    intercept = 51.75526315789472;
    std::array<float, 60> range
     = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5, 51.5, 52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5, 59.5};
    std::array<float, 60 > KE 
     = { 0.5, 13.5, 20.5, 25.5, 30.5, 33.5, 38.5, 42.5, 45.5, 48.5, 51.5, 54.5, 58.5, 60.5, 62.5, 67.5, 69.5, 72.5, 75.5, 77.5, 81.5, 82.5, 86.5, 89.5, 91.5, 92.5, 94.5, 98.5, 99.5, 103.5, 106.5, 108.5, 109.5, 113.5, 115.5, 118.5, 120.5, 123.5, 126.5, 128.5, 130.5, 132.5, 132.5, 138.5, 140.5, 140.5, 145.5, 147.5, 147.5, 151.5, 151.5, 156.5, 156.5, 160.5, 163.5, 165.5, 164.5, 171.5, 170.5, 172.5};
    
    _grange2T = new TGraph(60);
    for (int i=0; i<60; i++) {
        _grange2T->SetPoint( i, range[i], KE[i] );
    } 

}

pionRange2T::~pionRange2T() 
{
    delete _grange2T;
}

double pionRange2T::Eval(double length) 
{
    if ( length>=0.5 && length < cutoff )
        return _grange2T->Eval( length );
    if ( length>=cutoff )
        return intercept + slope*length;
    return 0.;
}

pionRange2T* pionRange2T::get() 
{
    if ( _gPionRange2T==nullptr  ) {
        _gPionRange2T = new pionRange2T;
    }

    return _gPionRange2T;
}

}