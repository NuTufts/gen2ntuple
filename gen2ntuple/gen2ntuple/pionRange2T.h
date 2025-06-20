#ifndef __GEN2NTUPLE_PIONRANGE2T_H__
#define __GEN2NTUPLE_PIONRANGE2T_H__

class TGraph;

namespace gen2ntuple {

    class pionRange2T {

    private:

        pionRange2T();
        ~pionRange2T();


    public:

        static pionRange2T* get();

        double Eval(double length);


    private:

        double cutoff;
        double slope;
        double intercept;
        TGraph* _grange2T;

        static pionRange2T* _gPionRange2T;

    };

}

#endif