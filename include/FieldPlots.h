#ifndef FIELDPLOTS_H_INCLUDED
#define FIELDPLOTS_H_INCLUDED

#include <vector>

#include "declarations.h"

class FieldPlots : public LandauLifshitzGilbert
{
    private:
        //protected default constructor
        FieldPlots () {}

        //protected variables
        int m_spatial_sizePrime;
        int m_numPositionsPrime;
        double m_phase, m_xpos0, m_ypos0;
        std::vector <double> m_xposFields;
        std::vector <double> m_yposFields;

    public:
        //constructor
        FieldPlots ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision );

        //public functions
        void PlotFT_Variables ();
        void Plot ();
        void PlotIntegrated ( int PartNum );
};

#endif // FIELDPLOTS_H_INCLUDED
