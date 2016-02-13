#ifndef RESPONSE_H_INCLUDED
#define RESPONSE_H_INCLUDED

#include <vector>

#include "declarations.h"

class Response : public LandauLifshitzGilbert
{
    private:
        //protected default constructor
        Response () {}

        //private variables
        time_t m_start, m_finish;
        double m_timeEst;
        std::string m_GxxString;
        std::string m_GxyString;
        std::string m_GyxString;
        std::string m_GyyString;
        PlotData m_sPlotData;
        std::vector <int> m_posInt;
        std::vector <double> m_xpos,      m_ypos,
                             m_xposPrime, m_yposPrime,
                             m_xposMod,   m_yposMod;

        Eigen::MatrixXcd  m_GxxMatrix, m_GxyMatrix,
                          m_GyxMatrix, m_GyyMatrix;

    public:
        //constructor
        Response ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision, PlotData sPlotData );

        //public functions
        void SetPositions ();
        void CalculateResponse ( int PartNum, double kxlstart, double kylstart, double kxlend, double kylend,
                                 double kMin, double kMax, int ExchField, int DipolarField, int ElectricField, int FreqGF, int kTotal );
        void WriteToFile ( int PartNum );
};

#endif // RESPONSE_H_INCLUDED
