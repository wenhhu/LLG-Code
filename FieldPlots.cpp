#include <cmath>
#include <iostream>

//#include <gsl/gsl_complex_math.h>

#include "constants.h"
#include "declarations.h"

using namespace std;

FieldPlots::FieldPlots ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision )
    : LandauLifshitzGilbert ( Mat_A, Mat_B, Crystal, sPrecision ),
    m_spatial_sizePrime ( 100 ),
    m_numPositionsPrime ( ( 2 * m_spatial_sizePrime + 1 ) * ( 2 * m_spatial_sizePrime + 1 ) )
    {
        m_xposFields.resize ( m_numPositionsPrime );
        m_yposFields.resize ( m_numPositionsPrime );

        for ( int j = 0; j < ( m_numPositionsPrime ); j++ )
        {
            m_xposFields[j] = 2 * ( m_Crystal.a / 2 ) * (fmod ( j, 2. * m_spatial_sizePrime + 1. ) - m_spatial_sizePrime) / m_spatial_sizePrime;
            m_yposFields[j] = 2 * ( m_Crystal.a / 2 ) * (floor ( j / (2. * m_spatial_sizePrime + 1.) ) - m_spatial_sizePrime) / m_spatial_sizePrime;
        }
    }

void FieldPlots::PlotFT_Variables ()
{
    int mode = 0;
    complex <double> M_sPlot, ExchLengthPlot, alphaPlot;

    ofstream FT_VariablesFile;
    FT_VariablesFile.open ( "Data/FT_Variables.csv", ios::out | ios::trunc );

    FT_VariablesFile << "Col. 1:   xpos / a" << endl;
    FT_VariablesFile << "Col. 2:   ypos / a" << endl;
    FT_VariablesFile << "Col. 3:   Re (M_s)" << endl;
    FT_VariablesFile << "Col. 4:   Im (M_s)" << endl;
    FT_VariablesFile << "Col. 5:   Re (ExchLength)" << endl;
    FT_VariablesFile << "Col. 6:   Im (ExchLength)" << endl;
    FT_VariablesFile << "Col. 7:   Re (alpha)" << endl;
    FT_VariablesFile << "Col. 8:   Im (alpha)" << endl << endl;

    FT_VariablesFile << "Nmax=" << m_sPrecision.Nmax << endl  <<
                        "spatial_sizePrime=" << m_spatial_sizePrime << endl <<
                        "stepsize=" << m_sPrecision.stepsize << endl <<
                        "a=" << m_Crystal.a << endl <<
                        "f=" << m_Crystal.f << endl <<
                        "mu0*H0=" << m_Crystal.mu0 * m_Crystal.H0 << endl <<
                        "Msa=" << m_Mat_A.Ms << endl <<
                        "Msb=" << m_Mat_B.Ms << endl <<
                        "alphaa=" << m_Mat_A.alpha << endl <<
                        "alphab=" << m_Mat_B.alpha << endl <<
                        "Aa=" << m_Mat_A.A << endl <<
                        "Ab=" << m_Mat_B.A << endl;

    FT_VariablesFile << endl << endl << endl;  //creates a new index for pyxplot

    FT_VariablesFile.setf ( ios::scientific );
    FT_VariablesFile.precision ( 15 );
    for ( int i = 0; i < m_numPositionsPrime; i++ )
    {
        cout << m_xposFields[i] / m_Crystal.a << " " << m_yposFields[i] / m_Crystal.a << endl;

        M_sPlot        = ComputeFields ( SATURATION_MAGNETIZATION, m_xposFields[i], m_yposFields[i], mode );
        ExchLengthPlot = ComputeFields ( EXCHANGE_LENGTH         , m_xposFields[i], m_yposFields[i], mode );
        alphaPlot      = ComputeFields ( GILBERT_DAMPING         , m_xposFields[i], m_yposFields[i], mode );

        FT_VariablesFile << m_xposFields[i] / m_Crystal.a << " " << m_yposFields[i] / m_Crystal.a << " " <<
                            M_sPlot.real()          << " " << M_sPlot.imag()          << " " <<
                            ExchLengthPlot.real()   << " " << ExchLengthPlot.imag()   << " " <<
                            alphaPlot.real()        << " " << alphaPlot.imag()        << endl;
    }
}

void FieldPlots::Plot ()
{
    int mode;
    complex <double> mx, my, Psi, hx, hy, H_x, H_y, mCurl_xWRTy, mCurl_yWRTx, hCurl_xWRTy, hCurl_yWRTx;

    ofstream FieldPlotsFile;
    FieldPlotsFile.open ( "Data/FieldPlots.csv", ios::out | ios::trunc );

    FieldPlotsFile << "Col. 1:   xpos / a" << endl;
    FieldPlotsFile << "Col. 2:   ypos / a" << endl;
    FieldPlotsFile << "Col. 3:   Re (mx)" << endl;
    FieldPlotsFile << "Col. 4:   Im (mx)" << endl;
    FieldPlotsFile << "Col. 5:   Re (my)" << endl;
    FieldPlotsFile << "Col. 6:   Im (my)" << endl;
    FieldPlotsFile << "Col. 7:   Re (Psi)" << endl;  //temporarily changed from Re (Psi) to look at demag field for slab geometry
    FieldPlotsFile << "Col. 8:   Im (Psi)" << endl;
    FieldPlotsFile << "Col. 9:   Re (hx)" << endl;
    FieldPlotsFile << "Col. 10:  Im (hx)" << endl;
    FieldPlotsFile << "Col. 11:  Re (hy)" << endl;
    FieldPlotsFile << "Col. 12:  Im (hy)" << endl;
    FieldPlotsFile << "Col. 13:  Re (H_ex_x)" << endl;
    FieldPlotsFile << "Col. 14:  Im (H_ex_x)" << endl;
    FieldPlotsFile << "Col. 15:  Re (H_ex_y)" << endl;
    FieldPlotsFile << "Col. 16:  Im (H_ex_y)" << endl;
    FieldPlotsFile << "Col. 17:  Re (mCurl_xWRTy)" << endl;
    FieldPlotsFile << "Col. 18:  Im (mCurl_xWRTy)" << endl;
    FieldPlotsFile << "Col. 19:  Re (mCurl_yWRTx)" << endl;
    FieldPlotsFile << "Col. 20:  Im (mCurl_yWRTx)" << endl;
    FieldPlotsFile << "Col. 21:  Re (hCurl_xWRTy)" << endl;
    FieldPlotsFile << "Col. 22:  Im (hCurl_xWRTy)" << endl;
    FieldPlotsFile << "Col. 23:  Re (hCurl_yWRTx)" << endl;
    FieldPlotsFile << "Col. 24:  Im (hCurl_yWRTx)" << endl;
    FieldPlotsFile << "Col. 25:  |h|/|H_ex|" << endl << endl;

    FieldPlotsFile << "Nmax=" << m_sPrecision.Nmax << endl  <<
                      "spatial_sizePrime=" << m_spatial_sizePrime << endl <<
                      "stepsize=" << m_sPrecision.stepsize << endl <<
                      "a=" << m_Crystal.a << endl <<
                      "f=" << m_Crystal.f << endl <<
                      "mu0*H0=" << m_Crystal.mu0 * m_Crystal.H0 << endl <<
                      "Msa=" << m_Mat_A.Ms << endl <<
                      "Msb=" << m_Mat_B.Ms << endl <<
                      "alphaa=" << m_Mat_A.alpha << endl <<
                      "alphab=" << m_Mat_B.alpha << endl <<
                      "Aa=" << m_Mat_A.A << endl <<
                      "Ab=" << m_Mat_B.A << endl;

    FieldPlotsFile << endl << endl << endl;  //creates a new index for pyxplot


    FieldPlotsFile.setf ( ios::scientific );
    FieldPlotsFile.precision ( 15 );
    for ( int i = 0; i < m_numPositionsPrime; i++ )
    {
        FieldPlotsFile << m_xposFields[i] / m_Crystal.a << " " << m_yposFields[i] / m_Crystal.a << " ";
        cout << m_xposFields[i] << " " << m_yposFields[i] << endl;
        for ( mode = 0; mode < 1; mode++ )
        {
            mx          = ComputeFields_kdotx ( MAGNETIZATION_x,             m_xposFields[i], m_yposFields[i], mode );
            my          = ComputeFields_kdotx ( MAGNETIZATION_y,             m_xposFields[i], m_yposFields[i], mode );
            Psi         = ComputeFields_kdotx ( MAGNETOSTATIC_POTENTIAL,     m_xposFields[i], m_yposFields[i], mode );
            hx          = ComputeFields_kdotx ( DIPOLAR_x,                   m_xposFields[i], m_yposFields[i], mode );
            hy          = ComputeFields_kdotx ( DIPOLAR_y,                   m_xposFields[i], m_yposFields[i], mode );
            H_x         = ComputeFields_kdotx ( EXCHANGE_x,                  m_xposFields[i], m_yposFields[i], mode );
            H_y         = ComputeFields_kdotx ( EXCHANGE_y,                  m_xposFields[i], m_yposFields[i], mode );
            mCurl_xWRTy = ComputeFields_kdotx ( CURL_OF_MAGNETIZATION_xWRTy, m_xposFields[i], m_yposFields[i], mode );
            mCurl_yWRTx = ComputeFields_kdotx ( CURL_OF_MAGNETIZATION_yWRTx, m_xposFields[i], m_yposFields[i], mode );
            hCurl_xWRTy = ComputeFields_kdotx ( CURL_OF_DIPOLAR_xWRTy,       m_xposFields[i], m_yposFields[i], mode );
            hCurl_yWRTx = ComputeFields_kdotx ( CURL_OF_DIPOLAR_yWRTx,       m_xposFields[i], m_yposFields[i], mode );
            FieldPlotsFile  << mx.real()          << " " << mx.imag()          << " " <<
                               my.real()          << " " << my.imag()          << " " <<
                               Psi.real()         << " " << Psi.imag()         << " " <<
                               hx.real()          << " " << hx.imag()          << " " <<
                               hy.real()          << " " << hy.imag()          << " " <<
                               H_x.real()         << " " << H_x.imag()         << " " <<
                               H_y.real()         << " " << H_y.imag()         << " " <<
                               mCurl_xWRTy.real() << " " << mCurl_xWRTy.imag() << " " <<
                               mCurl_yWRTx.real() << " " << mCurl_yWRTx.imag() << " " <<
                               hCurl_xWRTy.real() << " " << hCurl_xWRTy.imag() << " " <<
                               hCurl_yWRTx.real() << " " << hCurl_yWRTx.imag() << " " <<
                               ( norm ( hx )  + norm ( hy ) ) /
                               ( norm ( H_x ) + norm ( H_y ) );
        }
        FieldPlotsFile << endl;
    }
}

void FieldPlots::PlotIntegrated ( int PartNum )
{
    int mode;
    int const NumOfModes = 1;
    double MagVecRef = 0;
    double h_x_Integral         [ NumOfModes ] = { 0. }, H_ex_x_Integral         [ NumOfModes ] = { 0. },
           h_y_Integral         [ NumOfModes ] = { 0. }, H_ex_y_Integral         [ NumOfModes ] = { 0. },
           h_Magnitude_Integral [ NumOfModes ] = { 0. }, H_ex_Magnitude_Integral [ NumOfModes ] = { 0. },
           h_Over_H_ex_Integral [ NumOfModes ] = { 0. },
           MagVecSymm_Integral  [ NumOfModes ] = { 0. }, MagVecSymm_Max          [ NumOfModes ] = { 0. };
    complex <double> mx, my, mx0, my0, hx, hy, H_x, H_y;

    ofstream Integrated_Fields;
    Integrated_Fields.open ( NameFILE ( "Integrated_Fields/Integrated_Fields_a_Nmax", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".csv" ).c_str(), ios::out | ios::trunc );

    for ( int i = 0; i < m_sPrecision.numPositionsPrime; i++ )
    {
        cout << m_xposFields[i] << " " << m_yposFields[i] << endl;
        for ( mode = 0; mode < NumOfModes; mode++ )
        {
            hx  = ComputeFields_kdotx ( DIPOLAR_x,  m_xposFields[i], m_yposFields[i], mode );
            hy  = ComputeFields_kdotx ( DIPOLAR_y,  m_xposFields[i], m_yposFields[i], mode );
            H_x = ComputeFields_kdotx ( EXCHANGE_x, m_xposFields[i], m_yposFields[i], mode );
            H_y = ComputeFields_kdotx ( EXCHANGE_y, m_xposFields[i], m_yposFields[i], mode );

            h_x_Integral            [ mode ] += abs ( hx );
            h_y_Integral            [ mode ] += abs ( hy );
            H_ex_x_Integral         [ mode ] += abs ( H_x );
            H_ex_y_Integral         [ mode ] += abs ( H_y );
            h_Magnitude_Integral    [ mode ] += sqrt ( norm ( hx  ) + norm ( hy  ) );
            H_ex_Magnitude_Integral [ mode ] += sqrt ( norm ( H_x ) + norm ( H_y ) );
            h_Over_H_ex_Integral    [ mode ] += sqrt ( norm ( hx  ) + norm ( hy  ) ) /
                                                sqrt ( norm ( H_x ) + norm ( H_y ) );

            mx0 = ComputeFields_kdotx ( MAGNETIZATION_x,  m_xposFields[i], m_yposFields[i], mode );
            my0 = ComputeFields_kdotx ( MAGNETIZATION_y,  m_xposFields[i], m_yposFields[i], mode );
            mx  = ComputeFields_kdotx ( MAGNETIZATION_x, -m_xposFields[i], m_yposFields[i], mode );
            my  = ComputeFields_kdotx ( MAGNETIZATION_y, -m_xposFields[i], m_yposFields[i], mode );

            MagVecRef =      abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                   sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                   sqrt ( norm ( mx0 ) + norm ( my0 ) );


            if ( abs ( mx0 ) < 1e-10 )
            {
                mx = mx0;
                my = my0;
            }

            if ( MagVecRef > MagVecSymm_Max [ mode ] ) MagVecSymm_Max [ mode ] = MagVecRef;

            MagVecSymm_Integral [ mode ] += MagVecRef;

//            if (  m_xposFields[i] < 0 && m_yposFields[i] < m_xposFields[i] + m_Crystal.a / 200 )  //find max value of magnetization asymmetry
//            {
//                mx0 = ComputeFields_kdotx ( MAGNETIZATION_x,  m_xposFields[i], m_yposFields[i], mode );
//                my0 = ComputeFields_kdotx ( MAGNETIZATION_y,  m_xposFields[i], m_yposFields[i], mode );
//                mx  = ComputeFields_kdotx ( MAGNETIZATION_x, -m_xposFields[i], m_yposFields[i], mode );
//                my  = ComputeFields_kdotx ( MAGNETIZATION_y, -m_xposFields[i], m_yposFields[i], mode );
//
//                Max_MagVecSymm_Comp = abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
//                                            sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
//                                            sqrt ( norm ( mx0 ) + norm ( my0 ) );
//                if ( Max_MagVecSymm [ mode ] < Max_MagVecSymm_Comp )
//                {
//                    Max_MagVecSymm [ mode ] = Max_MagVecSymm_Comp;
//                }
//            }
//            Max_MagVecSymm_Comp = 0;
        }
    }

    if ( PartNum == 0 )
    {
        Integrated_Fields << "Col. 1:   Nmax" << endl;
        Integrated_Fields << "Col. 2:   ssp" << endl;
        Integrated_Fields << "Col. 3:   f" << endl;
        Integrated_Fields << "Col. 4:   a (nm)" << endl;
        Integrated_Fields << "Col. 5:   Int |h_x|" << endl;
        Integrated_Fields << "Col. 6:   Int |h_y|" << endl;
        Integrated_Fields << "Col. 7:   Int |H_ex_x|" << endl;
        Integrated_Fields << "Col. 8:   Int |H_ex_y|" << endl;
        Integrated_Fields << "Col. 9:   Int |h|" << endl;
        Integrated_Fields << "Col. 10:  Int |H_ex|" << endl;
        Integrated_Fields << "Col. 11:  Int |h|/|H_ex|" << endl;
        Integrated_Fields << "Col. 12:  Int |MagVecSymm|" << endl;
        Integrated_Fields << "Col. 13:  Max |MagVecSymm|" << endl;
        Integrated_Fields << "Col. 14:  Im (omega)" << endl;
        Integrated_Fields << "Col. 15:  Re (omega)" << endl << endl;

        Integrated_Fields << "Nmax=" << m_sPrecision.Nmax << endl  <<
                             "spatial_sizePrime=" << m_sPrecision.spatial_sizePrime << endl <<
                             "stepsize=" << m_sPrecision.stepsize << endl <<
                             "a=" << m_Crystal.a << endl <<
                             "f=" << m_Crystal.f << endl <<
                             "mu0*H0=" << m_Crystal.mu0 * m_Crystal.H0 << endl <<
                             "Msa=" << m_Mat_A.Ms << endl <<
                             "Msb=" << m_Mat_B.Ms << endl <<
                             "alphaa=" << m_Mat_A.alpha << endl <<
                             "alphab=" << m_Mat_B.alpha << endl <<
                             "Aa=" << m_Mat_A.A << endl <<
                             "Ab=" << m_Mat_B.A << endl;

        Integrated_Fields << endl << endl << endl;  //creates a new index for pyxplot
    }
    Integrated_Fields.setf ( ios::scientific );
    Integrated_Fields.precision ( 15 );

    Integrated_Fields << m_sPrecision.Nmax << " " <<
                         m_sPrecision.spatial_sizePrime << " " <<
                         m_Crystal.f << " " <<
                         m_Crystal.a;
    for ( mode = 0;  mode < NumOfModes; mode++ )
    {
        Integrated_Fields << " " << h_x_Integral            [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    h_y_Integral            [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    H_ex_x_Integral         [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    H_ex_y_Integral         [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    h_Magnitude_Integral    [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    H_ex_Magnitude_Integral [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    h_Over_H_ex_Integral    [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    MagVecSymm_Integral     [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    MagVecSymm_Max          [ mode ] * pow ( m_Crystal.a / ( 2 * m_sPrecision.spatial_sizePrime + 1 ), 2 ) << " " <<
                                    m_evals ( mode ).real() << " " <<
                                    m_evals ( mode ).real();
    }
    Integrated_Fields << endl;
}
