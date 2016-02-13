#include <cmath>
#include <cassert>
#include <iostream>

//#include <gsl/gsl_complex_math.h>

#include "constants.h"
#include "declarations.h"

using namespace std;
using namespace Eigen;

Response::Response ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision, PlotData sPlotData )
    : LandauLifshitzGilbert ( Mat_A, Mat_B, Crystal, sPrecision ),
    m_sPlotData ( sPlotData )
    {
        m_posInt.resize    ( m_sPlotData.numPositions, 0.0 );
        m_xpos.resize      ( m_sPlotData.numPositions, 0.0 );
        m_ypos.resize      ( m_sPlotData.numPositions, 0.0 );
        m_xposPrime.resize ( m_sPrecision.numPositionsPrime, 0.0 );
        m_yposPrime.resize ( m_sPrecision.numPositionsPrime, 0.0 );
        m_xposMod.resize   ( m_sPlotData.numPositions, 0.0 );
        m_yposMod.resize   ( m_sPlotData.numPositions, 0.0 );

        m_GxxMatrix.resize ( m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
        m_GxyMatrix.resize ( m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
        m_GyxMatrix.resize ( m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
        m_GyyMatrix.resize ( m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );

        m_GxxMatrix.setZero();
        m_GxyMatrix.setZero();
        m_GyxMatrix.setZero();
        m_GyyMatrix.setZero();
    }

void Response::SetPositions()
{
    //Create array for spatial grid.
    for ( int j = 0; j < ( m_sPlotData.numPositions ); j++ )
    {
        m_posInt[j] = -1;
    }

    //Create square grid in which Response function will be evaluated
    for ( int j = 0; j < ( m_sPlotData.numPositions ); j++ )
    {
        m_xpos[j] = m_sPlotData.xmax * m_Crystal.a * (fmod ( j, 2. * m_sPlotData.spatial_size + 1. ) - m_sPlotData.spatial_size) / m_sPlotData.spatial_size;
        m_ypos[j] = m_sPlotData.ymax * m_Crystal.a * (floor ( j / (2. * m_sPlotData.spatial_size + 1.) ) - m_sPlotData.spatial_size) / m_sPlotData.spatial_size;
    }

    //Create square grid within 1st unit cell. If source amplitude is not significantly small outside this area, then increase area size.
    for ( int j = 0; j < ( m_sPrecision.numPositionsPrime ); j++ )
    {
        m_xposPrime[j] = (m_Crystal.a/2) * (fmod ( j, 2. * m_sPrecision.spatial_sizePrime + 1. ) - m_sPrecision.spatial_sizePrime) / m_sPrecision.spatial_sizePrime;
        m_yposPrime[j] = (m_Crystal.a/2) * (floor ( j / (2. * m_sPrecision.spatial_sizePrime + 1.) ) - m_sPrecision.spatial_sizePrime) / m_sPrecision.spatial_sizePrime;
    }

    //reduce xpos/ypos to unit cell (xposprime/yposprime)
    for ( int j = 0; j < ( m_sPlotData.numPositions ); j++ )
    {
        if ( lattice == 1 )  //square lattice
        {
            if ( m_xpos[j] < 0 ) m_xposMod[j] = fmod ( m_xpos[j], -m_Crystal.a );
            else if ( m_xpos[j] > 0 ) m_xposMod[j] = fmod ( m_xpos[j], m_Crystal.a );
            else m_xposMod[j] = m_xpos[j];
            if ( m_xposMod[j] > m_Crystal.a / 2 ) m_xposMod[j] = m_xposMod[j] - m_Crystal.a;
            else if ( m_xposMod[j] < -m_Crystal.a / 2 ) m_xposMod[j] = m_xposMod[j] + m_Crystal.a;

            if ( m_ypos[j] < 0 ) m_yposMod[j] = fmod ( m_ypos[j], -m_Crystal.a );
            else if ( m_ypos[j] > 0 ) m_yposMod[j] = fmod ( m_ypos[j], m_Crystal.a );
            else m_yposMod[j] = m_ypos[j];
            if ( m_yposMod[j] > m_Crystal.a / 2 ) m_yposMod[j] = m_yposMod[j] - m_Crystal.a;
            else if ( m_yposMod[j] < -m_Crystal.a / 2 ) m_yposMod[j] = m_yposMod[j] + m_Crystal.a;
        }
        else if ( lattice == 2 )  //hexagonal lattice
        {
            // Remap positions to -3a/sqrt(3) < m_xpos < 3a/sqrt(3) and -a < m_ypos < a
            if ( m_ypos[j] < 0 ) m_yposMod[j] = fmod ( m_ypos[j], -3 * m_Crystal.a / sqrt(3) );
            else if ( m_ypos[j] > 0 ) m_yposMod[j] = fmod ( m_ypos[j], 3 * m_Crystal.a / sqrt(3) );
            else m_yposMod[j] = m_ypos[j];
            if ( m_xpos[j] < 0 ) m_xposMod[j] = fmod ( m_xpos[j], -m_Crystal.a );
            else if ( m_xpos[j] > 0 ) m_xposMod[j] = fmod ( m_xpos[j], m_Crystal.a );
            else m_xposMod[j] = m_xpos[j];

            //Map x & y positions to 1st unit cell
            //center x region
            if ( ( abs ( m_yposMod[j] ) < m_Crystal.a / ( 2 * sqrt(3) ) || abs ( m_yposMod[j] ) > 5 * m_Crystal.a / ( 2 * sqrt(3) ) ) && m_xposMod[j] < -m_Crystal.a / 2 ) m_xposMod[j] = m_xposMod[j] + m_Crystal.a;
            else if ( ( abs ( m_yposMod[j] ) < m_Crystal.a / ( 2 * sqrt(3) ) || abs ( m_yposMod[j] ) > 5 * m_Crystal.a / ( 2 * sqrt(3) ) ) && m_xposMod[j] > m_Crystal.a / 2 ) m_xposMod[j] = m_xposMod[j] - m_Crystal.a;
            if ( m_yposMod[j] < -5 * m_Crystal.a / ( 2 * sqrt(3) ) && abs ( m_xposMod[j] ) <= m_Crystal.a / 2 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / sqrt(3);
            else if ( m_yposMod[j] > 5 * m_Crystal.a / ( 2 * sqrt(3) ) && abs ( m_xposMod[j] ) <= m_Crystal.a / 2 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / sqrt(3);
            //outer x regions (pos. and neg.)
            else if ( abs ( m_yposMod[j] ) > m_Crystal.a / sqrt(3) && abs ( m_yposMod[j] ) < 2 * m_Crystal.a / sqrt(3) && m_xposMod[j] <= 0 )
            {
                m_xposMod[j] = m_xposMod[j] + m_Crystal.a / 2;
                if ( m_yposMod[j] < 0 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / ( 2 * sqrt(3) );
                else if ( m_yposMod[j] > 0 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / ( 2 * sqrt(3) );
            }
            else if ( abs ( m_yposMod[j] ) > m_Crystal.a / sqrt(3) && abs ( m_yposMod[j] ) < 2 * m_Crystal.a / sqrt(3) && m_xposMod[j] > 0 )
            {
                m_xposMod[j] = m_xposMod[j] - m_Crystal.a / 2;
                if ( m_yposMod[j] < 0 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / ( 2 * sqrt(3) );
                else if ( m_yposMod[j] > 0 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / ( 2 * sqrt(3) );
            }
            //remaining x regions
            else if ( abs ( m_yposMod[j] ) > m_Crystal.a / ( 2 * sqrt(3) ) && m_xposMod[j] /  abs ( m_yposMod[j] ) > sqrt(3) ) m_xposMod[j] = m_xposMod[j] - m_Crystal.a;
            else if ( abs ( m_yposMod[j] ) > m_Crystal.a / ( 2 * sqrt(3) ) && m_xposMod[j] /  abs ( m_yposMod[j] ) < -sqrt(3) ) m_xposMod[j] = m_xposMod[j] + m_Crystal.a;
            else if ( abs ( m_yposMod[j] ) > m_Crystal.a / ( 2 * sqrt(3) ) && abs ( m_yposMod[j] ) < m_Crystal.a / sqrt(3)
                      && m_xposMod[j] /  abs ( m_yposMod[j] ) < sqrt(3) && ( m_xposMod[j] - m_Crystal.a ) /  abs ( m_yposMod[j] ) > -sqrt(3) )
            {
                m_xposMod[j] = m_xposMod[j] - m_Crystal.a / 2;
                if ( m_yposMod[j] < 0 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / ( 2 * sqrt(3) );
                else if ( m_yposMod[j] > 0 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / ( 2 * sqrt(3) );
            }
            else if ( abs ( m_yposMod[j] ) > m_Crystal.a / ( 2 * sqrt(3) ) && abs ( m_yposMod[j] ) < m_Crystal.a / sqrt(3)
                      && m_xposMod[j] /  abs ( m_yposMod[j] ) > -sqrt(3) && ( m_xposMod[j] + m_Crystal.a ) /  abs ( m_yposMod[j] ) < sqrt(3) )
            {
                m_xposMod[j] = m_xposMod[j] + m_Crystal.a / 2;
                if ( m_yposMod[j] < 0 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / ( 2 * sqrt(3) );
                else if ( m_yposMod[j] > 0 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / ( 2 * sqrt(3) );
            }
            else if ( abs( m_yposMod[j] ) > 2 * m_Crystal.a / sqrt(3) && abs( m_yposMod[j] ) < 5 * m_Crystal.a / ( 2 * sqrt(3) ) && m_xposMod[j] / ( abs ( m_yposMod[j] ) - 2 * m_Crystal.a / sqrt(3) ) > sqrt(3)
                      && ( m_xposMod[j] - m_Crystal.a ) / ( abs ( m_yposMod[j] ) - 2 * m_Crystal.a / sqrt(3) ) < -sqrt(3) )
            {
                m_xposMod[j] = m_xposMod[j] - m_Crystal.a / 2;
                if ( m_yposMod[j] < 0 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / (2 * sqrt(3) );
                else if ( m_yposMod[j] > 0 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / ( 2 * sqrt(3) );
            }
            else if ( abs( m_yposMod[j] ) > 2 * m_Crystal.a / sqrt(3) && abs( m_yposMod[j] ) < 5 * m_Crystal.a / ( 2 * sqrt(3) ) && m_xposMod[j] / ( abs ( m_yposMod[j] ) - 2 * m_Crystal.a / sqrt(3) ) < -sqrt(3)
                      && ( m_xposMod[j] + m_Crystal.a ) / ( abs ( m_yposMod[j] ) - 2 * m_Crystal.a / sqrt(3) ) > sqrt(3) )
            {
                m_xposMod[j] = m_xposMod[j] + m_Crystal.a / 2;
                if ( m_yposMod[j] < 0 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / ( 2 * sqrt(3) );
                else if ( m_yposMod[j] > 0 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / ( 2 * sqrt(3) );
            }
            else if ( m_yposMod[j] < -2 * m_Crystal.a / sqrt(3) && abs ( m_xposMod[j] ) < m_Crystal.a / 2 ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / sqrt(3);
            else if ( m_yposMod[j] > 2 * m_Crystal.a / sqrt(3) && abs ( m_xposMod[j] ) < m_Crystal.a / 2 ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / sqrt(3);
            else if ( m_xposMod[j] < -m_Crystal.a / 2 )
            {
                m_xposMod[j] = m_xposMod[j] + m_Crystal.a;
                if ( ( m_xposMod[j] + m_Crystal.a ) / ( abs ( m_yposMod[j] ) - 2 * m_Crystal.a / sqrt(3) ) < sqrt(3) ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / sqrt(3);
                else if ( ( m_xposMod[j] + m_Crystal.a ) / ( abs ( m_yposMod[j] ) + 2 * m_Crystal.a / sqrt(3) ) > -sqrt(3) ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / sqrt(3);
            }
            else if ( m_xposMod[j] > m_Crystal.a / 2 )
            {
                m_xposMod[j] = m_xposMod[j] - m_Crystal.a;
                if ( ( m_xposMod[j] + m_Crystal.a ) / ( abs ( m_yposMod[j] ) - 2 * m_Crystal.a / sqrt(3) ) < sqrt(3) ) m_yposMod[j] = m_yposMod[j] - 3 * m_Crystal.a / sqrt(3);
                else if ( ( m_xposMod[j] + m_Crystal.a ) / ( abs ( m_yposMod[j] ) + 2 * m_Crystal.a / sqrt(3) ) > -sqrt(3) ) m_yposMod[j] = m_yposMod[j] + 3 * m_Crystal.a / sqrt(3);
            }
        }
    }

    for ( int i = 0; i < m_sPlotData.numPositions; i++ )
    {
        for ( int j = 0; j < m_sPrecision.numPositionsPrime; j++ )
        {
            if ( ( abs ( m_xposMod[i] - m_xposPrime[j] ) < m_Crystal.a * 1e-6 ) && ( abs ( m_yposMod[i] - m_yposPrime[j] ) < m_Crystal.a * 1e-6 ) ) m_posInt[i] = j;
        }
        assert ( m_posInt[i] >= 0 || lattice == 2 );
    }
}

void Response::CalculateResponse ( int PartNum, double kxlstart, double kylstart, double kxlend, double kylend,
                                   double kMin, double kMax, int ExchField, int DipolarField, int ElectricField, int FreqGF, int kTotal )
{
    //Declare variables
    //int const numPositionsPrime = sPrecision.numPositionsPrime;
    int mode, kCount = 0; /*MagicNum = 2 * spatial_sizePrime + 1;*/
    complex <double> mx, my, mxprime, myprime, Mxx0, Mxy0, Myx0, Myy0, Mxx, Mxy, Myx, Myy,
                     expkdotx, EigenOmega;
    VectorXcd mxvector = VectorXcd::Zero ( m_sPlotData.numPositions ),
              myvector = VectorXcd::Zero ( m_sPlotData.numPositions ),
              mxvectorPrime = VectorXcd::Zero ( m_sPrecision.numPositionsPrime ),
              myvectorPrime = VectorXcd::Zero ( m_sPrecision.numPositionsPrime ),
              GxxVector = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              GxyVector = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              GyxVector = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              GyyVector = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              GxxVectorTemp = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              GxyVectorTemp = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              GyxVectorTemp = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              GyyVectorTemp = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              Gxx = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              Gxy = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              Gyx = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              Gyy = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq ),
              Integral = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
    double kxl, kyl, Distribution;
    ofstream ERRORfile ( "ERROR_Info.dat", ios::out | ios::trunc );
    ofstream ProgressFile ( NameFILE ( "Info/Progress", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".dat" ).c_str(), ios::out | ios::trunc );

    //Check for existence of Progress files in Progress folder for continuing after stopping the program:
    string ProgressCheckStr = NameFILE ( "Progress/Progress", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".dat" );
    fstream ProgressCheckFile ( ProgressCheckStr.c_str(), ios::in | ios::out | ios::app );
    ProgressCheckFile.seekg ( 0, ios::beg );


    while ( !ProgressCheckFile.eof() )
    {
        ProgressCheckFile >> kxlstart;
        ProgressCheckFile.get();
        ProgressCheckFile >> kylstart;
        ProgressCheckFile.get();
    }
    ProgressCheckFile.close();
    ProgressCheckFile.open ( ProgressCheckStr.c_str(), ios::in | ios::out | ios::app );

    //Name green function files
    string m_GxxString = NameFILE ( "GreenFxn/Gxx", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".bin" );
    string m_GxyString = NameFILE ( "GreenFxn/Gxy", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".bin" );
    string m_GyxString = NameFILE ( "GreenFxn/Gyx", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".bin" );
    string m_GyyString = NameFILE ( "GreenFxn/Gyy", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".bin" );

    ifstream GxxFileRead ( m_GxxString.c_str(), ios::in );

    if ( GxxFileRead.is_open() )  //if above files exist, then read in old calculation data
    {
        GxxFileRead.close();

        ReadMatrix ( m_GxxString, m_GxxMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
        ReadMatrix ( m_GxyString, m_GxyMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
        ReadMatrix ( m_GyxString, m_GyxMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
        ReadMatrix ( m_GyyString, m_GyyMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
    }
    else GxxFileRead.close();


//The long part.
    //Name & open files for k-space contours:
    ofstream DispersionFile ( NameFILE ( "Contours/Dispersion", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".dat" ).c_str(), ios::out | ios::trunc );
    ofstream LinewidthFile ( NameFILE ( "Contours/Linewidth", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".dat" ).c_str(), ios::out | ios::trunc );

    if ( PartNum == 0 )
    {
        DispersionFile << "Col. 1: kxl" << endl <<
                          "Col. 2: kyl" << endl <<
                          "Col. 3: Re (OMEGA)" << endl << endl << endl;
        LinewidthFile  << "Col. 1: kxl" << endl <<
                          "Col. 2: kyl" << endl <<
                          "Col. 3: Im (OMEGA)" << endl << endl << endl;
    }

    //removed in favor of using actual frequency instead of dimensionless frequency
//    for ( int n = 0; n < m_sPlotData.NumFreq; n++ )
//    {
//        m_sPlotData.OMEGA0 ( n ) = m_Crystal.gyroRatio * m_Crystal.mu0 * m_Crystal.H0 * m_sPlotData.OMEGA0_nodim[n];
//        //cout << m_sPlotData.OMEGA0_nodim ( n ) << endl;
//    }
    time ( &m_start );
    for ( kxl = kxlstart; kxl < kxlend + 1e-6; kxl += m_sPrecision.stepsize )
    {
        //if ( abs ( kxl ) > 1e-6 ) continue;
        for ( kyl = kylstart; kyl < kylend + 1e-6; kyl += m_sPrecision.stepsize )
        {
            //The following if statement is only to be used for hexagonal lattice
            if ( lattice == 2 )
            {
                if ( abs ( kxl ) > 1 / sqrt(3) )
                {
                    if ( kxl * kyl > 0 && abs ( rotateycomp ( kxl, kyl, m_Crystal.Pi / 3 ) ) > 1. )
                    {
                        continue;
                    }
                    else if ( kxl * kyl < 0 && abs ( rotateycomp ( kxl, kyl, -m_Crystal.Pi / 3 ) ) > 1. )
                    {
                        continue;
                    }
                }
            }
            kCount++;

            if ( kCount <= kMin ) continue;    // if true, then this k-value was calcualated in a different part
            else if ( kCount > kMax ) break;   // don't go outside the BZ
            ProgressCheckFile << kxl << " " << kyl << endl;

            ///////////////Don't forget about this*****************************************:
            //kxl = 1e-50; kyl = 1e-50;
            double kzl = 0.0;
            Solve ( kxl, kyl, kzl, ExchField, DipolarField, ElectricField );

            //test vailidity of eigenvectors:
//            int modechoice = 0;
//            gsl_matrix_complex *MatProdComplex = gsl_matrix_complex_alloc ( 2 * numvectors, 2 * numvectors );
//            gsl_vector_complex *Result = gsl_vector_complex_alloc ( 2 * numvectors ), *evector = gsl_vector_complex_alloc ( 2 * numvectors );
//            for ( int i = 0; i < 2 * numvectors; i++ )
//            {
//                for ( int j = 0; j < 2 * numvectors; j++ )
//                {
//                    gsl_matrix_complex_set ( MatProdComplex, i, j, gsl_complex_rect ( gsl_matrix_get ( MatProd, i, j ), 0.0 ) );
//                }
//            }
//
//            double RealEvec, OtherEvec;
//            for ( modechoice = 0; modechoice < 2 * numvectors; modechoice++)
//            {
//                gsl_matrix_complex_get_col ( evector, evecs,  modechoice );
//                gsl_blas_zgemv ( CblasNoTrans, gsl_complex_inverse ( gsl_vector_complex_get ( evals, modechoice ) ), MatProdComplex, evector, gsl_complex_rect ( 0.0, 0.0 ), Result );
//                for ( int i = 0; i < 2 * numvectors; i++ )
//                {
//                    RealEvec = GSL_REAL ( gsl_matrix_complex_get ( evecs, i, modechoice ) );
//                    OtherEvec = GSL_REAL ( gsl_vector_complex_get ( evector, i ) );
//                    if ( abs ( RealEvec - OtherEvec ) > abs ( 0.001 * RealEvec ) )
//                    {
//                        cout << i << " " << RealEvec << " " << OtherEvec << endl;
//                    }
//                }
//            }
//
//            for ( int i = 0; i < numvectors; i++ )
//            {
//                ProgressFile << Gvx[i] * a / ( 2 * Pi ) << " " << Gvy[i] * a / ( 2 * Pi ) << " " << GSL_REAL ( gsl_matrix_complex_get ( evecs, i, 0 ) ) << " " << GSL_IMAG ( gsl_matrix_complex_get ( evecs, i, 0 ) ) << " "
//                             << GSL_REAL ( gsl_matrix_complex_get ( evecs, i + numvectors, 0 ) ) << " " << GSL_IMAG ( gsl_matrix_complex_get ( evecs, i + numvectors, 0 ) ) << endl;
//            }

            //return 50;
            //end of validity test

//            Output k-space contour data:
            DispersionFile << kxl << " " << kyl;
            LinewidthFile << kxl << " " << kyl;
            for ( mode = 0; mode < 2 * m_sPrecision.nummodes; mode++ )
            {
                DispersionFile << " " << m_evals ( m_evalsSortedIndex.at(mode) ).imag();
                LinewidthFile  << " " << m_evals ( m_evalsSortedIndex.at(mode) ).real();
            }
            DispersionFile << endl;
            LinewidthFile << endl;

            for ( mode = 0; mode < 2 * m_sPrecision.nummodes; mode++ )
            {
                //set phase factor for this mode
//                mx = spatial_mx ( 0.0, 0.0, mode );
//                phase = atan2 ( -GSL_IMAG ( mx ), GSL_REAL ( mx ) );
//                gsl_matrix_complex_scale ( evecs, gsl_complex_exp ( gsl_complex_rect ( 0.0, phase ) ) );

                //get frequency for this mode and rescale it
                EigenOmega = iii * m_evals ( m_evalsSortedIndex.at(mode) ) * -m_Crystal.gyroRatio * m_Crystal.mu0 * m_Crystal.H0; //( eval is dimensionless i*Omega )
                if ( EigenOmega.imag() > 0 ) ERRORfile << "Negative imaginary frequency at mode #" << kxl << " " << kyl << " " << mode << endl;
                //time ( &m_finish );
                //cout << m_finish - m_start << endl;
                for ( int n = 0; n < m_sPlotData.NumFreq; n++ )
                {
                    for ( int m = 0; m < m_sPlotData.temporalSteps; m++ )
                    {
                        Integral ( m + n * m_sPlotData.temporalSteps ) = FourierIntegral ( EigenOmega, ( m + 1 ) * m_sPlotData.temporal_stepsize, n, m_sPlotData.EndTime, m_sPlotData.OMEGA0, m_sPlotData, m_Crystal.Pi );
                    }
                }

                for ( int i = 0; i < m_sPrecision.numPositionsPrime; i++ )
                {
                    mxvectorPrime ( i ) = ComputeFields ( MAGNETIZATION_x, m_xposPrime[i], m_yposPrime[i], mode );
                    myvectorPrime ( i ) = ComputeFields ( MAGNETIZATION_y, m_xposPrime[i], m_yposPrime[i], mode );
                }

                if ( lattice == 2 )
                {
                    for ( int i = 0; i < m_sPlotData.numPositions; i++ )
                    {
                        if ( m_posInt[i] > 0 )
                        {
                            mxvector ( i ) = mxvectorPrime ( m_posInt[i] );
                            myvector ( i ) = myvectorPrime ( m_posInt[i] );
                        }
                        else if ( m_posInt[i] < 0 )
                        {
                            for ( int j = 0; j < i; j++ )
                            {
                                if ( ( abs ( m_xposMod[i] - m_xposMod[j] ) < m_Crystal.a * 1e-6 ) && ( abs ( m_yposMod[i] - m_yposMod[j] ) < m_Crystal.a * 1e-6 ) && m_posInt[i] < 0 )
                                {
                                    mxvector ( i ) = mxvector ( j );
                                    myvector ( i ) = myvector ( j );

                                    m_posInt[i] = 0;
                                }
                            }
                            if ( m_posInt[i] < 0 )
                            {
                                mxvector ( i ) = ComputeFields ( MAGNETIZATION_x, m_xposMod[i], m_yposMod[i], mode );
                                myvector ( i ) = ComputeFields ( MAGNETIZATION_y, m_xposMod[i], m_yposMod[i], mode );
                            }
                            m_posInt[i] = -1;
                        }
                    }
                }

                for ( int j = 0; j < ( m_sPlotData.numPositions ); j++ )           //x
                {
                    if ( lattice == 1 )
                    {
                        mx = mxvectorPrime ( m_posInt[j] );
                        my = myvectorPrime ( m_posInt[j] );
                    }
                    else if ( lattice == 2 )
                    {
                        mx = mxvector ( j );
                        my = myvector ( j );
                    }

                    if ( FreqGF == 1 )  //compute frequency dependent green's function ( set temporalsteps to zero )
                    {
                        expkdotx = exp ( iii * ( kx(kxl, m_Crystal) * m_xpos[j] + ky(kyl, m_Crystal) * m_ypos[j] ) );
                        mxprime = expkdotx * conj ( ComputeFields ( MAGNETIZATION_x, 0.0, 0.0, mode ) );
                        myprime = expkdotx * conj ( ComputeFields ( MAGNETIZATION_y, 0.0, 0.0, mode ) );

                        //linear basis:
                        Mxx0 = mx * mxprime;
                        Mxy0 = mx * myprime;
                        Myx0 = my * mxprime;
                        Myy0 = my * myprime;

                        for ( int n = 0; n < m_sPlotData.NumFreq; n++ )
                        {
                            Mxx = Mxx0 / ( m_sPlotData.OMEGA0 ( n ) - EigenOmega );
                            Mxy = Mxy0 / ( m_sPlotData.OMEGA0 ( n ) - EigenOmega );
                            Myx = Myx0 / ( m_sPlotData.OMEGA0 ( n ) - EigenOmega );
                            Myy = Myy0 / ( m_sPlotData.OMEGA0 ( n ) - EigenOmega );

                            m_GxxMatrix ( j, n ) += Mxx;
                            m_GxyMatrix ( j, n ) += Mxy;
                            m_GyxMatrix ( j, n ) += Myx;
                            m_GyyMatrix ( j, n ) += Myy;
                        }
                        continue;
                    }

                    for ( int i = 0; i < m_sPrecision.numPositionsPrime; i++ )       //xprime
                    {
                        //for Gaussian distribution
                        if ( m_xposPrime[i] * m_xposPrime[i] + m_yposPrime[i] * m_yposPrime[i] > m_Crystal.R * m_Crystal.R ) continue;
                        Distribution = DistributionFunction ( m_xposPrime[i], m_yposPrime[i], m_Crystal );

                        //for point source
//                        if ( ( abs ( m_xposPrime[i] ) < m_Crystal.a * 1e-6 ) && ( abs ( m_yposPrime[i] ) < m_Crystal.a * 1e-6 ) ) Distribution = 1.0;
//                        else continue;

                        expkdotx = exp ( iii * ( kx(kxl, m_Crystal) * ( m_xpos[j] - m_xposPrime[i] ) + ky(kyl, m_Crystal) * ( m_ypos[j] - m_yposPrime[i] ) ) );
                        mxprime = expkdotx * conj ( mxvectorPrime ( i ) );
                        myprime = expkdotx * conj ( myvectorPrime ( i ) );

                        //linear basis ( circular basis is done in Add_GF_VTK program ):
                        Mxx = mx * mxprime;
                        Mxy = mx * myprime;
                        Myx = my * mxprime;
                        Myy = my * myprime;

                        Gxx += Mxx * Distribution * Integral;
                        Gxy += Mxy * Distribution * Integral;
                        Gyx += Myx * Distribution * Integral;
                        Gyy += Myy * Distribution * Integral;
                    }
                    m_GxxMatrix.row ( j ) += Gxx;
                    m_GxyMatrix.row ( j ) += Gxy;
                    m_GyxMatrix.row ( j ) += Gyx;
                    m_GyyMatrix.row ( j ) += Gyy;
                    Gxx = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
                    Gxy = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
                    Gyx = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
                    Gyy = VectorXcd::Zero ( m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
                }
                ProgressFile << 50. * mode / m_sPrecision.nummodes << endl;
            }
            ProgressFile << endl << "Writing to file. Do not stop until ETA is shown." << endl;

            WriteMatrix ( m_GxxString, m_GxxMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
            WriteMatrix ( m_GxyString, m_GxyMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
            WriteMatrix ( m_GyxString, m_GyxMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
            WriteMatrix ( m_GyyString, m_GyyMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );

            time ( &m_finish );
            m_timeEst = difftime ( m_finish, m_start ) * ( kTotal * ( PartNum + 1 )  - kCount ) / ( kCount - kTotal * PartNum ) / 3600.0;
            if ( m_timeEst > 24 ) ProgressFile << kxl << " " << kyl << " ETA: " << m_timeEst / 24 << " days" << endl;
            else ProgressFile << endl << kxl << " " << kyl << " ETA: " << m_timeEst << " hours" << endl;
        }
    }
    ProgressFile.close(); ERRORfile.close(); DispersionFile.close(); LinewidthFile.close();
}

void Response::WriteToFile ( int PartNum )
{
    complex <double> scaleFactor;
    ofstream Info ( NameFILE ( "Info/Info", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".dat" ).c_str(), ios::out | ios::trunc );

    scaleFactor = pow ( m_sPrecision.spatial_stepsizePrime * m_sPrecision.stepsize * m_Crystal.Pi, 2 ) / ( 2 * m_Crystal.Pi * m_Crystal.sigmaSpace );
    m_GxxMatrix *= scaleFactor;
    m_GxyMatrix *= scaleFactor;
    m_GyxMatrix *= scaleFactor;
    m_GyyMatrix *= scaleFactor;

    time ( &m_finish );
    m_timeEst = difftime ( m_finish, m_start ) / 3600.0;

    //Output Other Information
    Info << "Nmax=" << m_sPrecision.Nmax << endl  <<
            "stepsize=" << m_sPrecision.stepsize << endl <<
            "f=" << m_Crystal.f << endl <<
            "a=" << m_Crystal.a << endl <<
            "mu0*H0=" << m_Crystal.mu0 * m_Crystal.H0 << endl <<
            "Msa=" << m_Mat_A.Ms << endl <<
            "Msb=" << m_Mat_B.Ms << endl <<
            "alphaa=" << m_Mat_A.alpha << endl <<
            "alphab=" << m_Mat_B.alpha << endl <<
            "Aa=" << m_Mat_A.A << endl <<
            "Ab=" << m_Mat_B.A << endl <<
            "spatial_size=" << m_sPlotData.spatial_size << endl <<
            "xmax=" << m_sPlotData.xmax << endl <<
            "ymax=" << m_sPlotData.ymax << endl <<
            "temporal_stepsize=" << m_sPlotData.temporal_stepsize << endl <<
            "temporalSteps=" << m_sPlotData.temporalSteps << endl << "time=" << m_timeEst << " hours";

    if ( PartNum == 0 )
    {
        ofstream Info0 ( "Info.dat", ios::out | ios::trunc );
        Info0 << "Nmax=" << m_sPrecision.Nmax << endl <<
                 "temporalSteps=" << m_sPlotData.temporalSteps << endl <<
                 "temporal_stepsize=" << m_sPlotData.temporal_stepsize << endl <<
                 "spatial_size=" << m_sPlotData.spatial_size << endl <<
                 "spatial_sizePrime=" << m_sPrecision.spatial_sizePrime << endl <<
                 "stepsize=" << m_sPrecision.stepsize << endl <<
                 "a=" << m_Crystal.a << endl <<
                 "xmax=" << m_sPlotData.xmax << endl <<
                 "ymax=" << m_sPlotData.ymax << endl <<
                 "NumFreq=" << m_sPlotData.NumFreq << endl;
        for ( int i = 0; i < m_sPlotData.NumFreq; i++ )
        {
//            if ( i == 0 ) Info0 << "omega= ( " << m_Crystal.gyroRatio * m_Crystal.mu0 * m_Crystal.H0 * m_sPlotData.OMEGA0_nodim [i];
//            else Info0 << ", " << m_Crystal.gyroRatio * m_Crystal.mu0 * m_Crystal.H0 * m_sPlotData.OMEGA0_nodim [i];

            //dimensionless omega
            if ( i == 0 ) Info0 << "OMEGA= ( " << m_sPlotData.OMEGA0 [i].real();
            else Info0 << ", " << m_sPlotData.OMEGA0 [i].real();
        }

        Info0 << " )" << endl <<
                 "f=" << m_Crystal.f << endl <<
                 "mu0*H0=" << m_Crystal.mu0 * m_Crystal.H0 << endl <<
                 "Msa=" << m_Mat_A.Ms << endl <<
                 "Msb=" << m_Mat_B.Ms << endl <<
                 "alphaa=" << m_Mat_A.alpha << endl <<
                 "alphab=" << m_Mat_B.alpha << endl <<
                 "Aa=" << m_Mat_A.A << endl <<
                 "Ab=" << m_Mat_B.A << endl <<
                 "temporal_stepsize=" << m_sPlotData.temporal_stepsize << endl <<
                 "endtime = " << m_sPlotData.EndTime << endl <<
                 "interface thickness = " << m_Crystal.t_Interface << endl <<
                 "interface layers = " << m_sPrecision.N_Interface << endl <<
                 "time=" << m_timeEst / 24 << " days";
        Info0.close();
    }

//Output for multiple parts:
    WriteMatrix ( m_GxxString.c_str(), m_GxxMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
    WriteMatrix ( m_GxyString.c_str(), m_GxyMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
    WriteMatrix ( m_GyxString.c_str(), m_GyxMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
    WriteMatrix ( m_GyyString.c_str(), m_GyyMatrix, m_sPlotData.numPositions, m_sPlotData.temporalSteps * m_sPlotData.NumFreq );
    Info.close();
}
