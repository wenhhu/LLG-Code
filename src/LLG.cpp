#include <cmath>
#include <iostream>
#include <cassert>

//gsl headers
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "constants.h"
#include "declarations.h"

using namespace std;

//constructor
//LandauLifshitzGilbert::LandauLifshitzGilbert ( Material Mat_A, Material Mat_B, Structure Crystal, int Nmax, int N_Interface ):
//    m_cPrecision ()
//{
//    SetParameters ( Mat_A, Mat_B, Crystal, Nmax, N_Interface );
//}

LandauLifshitzGilbert::LandauLifshitzGilbert ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision )
{
    SetParameters ( Mat_A, Mat_B, Crystal, sPrecision );
}

void LandauLifshitzGilbert::SetParameters ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision )
{
    m_kzl         = 0;
    m_Mat_A       = Mat_A;
    m_Mat_B       = Mat_B;
    m_Crystal     = Crystal;
    m_sPrecision  = sPrecision;
    m_sPrecision.numPositionsPrime     = ( 2 * m_sPrecision.spatial_sizePrime + 1 ) * ( 2 * m_sPrecision.spatial_sizePrime + 1 );
    m_sPrecision.spatial_stepsizePrime = 1 / m_sPrecision.spatial_sizePrime;

//    m_Hex_Gv_indexer <<  1, -1,
//                         0, -2,
//                        -1, -1,
//                        -1,  1,
//                         0,  2,
//                         1,  1;

    if ( m_Mat_A.Ms < 1e-8 )
    {
        m_Mat_A.Q = 0;
        m_Mat_A.ExchLength = 0;
    }
    else
    {
        m_Mat_A.Q = ( 2 * m_Mat_A.A ) / ( m_Mat_A.Ms * m_Crystal.mu0 * m_Crystal.H0 );
        m_Mat_A.ExchLength = 2 * m_Mat_A.A / ( m_Crystal.mu0 * m_Mat_A.Ms * m_Mat_A.Ms );
    }

    if ( m_Mat_B.Ms < 1e-8 )
    {
        m_Mat_B.Q = 0;
        m_Mat_B.ExchLength = 0;
    }
    else
    {
        m_Mat_B.Q = ( 2 * m_Mat_B.A ) / ( m_Mat_B.Ms * m_Crystal.mu0 * m_Crystal.H0 );
        m_Mat_B.ExchLength = 2 * m_Mat_B.A / ( m_Crystal.mu0 * m_Mat_B.Ms * m_Mat_B.Ms );
    }

    m_Mat_A.DMvector = 2 * m_Mat_A.A * m_Crystal.eCharge * m_Crystal.eField / ( m_Crystal.E_so * m_Mat_A.Ms * m_Mat_A.Ms );
    m_Mat_B.DMvector = 2 * m_Mat_B.A * m_Crystal.eCharge * m_Crystal.eField / ( m_Crystal.E_so * m_Mat_B.Ms * m_Mat_B.Ms );

    m_sPrecision.Nmax_Sum = m_sPrecision.Nmax + m_sPrecision.NmaxLimiter;


//    m_Mat_A.Q = ( 2 * m_Mat_A.A ) / ( m_Mat_A.Ms * m_Crystal.mu0 * m_Crystal.H0 );
//    m_Mat_B.Q = ( 2 * m_Mat_B.A ) / ( m_Mat_B.Ms * m_Crystal.mu0 * m_Crystal.H0 );
//    m_Mat_A.ExchLength = 2 * m_Mat_A.A / ( Crystal.mu0 * m_Mat_A.Ms * m_Mat_A.Ms );
//    m_Mat_B.ExchLength = 2 * m_Mat_B.A / ( Crystal.mu0 * m_Mat_B.Ms * m_Mat_B.Ms );

    if ( lattice == 1 )
    {
        m_numvectors = ( 2 * m_sPrecision.Nmax + 1 ) * ( 2 * m_sPrecision.Nmax + 1 );
        m_numvectors_Sum = ( 2 * m_sPrecision.Nmax_Sum + 1 ) * ( 2 * m_sPrecision.Nmax_Sum + 1 );
        m_Crystal.R = sqrt( m_Crystal.Ratio * ( m_Crystal.f * m_Crystal.a * m_Crystal.a ) / m_Crystal.Pi ) - 0.5 * m_Crystal.t_Interface;
    }
    else if ( lattice == 2 )
    {
        m_numvectors = ( 6 * ( m_sPrecision.Nmax * m_sPrecision.Nmax + m_sPrecision.Nmax ) ) / 2 + 1;
        m_numvectors_Sum = ( 6 * ( m_sPrecision.Nmax_Sum * m_sPrecision.Nmax_Sum + m_sPrecision.Nmax_Sum ) ) / 2 + 1;
        m_Crystal.R = m_Crystal.a * sqrt ( m_Crystal.Ratio * ( m_Crystal.f / m_Crystal.Pi ) * 3 / ( 2 * sqrt (3) ) ) - 0.5 * m_Crystal.t_Interface;
    }
    m_Crystal.sigmaSpace = m_Crystal.R / 4;

    //make sure cylinders are not touching:
    assert ( m_Crystal.R  + 0.5 * m_Crystal.t_Interface < m_Crystal.a / 2 );  //x-axis
    assert ( m_Crystal.R / m_Crystal.Ratio  + 0.5 * m_Crystal.t_Interface < m_Crystal.a / 2 );  //y-axis

    m_evalsSortedIndex.resize ( 2 * m_numvectors, 0 );
    m_Gvx.resize     ( m_numvectors, 0 );
    m_Gvy.resize     ( m_numvectors, 0 );
    m_Gvx_Sum.resize ( m_numvectors_Sum, 0 );
    m_Gvy_Sum.resize ( m_numvectors_Sum, 0 );
    m_Gv1_index.resize ( m_numvectors_Sum );
    m_Gv2_index.resize ( m_numvectors_Sum );
    m_Id.resize                  ( m_numvectors, vector <double> ( m_numvectors, 0 ) );
    m_Qarray.resize              ( m_numvectors, vector <double> ( m_numvectors, 0 ) );
    m_Msarray.resize             ( m_numvectors, vector <double> ( m_numvectors, 0 ) );
    m_ExchLengthArray.resize     ( m_numvectors, vector <double> ( m_numvectors, 0 ) );
    m_Msarray_Sum.resize         ( m_numvectors_Sum, vector <double> ( m_numvectors, 0 ) );
    m_ExchLengthArray_Sum.resize ( m_numvectors_Sum, vector <double> ( m_numvectors, 0 ) );
    m_DMarray_Sum.resize         ( m_numvectors_Sum, vector <double> ( m_numvectors, 0 ) );

    m_Mmatrix.resize                   ( 2 * m_numvectors, 2 * m_numvectors );
    m_alphaMatrix.resize               ( 2 * m_numvectors, 2 * m_numvectors );
    m_alphaMatrixInverse.resize        ( 2 * m_numvectors, 2 * m_numvectors );
    m_MatProd.resize                   ( 2 * m_numvectors, 2 * m_numvectors );
    m_MmatrixComplex.resize            ( 2 * m_numvectors, 2 * m_numvectors );
    m_alphaMatrixInverseComplex.resize ( 2 * m_numvectors, 2 * m_numvectors );

    m_evals.resize ( 2 * m_numvectors );
    m_evecs.resize ( 2 * m_numvectors, 2 * m_numvectors );

    CreateGvectors ();
    BuildMatrices ();
}

void LandauLifshitzGilbert::ExportData ( ofstream& filename )
{
    filename << "Nmax=" << m_sPrecision.Nmax << endl  <<
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
                "Ab=" << m_Mat_B.A << endl << endl << endl;
}

void LandauLifshitzGilbert::BubbleSort( int ElectricField )
{
    int temp, flag = 1;

    for ( int i = 0; i < 2 * m_numvectors; i++ )
    {
        m_evalsSortedIndex.at(i) = i;
    }

//For sorting by absolute value:
    while ( flag )
    {
        flag = 0;
        for ( int i = 0; i < 2 * m_numvectors - 1; i++ )
        {
            if ( abs ( m_evals ( m_evalsSortedIndex.at(i) ).imag() ) > abs ( m_evals ( m_evalsSortedIndex.at(i + 1) ).imag() ) )
            {
                temp = m_evalsSortedIndex.at(i);
                m_evalsSortedIndex.at(i) = m_evalsSortedIndex.at(i + 1);
                m_evalsSortedIndex.at(i + 1) = temp;
                flag = 1;
            }
        }
    }

//For sorting only positive frequencies:
    if ( ElectricField == 1 )
    {
        flag = 1;
        while ( flag )
        {
            flag = 0;
            for ( int i = 0; i < 2 * m_numvectors - 1; i++ )
            {
                if ( m_evals ( m_evalsSortedIndex.at(i) ).imag() < 0 && m_evals ( m_evalsSortedIndex.at(i + 1) ).imag() > 0 )
                {
                    temp = m_evalsSortedIndex.at(i);
                    m_evalsSortedIndex.at(i) = m_evalsSortedIndex.at(i + 1);
                    m_evalsSortedIndex.at(i + 1) = temp;
                    flag = 1;
                }
            }
        }
    }

//    for ( int i = 0; i < 2 * m_numvectors - 1; i+=2 )
//    {
//        if ( m_evals ( m_evalsSortedIndex.at(i) ).imag() < 0 &&
//             abs ( abs ( m_evals ( m_evalsSortedIndex.at(i) ).imag() ) -
//                   abs ( m_evals ( m_evalsSortedIndex.at(i + 1) ).imag() ) ) <=
//             1e-3 * abs ( m_evals ( m_evalsSortedIndex.at(i) ).imag() ) )
//        {
//            temp = m_evalsSortedIndex.at(i);
//            m_evalsSortedIndex.at(i) = m_evalsSortedIndex.at(i + 1);
//            m_evalsSortedIndex.at(i + 1) = temp;
//        }
//    }
}

//Reciprocal Lattice Vectors and k-vector
void LandauLifshitzGilbert::CreateGvectors (  )       //This is designed so that the summing of recip. lattice vectors is symmetric.
{
    if ( lattice == 1 )
    {
        for ( int l = 0; l < m_numvectors; l++ )
        {
            m_Gvx.at(l) = (( 2. * m_Crystal.Pi) / m_Crystal.a) * (fmod ( l, 2. * m_sPrecision.Nmax + 1. ) - m_sPrecision.Nmax);
            m_Gvy.at(l) = (( 2. * m_Crystal.Pi) / m_Crystal.a) * (floor ( l / (2. * m_sPrecision.Nmax + 1.) ) - m_sPrecision.Nmax);
        }

        for ( int l = 0; l < m_numvectors_Sum; l++ )
        {
            m_Gvx_Sum.at(l) = (( 2. * m_Crystal.Pi) / m_Crystal.a) * (fmod ( l, 2. * m_sPrecision.Nmax_Sum + 1. ) - m_sPrecision.Nmax_Sum);
            m_Gvy_Sum.at(l) = (( 2. * m_Crystal.Pi) / m_Crystal.a) * (floor ( l / (2. * m_sPrecision.Nmax_Sum + 1.) ) - m_sPrecision.Nmax_Sum);
        }
    }
    else if ( lattice == 2 )   //ensures # of recip. vectors chosen is hexagonally symmetric
    {
        double recipvector1y =   4 * m_Crystal.Pi / ( m_Crystal.a * sqrt(3) );            //Note: recipvector1x is zero
        double recipvector2x = ( 4 * m_Crystal.Pi / ( m_Crystal.a * sqrt(3) ) ) * 0.5 * sqrt(3);
        double recipvector2y = ( 4 * m_Crystal.Pi / ( m_Crystal.a * sqrt(3) ) ) * ( -0.5 );

//        double recipvectorx = ( 4 * m_Crystal.Pi / ( m_Crystal.a * sqrt(3) ) ) * 0.5 * sqrt(3),
//               recipvectory = ( 4 * m_Crystal.Pi / ( m_Crystal.a * sqrt(3) ) ) * 0.5;

        //using a different indexing scheme below for hex lattice to hoepfully work better with Lanczos Sigma factor
        m_Gvx[0] = 0;  m_Gvy[0] = 0;
        int vecCount1 = 1, vecCount2 = 0, count = 1, flag = 1;
        for ( int i = 1; i <= ( m_numvectors - 1 ) / 6; i++ )
        {
            for ( int j = 0; j < 6; j++ )
            {
                m_Gvx.at(count) = rotatexcomp ( vecCount2 * recipvector2x, vecCount1 * recipvector1y + vecCount2 * recipvector2y, j * m_Crystal.Pi / 3 );
                m_Gvy.at(count) = rotateycomp ( vecCount2 * recipvector2x, vecCount1 * recipvector1y + vecCount2 * recipvector2y, j * m_Crystal.Pi / 3 );
                count++;
            }
            vecCount2++;
            flag--;
            if ( flag == 0 )
            {
                vecCount1++;
                vecCount2 = 0;
                flag = vecCount1;
            }
        }

        m_Gvx_Sum[0] = 0;  m_Gvy_Sum[0] = 0;
        vecCount1 = 1; vecCount2 = 0; count = 1; flag = 1;
        for ( int i = 1; i <= ( m_numvectors_Sum - 1 ) / 6; i++ )
        {
            for ( int j = 0; j < 6; j++ )
            {
                m_Gvx_Sum.at(count) = rotatexcomp ( vecCount2 * recipvector2x, vecCount1 * recipvector1y + vecCount2 * recipvector2y, j * m_Crystal.Pi / 3 );
                m_Gvy_Sum.at(count) = rotateycomp ( vecCount2 * recipvector2x, vecCount1 * recipvector1y + vecCount2 * recipvector2y, j * m_Crystal.Pi / 3 );
                count++;
            }
            vecCount2++;
            flag--;
            if ( flag == 0 )
            {
                vecCount1++;
                vecCount2 = 0;
                flag = vecCount1;
            }
        }

        int m_Gv1_temp = 0, m_Gv2_temp = 1, index = 1;
        m_Gv1_index(0) = 0;
        m_Gv2_index(0) = 0;
        for ( int i = 0; i < m_sPrecision.Nmax_Sum; i++ )
        {
            m_Gv1_temp = 0;
            m_Gv2_temp = i + 1;
            for ( int m = 0; m <= i; m++ )
            {
                for ( int j = 0; j < 6; j++ )
                {
                    m_Gv1_index(index) = m_Gv1_temp;
                    m_Gv2_index(index) = m_Gv2_temp;

                    index++;
                }
                m_Gv1_temp++;
                m_Gv2_temp--;
            }
        }


//        int m_Gv1_temp = 0, m_Gv2_temp = 2, index = 1;
//        m_Gv1_index(0) = 0;
//        m_Gv2_index(0) = 0;
//        for ( int i = 0; i < m_sPrecision.Nmax_Sum; i++ )
//        {
//            for ( int j = 0; j < 6; j++ )
//            {
//                for ( int m = 0; m <= i; m++ )
//                {
//                    m_Gv1_index(index) = m_Gv1_temp;
//                    m_Gv2_index(index) = m_Gv2_temp;
//
//                    index++;
//                    m_Gv1_temp += m_Hex_Gv_indexer( j, 0 );
//                    m_Gv2_temp += m_Hex_Gv_indexer( j, 1 );
//                }
//            }
//            m_Gv2_temp += 2;
//        }

//        m_Gvx[0] = 0;  m_Gvy[0] = 0;
//        for ( int i = 1; i < m_numvectors; i++ )
//        {
//            m_Gvx.at(i) = rotatexcomp ( m_Gv1_index(i) * recipvectorx, m_Gv2_index(i) * recipvectory, 0 );
//            m_Gvy.at(i) = rotateycomp ( m_Gv1_index(i) * recipvectorx, m_Gv2_index(i) * recipvectory, 0 );
//        }
//
//        m_Gvx_Sum[0] = 0;  m_Gvy_Sum[0] = 0;
//        for ( int i = 1; i < m_numvectors_Sum; i++ )
//        {
//            m_Gvx_Sum.at(i) = rotatexcomp ( m_Gv1_index(i) * recipvectorx, m_Gv2_index(i) * recipvectory, 0 );
//            m_Gvy_Sum.at(i) = rotateycomp ( m_Gv1_index(i) * recipvectorx, m_Gv2_index(i) * recipvectory, 0 );
//        }
    }
}

//Fourier Transformed Variables
double LandauLifshitzGilbert::Structure_Factor ( int l, int j, int n )
{
    double RadiusShift = m_Crystal.t_Interface * ( m_sPrecision.N_Interface + 1 - n ) / m_sPrecision.N_Interface;
    switch ( m_Crystal.eInterface )
    {
        case INTERFACE_IDEAL:
            //ideal interface - not neccessary since R_N+1 = R?
            if ( j == -1 )
            {
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx.at(l) * m_Crystal.R, m_Gvy.at(l) * m_Crystal.R / m_Crystal.Ratio ) ) /
                                                           ( norm ( m_Gvx.at(l) * m_Crystal.R, m_Gvy.at(l) * m_Crystal.R / m_Crystal.Ratio ) );
            }
            else
            {
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( ( m_Gvx.at(l) - m_Gvx.at(j) ) * m_Crystal.R, ( m_Gvy.at(l) - m_Gvy.at(j) ) * m_Crystal.R / m_Crystal.Ratio ) ) /
                                                           ( norm ( ( m_Gvx.at(l) - m_Gvx.at(j) ) * m_Crystal.R, ( m_Gvy.at(l) - m_Gvy.at(j) ) * m_Crystal.R / m_Crystal.Ratio ) );
            }
        case INTERFACE_DIFFUSE:
            if ( j == -1 )
            {
//                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx.at(l), m_Gvy.at(l) ) * ( RadiusShift + m_Crystal.R ) ) /
//                                                           ( norm ( m_Gvx.at(l), m_Gvy.at(l) ) * ( RadiusShift + m_Crystal.R ) );
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx.at(l) * ( RadiusShift + m_Crystal.R ), m_Gvy.at(l) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) ) /
                                                           ( norm ( m_Gvx.at(l) * ( RadiusShift + m_Crystal.R ), m_Gvy.at(l) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) );
            }
            else
            {
//                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx.at(l) - m_Gvx.at(j), m_Gvy.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R ) ) /
//                                                           ( norm ( m_Gvx.at(l) - m_Gvx.at(j), m_Gvy.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R ) );
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( ( m_Gvx.at(l) - m_Gvx.at(j) ) * ( RadiusShift + m_Crystal.R ), ( m_Gvy.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) ) /
                                                           ( norm ( ( m_Gvx.at(l) - m_Gvx.at(j) ) * ( RadiusShift + m_Crystal.R ), ( m_Gvy.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) );
            }
        default:
            cerr << "Invalid case in Structure_Factor.";
            exit ( EXIT_FAILURE );
    }
}

double LandauLifshitzGilbert::C_n ( double C_a, double C_b, int n )
{
    //linear interface
    return C_b * ( m_sPrecision.N_Interface + 1 - n ) / ( m_sPrecision.N_Interface + 1 ) + C_a * n / ( m_sPrecision.N_Interface + 1 );

    //sinusoidal interface
    //return ( C_a + C_b ) / 2 - ( C_b - C_a ) / 2 * cos ( Crystal.Pi * ( m_sPrecision.N_Interface + 1 - n ) / ( m_sPrecision.N_Interface + 1 ) ) );
}

double LandauLifshitzGilbert::LanczosSigmaFactor ( int l, int j )
{
    double x=0, y=0;
    int m=0, n=0;

    if ( j == -1 )
    {
        if ( lattice == 1 )
        {
            m = (fmod ( l, 2. * m_sPrecision.Nmax + 1. ) - m_sPrecision.Nmax);
            n = (floor ( l / (2. * m_sPrecision.Nmax + 1.) ) - m_sPrecision.Nmax);
        }
        else if ( lattice == 2 )
        {
            m = m_Gv1_index(l);
            n = m_Gv2_index(l);
        }
    }
    else
    {
        if ( lattice == 1 )
        {
            m = (fmod ( l, 2. * m_sPrecision.Nmax + 1. ) - m_sPrecision.Nmax) - (fmod ( j, 2. * m_sPrecision.Nmax + 1. ) - m_sPrecision.Nmax);
            n = (floor ( l / (2. * m_sPrecision.Nmax + 1.) ) - m_sPrecision.Nmax) - (floor ( j / (2. * m_sPrecision.Nmax + 1.) ) - m_sPrecision.Nmax);
        }
        else if ( lattice == 2 )
        {
            m = m_Gv1_index(l) - m_Gv1_index(j);
            n = m_Gv2_index(l) - m_Gv2_index(j);
        }
    }
    //change Nmax to 2*Nmax to account for G_i-G_j in LLG matrix equation 2014-11-03
    if ( lattice == 1 )
    {
        x = m_Crystal.Pi * ( 1. * m / ( m_sPrecision.NmaxLimiter + 1. ) );
        y = m_Crystal.Pi * ( 1. * n / ( m_sPrecision.NmaxLimiter + 1. ) );
    }
    else if ( lattice == 2 )
    {
        x = m_Crystal.Pi * ( 1. * m / ( m_sPrecision.NmaxLimiter + 1. ) );
        y = m_Crystal.Pi * ( 1. * n / ( m_sPrecision.NmaxLimiter + 1. ) );
    }

    if      ( abs(x) > m_Crystal.Pi || abs(y) > m_Crystal.Pi ) return 0;
    else if ( m == 0 && n == 0 ) return 1;
    else if ( m == 0 )           return sin ( y ) / y;
    else if ( n == 0 )           return sin ( x ) / x;
    else                         return sin ( x ) * sin ( y ) / ( x * y );
}

double LandauLifshitzGilbert::FT_Variable ( double C_a, double C_b, int l, int j )
{
    double sum = 0;

    switch ( m_Crystal.eInterface )
    {
        case INTERFACE_IDEAL:
            //ideal interface
            if ( ( l == j ) || ( abs (m_Gvx.at(l)) < 0.001 / m_Crystal.a && abs (m_Gvy.at(l)) < 0.001 / m_Crystal.a && j == -1 ) )
                return LanczosSigmaFactor ( l, j ) * ( C_a * m_Crystal.f + C_b * (1. - m_Crystal.f ) );
            else
            {
                return LanczosSigmaFactor ( l, j ) * (C_a - C_b) * Structure_Factor ( l, j, m_sPrecision.N_Interface + 1 );
            }

        case INTERFACE_DIFFUSE:
            //linear interface
            if ( ( l == j ) || ( abs (m_Gvx.at(l)) < 0.001 / m_Crystal.a && abs (m_Gvy.at(l)) < 0.001 / m_Crystal.a && j == -1 ) )
            {
                for ( int n = 1; n < m_sPrecision.N_Interface + 1; n++ )
                {
                    sum += C_n ( C_a, C_b, n ) * ( m_Crystal.Pi / ( m_Crystal.a * m_Crystal.a ) ) * ( pow ( m_Crystal.t_Interface * ( m_sPrecision.N_Interface + 1 - n ) / m_sPrecision.N_Interface + m_Crystal.R, 2 ) -
                                                                                                      pow ( m_Crystal.t_Interface * ( m_sPrecision.N_Interface -     n ) / m_sPrecision.N_Interface + m_Crystal.R, 2 ) );
                }
                sum += C_a * m_Crystal.f + C_b * ( 1 - m_Crystal.Pi / ( m_Crystal.a * m_Crystal.a ) * pow ( m_Crystal.R + m_Crystal.t_Interface, 2 ) );
                return sum;
            }
            else
            {
                for ( int n = 1; n < m_sPrecision.N_Interface + 2; n++ )
                {
                    sum += ( C_n ( C_a, C_b, n ) - C_n ( C_a, C_b, n - 1 ) ) * Structure_Factor ( l, j, n );
                }
                return sum;
            }
        default:
            cerr << "Invalid case in FT_Variable.";
            exit ( EXIT_FAILURE );
    }

    //sinusoidal interface (to be implemented)
}

//Fourier transformed terms for use in exchange term summation
double LandauLifshitzGilbert::Structure_Factor_Sum ( int l, int j, int n )
{
    double RadiusShift = m_Crystal.t_Interface * ( m_sPrecision.N_Interface + 1 - n ) / m_sPrecision.N_Interface;
    switch ( m_Crystal.eInterface )
    {
        case INTERFACE_IDEAL:
            //ideal interface - not neccessary since R_N+1 = R?
            if ( j == -1 )
            {
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx.at(l) * m_Crystal.R, m_Gvy.at(l) * m_Crystal.R / m_Crystal.Ratio ) ) /
                                                           ( norm ( m_Gvx.at(l) * m_Crystal.R, m_Gvy.at(l) * m_Crystal.R / m_Crystal.Ratio ) );
            }
            else
            {
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( ( m_Gvx_Sum.at(l) - m_Gvx.at(j) ) * m_Crystal.R, ( m_Gvy_Sum.at(l) - m_Gvy.at(j) ) * m_Crystal.R / m_Crystal.Ratio ) ) /
                                                           ( norm ( ( m_Gvx_Sum.at(l) - m_Gvx.at(j) ) * m_Crystal.R, ( m_Gvy_Sum.at(l) - m_Gvy.at(j) ) * m_Crystal.R / m_Crystal.Ratio ) );
            }
        case INTERFACE_DIFFUSE:
            if ( j == -1 )
            {
//                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx.at(l), m_Gvy.at(l) ) * ( RadiusShift + m_Crystal.R ) ) /
//                                                           ( norm ( m_Gvx.at(l), m_Gvy.at(l) ) * ( RadiusShift + m_Crystal.R ) );
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx_Sum.at(l) * ( RadiusShift + m_Crystal.R ), m_Gvy_Sum.at(l) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) ) /
                                                           ( norm ( m_Gvx_Sum.at(l) * ( RadiusShift + m_Crystal.R ), m_Gvy_Sum.at(l) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) );
            }
            else
            {
//                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( m_Gvx.at(l) - m_Gvx.at(j), m_Gvy.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R ) ) /
//                                                           ( norm ( m_Gvx.at(l) - m_Gvx.at(j), m_Gvy.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R ) );
                return 2. * m_Crystal.f * gsl_sf_bessel_J1 ( norm ( ( m_Gvx_Sum.at(l) - m_Gvx.at(j) ) * ( RadiusShift + m_Crystal.R ), ( m_Gvy_Sum.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) ) /
                                                           ( norm ( ( m_Gvx_Sum.at(l) - m_Gvx.at(j) ) * ( RadiusShift + m_Crystal.R ), ( m_Gvy_Sum.at(l) - m_Gvy.at(j) ) * ( RadiusShift + m_Crystal.R / m_Crystal.Ratio ) ) );
            }
        default:
            cerr << "Invalid case in Structure_Factor.";
            exit ( EXIT_FAILURE );
    }
}

double LandauLifshitzGilbert::C_n_Sum ( double C_a, double C_b, int n )
{
    //linear interface
    return C_b * ( m_sPrecision.N_Interface + 1 - n ) / ( m_sPrecision.N_Interface + 1 ) + C_a * n / ( m_sPrecision.N_Interface + 1 );

    //sinusoidal interface
    //return ( C_a + C_b ) / 2 - ( C_b - C_a ) / 2 * cos ( Crystal.Pi * ( m_sPrecision.N_Interface + 1 - n ) / ( m_sPrecision.N_Interface + 1 ) ) );
}

double LandauLifshitzGilbert::LanczosSigmaFactor_Sum ( int l, int j )
{
    double x=0, y=0;
    int m=0, n=0;

    if ( j == -1 )
    {
        if ( lattice == 1 )
        {
            m = (fmod ( l, 2. * m_sPrecision.Nmax + 1. ) - m_sPrecision.Nmax);
            n = (floor ( l / (2. * m_sPrecision.Nmax + 1.) ) - m_sPrecision.Nmax);
        }
        else if ( lattice == 2 )
        {
            m = m_Gv1_index(l);
            n = m_Gv2_index(l);
        }
    }
    else
    {
        if ( lattice == 1 )
        {
            m = (fmod ( l, 2. * m_sPrecision.Nmax_Sum + 1. ) - m_sPrecision.Nmax_Sum) - (fmod ( j, 2. * m_sPrecision.Nmax + 1. ) - m_sPrecision.Nmax);
            n = (floor ( l / (2. * m_sPrecision.Nmax_Sum + 1.) ) - m_sPrecision.Nmax_Sum) - (floor ( j / (2. * m_sPrecision.Nmax + 1.) ) - m_sPrecision.Nmax);
        }
        else if ( lattice == 2 )
        {
            m = m_Gv1_index(l) - m_Gv1_index(j);
            n = m_Gv2_index(l) - m_Gv2_index(j);
        }
    }
    //change Nmax to 2*Nmax to account for G_i-G_j in LLG matrix equation 2014-11-03
    if ( lattice == 1 )
    {
        x = m_Crystal.Pi * ( 1. * m / ( m_sPrecision.NmaxLimiter + 1. ) );
        y = m_Crystal.Pi * ( 1. * n / ( m_sPrecision.NmaxLimiter + 1. ) );
    }
    else if ( lattice == 2 )
    {
        x = m_Crystal.Pi * ( 1. * m / ( m_sPrecision.NmaxLimiter + 1. ) );
        y = m_Crystal.Pi * ( 1. * n / ( m_sPrecision.NmaxLimiter + 1. ) );
    }

    if      ( abs (x) > m_Crystal.Pi || abs(y) > m_Crystal.Pi ) return 0;
    else if ( m == 0 && n == 0 ) return 1;
    else if ( m == 0 )           return sin ( y ) / y;
    else if ( n == 0 )           return sin ( x ) / x;
    else                         return sin ( x ) * sin ( y ) / ( x * y );
}

double LandauLifshitzGilbert::FT_Variable_Sum ( double C_a, double C_b, int l, int j )
{
    double sum = 0;

    switch ( m_Crystal.eInterface )
    {
        case INTERFACE_IDEAL:
            //ideal interface
            if ( ( abs ( m_Gvx_Sum.at(l) - m_Gvx.at(j) ) < 0.001 / m_Crystal.a &&
                   abs ( m_Gvy_Sum.at(l) - m_Gvy.at(j) ) < 0.001 / m_Crystal.a ) ||
                 ( abs ( m_Gvx_Sum.at(l) ) < 0.001 / m_Crystal.a && abs (m_Gvy_Sum.at(l)) < 0.001 / m_Crystal.a && j == -1 ) )
                return LanczosSigmaFactor_Sum ( l, j ) * ( C_a * m_Crystal.f + C_b * (1. - m_Crystal.f ) );
            else
            {
                return LanczosSigmaFactor_Sum ( l, j ) * (C_a - C_b) * Structure_Factor_Sum ( l, j, m_sPrecision.N_Interface + 1 );
            }

        case INTERFACE_DIFFUSE:
            //linear interface
            if ( ( abs ( m_Gvx_Sum.at(l) - m_Gvx.at(j) ) < 0.001 / m_Crystal.a &&
                   abs ( m_Gvy_Sum.at(l) - m_Gvy.at(j) ) < 0.001 / m_Crystal.a  ) ||
                 ( abs ( m_Gvx_Sum.at(l) ) < 0.001 / m_Crystal.a && abs (m_Gvy_Sum.at(l)) < 0.001 / m_Crystal.a && j == -1 ) )
            {
                for ( int n = 1; n < m_sPrecision.N_Interface + 1; n++ )
                {
                    sum += C_n_Sum ( C_a, C_b, n ) * ( m_Crystal.Pi / ( m_Crystal.a * m_Crystal.a ) ) * ( pow ( m_Crystal.t_Interface * ( m_sPrecision.N_Interface + 1 - n ) / m_sPrecision.N_Interface + m_Crystal.R, 2 ) -
                                                                                                          pow ( m_Crystal.t_Interface * ( m_sPrecision.N_Interface -     n ) / m_sPrecision.N_Interface + m_Crystal.R, 2 ) );
                }
                sum += C_a * m_Crystal.f + C_b * ( 1 - m_Crystal.Pi / ( m_Crystal.a * m_Crystal.a ) * pow ( m_Crystal.R + m_Crystal.t_Interface, 2 ) );
                return sum;
            }
            else
            {
                for ( int n = 1; n < m_sPrecision.N_Interface + 2; n++ )
                {
                    sum += ( C_n_Sum ( C_a, C_b, n ) - C_n_Sum ( C_a, C_b, n - 1 ) ) * Structure_Factor_Sum ( l, j, n );
                }
                return sum;
            }
        default:
            cerr << "Invalid case in FT_Variable.";
            exit ( EXIT_FAILURE );
    }

    //sinusoidal interface (to be implemented)
}

void LandauLifshitzGilbert::BuildMatrices ()
{
    //gsl_matrix *alphaMatrix = gsl_matrix_alloc ( 2 * m_numvectors, 2 * m_numvectors );

    //CreateGvectors ( Nmax, numvectors );

//Create Matrices independent of kx,ky:
    for ( int i = 0; i < m_numvectors ; i++)
    {
        for ( int j = 0; j < m_numvectors ; j++)
        {
            if ( i == j )
                m_Id[i][j] = 1;
            else
                m_Id[i][j] = 0;
            m_Qarray[i][j] = FT_Variable ( m_Mat_A.Q, m_Mat_B.Q, i, j );
            m_Msarray[i][j] = FT_Variable ( m_Mat_A.Ms, m_Mat_B.Ms, i, j ) / m_Crystal.H0;
            m_ExchLengthArray[i][j] = FT_Variable ( m_Mat_A.ExchLength, m_Mat_B.ExchLength, i, j );
        }
    }
    for ( int i = 0; i < m_numvectors_Sum; i++)
    {
        for ( int j = 0; j < m_numvectors; j++)
        {
            m_Msarray_Sum[i][j] = FT_Variable_Sum ( m_Mat_A.Ms, m_Mat_B.Ms, i, j ) / m_Crystal.H0;
            m_ExchLengthArray_Sum[i][j] = FT_Variable_Sum ( m_Mat_A.ExchLength, m_Mat_B.ExchLength, i, j );
            m_DMarray_Sum[i][j] = FT_Variable_Sum ( m_Mat_A.DMvector, m_Mat_B.DMvector, i, j );
        }
    }

//Create Alpha Matrix
    //set upper left and lower right blocks to identity matrix
    m_alphaMatrix.topLeftCorner     ( m_numvectors, m_numvectors ).setIdentity();
    m_alphaMatrix.bottomRightCorner ( m_numvectors, m_numvectors ).setIdentity();
    for ( int i = 0; i < m_numvectors ; i++)
    {
        for ( int j = 0; j < m_numvectors ; j++)
        {
            m_alphaMatrix.topRightCorner ( m_numvectors, m_numvectors ) ( i, j ) = FT_Variable ( m_Mat_A.alpha, m_Mat_B.alpha, i, j );
        }
    }
    //set lower left block to negative of upper right block
    m_alphaMatrix.bottomLeftCorner ( m_numvectors, m_numvectors ) = -m_alphaMatrix.topRightCorner ( m_numvectors, m_numvectors );

//Invert Alpha matrix
    if ( ( m_Mat_A.alpha < 1e-12 ) && ( m_Mat_B.alpha < 1e-12 ) )
    {
        m_alphaMatrixInverse.setIdentity();
    }
    else
    {
        m_alphaMatrixInverse = m_alphaMatrix.inverse();
    }
}

//Exchange term compuatation for matrix, following Krawczyk and Puszkarski: ( obsolete - ComputeSumTerms is used in place of this function )
double LandauLifshitzGilbert::ComputeExchange ( int i, int j )
{
    double sum = 0;

    for ( int l = 0; l < m_numvectors_Sum; l++ )
    {
        sum += ( ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) * ( kx(m_kxl, m_Crystal) + m_Gvx_Sum.at(l) ) +
                   ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) * ( ky(m_kyl, m_Crystal) + m_Gvy_Sum.at(l) ) ) -
                 ( ( m_Gvx.at(i) - m_Gvx_Sum.at(l) ) * ( m_Gvx.at(i) - m_Gvx.at(j) ) +
                   ( m_Gvy.at(i) - m_Gvy_Sum.at(l) ) * ( m_Gvy.at(i) - m_Gvy.at(j) ) ) ) *
               m_ExchLengthArray_Sum [l][j] * m_Msarray_Sum [l][i];
    }
    return sum;
}

//term arising from the DM interaction with electric field  ( obsolete - ComputeSumTerms is used in place of this function )
complex <double> LandauLifshitzGilbert::ComputeDMterm ( int i, int j, int ExchField )
{
    complex <double> sum1 = 0;

    if ( ExchField == 0 || ExchField == 2 )
    {
        for ( int l = 0; l < m_numvectors_Sum; l++ )
        {
            sum1 += m_Msarray_Sum [l][i] * m_DMarray_Sum [l][j] *
                   ( m_Gvx.at(j) + kx ( m_kxl, m_Crystal ) ); //E in y direction, e_ij in x direction
        }
    }

    else if ( ExchField == 1 )
    {
        for ( int l = 0; l < m_numvectors_Sum; l++ )
        {
            sum1 += m_Msarray_Sum [l][i] * m_DMarray_Sum [l][j] *
                   ( m_Gvx.at(j) + kx ( m_kxl, m_Crystal ) ); //E in y direction, e_ij in x direction

//            sum2 += ( ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) * ( kx(m_kxl, m_Crystal) + m_Gvx_Sum.at(l) ) +
//                   ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) * ( ky(m_kyl, m_Crystal) + m_Gvy_Sum.at(l) ) ) -
//                 ( ( m_Gvx.at(i) - m_Gvx_Sum.at(l) ) * ( m_Gvx.at(i) - m_Gvx.at(j) ) +
//                   ( m_Gvy.at(i) - m_Gvy_Sum.at(l) ) * ( m_Gvy.at(i) - m_Gvy.at(j) ) ) ) *
//               m_ExchLengthArray_Sum [l][j] * m_Msarray_Sum [l][i];
        }
    }
    return 2.0 * iii * sum1 / m_Crystal.mu0;
}

void LandauLifshitzGilbert::ComputeSumTerms ( int i, int j, int ExchField, int ElectricField )  //computes both exchange and DM terms at once to save computation time
{
    m_Crystal.DMTerm = 0;
    m_Crystal.ExchangeTerm = 0;

    if ( ExchField == 1 && ElectricField == 0 )  // H_ex !=0 and E = 0
    {
        for ( int l = 0; l < m_numvectors_Sum; l++ )
        {
            m_Crystal.ExchangeTerm += ( ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) * ( kx(m_kxl, m_Crystal) + m_Gvx_Sum.at(l) ) +
                                          ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) * ( ky(m_kyl, m_Crystal) + m_Gvy_Sum.at(l) ) +
                                          kz(m_kzl, m_Crystal) * kz(m_kzl, m_Crystal)                                         ) -
                                        ( ( m_Gvx.at(i) - m_Gvx_Sum.at(l) ) * ( m_Gvx.at(i) - m_Gvx.at(j) ) +
                                          ( m_Gvy.at(i) - m_Gvy_Sum.at(l) ) * ( m_Gvy.at(i) - m_Gvy.at(j) ) ) ) *
                                      m_ExchLengthArray_Sum [l][j] * m_Msarray_Sum [l][i];
        }
    }

    else if ( ( ExchField == 0 || ExchField == 2 ) && ElectricField == 1 )  // H_ex = 0 and E != 0
    {
        for ( int l = 0; l < m_numvectors_Sum; l++ )
        {
            m_Crystal.DMTerm += m_Msarray_Sum [l][i] * m_DMarray_Sum [l][j] *
                   ( m_Gvx.at(j) + kx ( m_kxl, m_Crystal ) ); //E in y direction, e_ij in x direction
        }
    }

    else if ( ExchField == 1 && ElectricField == 1 ) // H_ex != 0 and E != 0
    {
        for ( int l = 0; l < m_numvectors_Sum; l++ )
        {
            m_Crystal.ExchangeTerm += ( ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) * ( kx(m_kxl, m_Crystal) + m_Gvx_Sum.at(l) ) +
                                          ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) * ( ky(m_kyl, m_Crystal) + m_Gvy_Sum.at(l) ) +
                                          kz(m_kzl, m_Crystal) * kz(m_kzl, m_Crystal)                                         ) -
                                        ( ( m_Gvx.at(i) - m_Gvx_Sum.at(l) ) * ( m_Gvx.at(i) - m_Gvx.at(j) ) +
                                          ( m_Gvy.at(i) - m_Gvy_Sum.at(l) ) * ( m_Gvy.at(i) - m_Gvy.at(j) ) ) ) *
                                      m_ExchLengthArray_Sum [l][j] * m_Msarray_Sum [l][i];

            m_Crystal.DMTerm += m_Msarray_Sum [l][i] * m_DMarray_Sum [l][j] *
                   ( m_Gvx.at(j) + kx ( m_kxl, m_Crystal ) ); //E in y direction, e_ij in x direction
        }
    }

    m_Crystal.DMTerm = 2.0 * iii * m_Crystal.DMTerm / m_Crystal.mu0;
}

double LandauLifshitzGilbert::ComputeSlabDemagField ( int i, int j )
{
    switch ( m_Crystal.eSlabType )
    {
        case ( FINITE ):
            return -4 * m_Crystal.Pi * m_Msarray[i][j] / ( cosh ( m_Crystal.SlabWidth * norm ( ( m_Gvx.at(i) - m_Gvx.at(j) ), ( m_Gvy.at(i) - m_Gvy.at(j) ) ) ) + sinh ( m_Crystal.SlabWidth * norm ( ( m_Gvx.at(i) - m_Gvx.at(j) ), ( m_Gvy.at(i) - m_Gvy.at(j) ) ) ) );

        case ( INFINITE ):
            return 0;
        default:
            cerr << "Invalid case in function ComputeSlabDemagField";
            exit ( EXIT_FAILURE );
    }
}

//Alternate matrix, computed using exchange length instead of Q:
//Now contains both forms.  Set using ExchField=1 for derived and 2 for alternate.

void LandauLifshitzGilbert::CreateMmatrix ( int ExchField, int DipolarField, int ElectricField )
{
    //complex <double> DMterm;

    if ( ExchField == 1 || ExchField == 0 )
    {
        //double ExchangeTerm;
        for ( int i = 0; i < m_numvectors; i++)
        {
            for ( int j = 0; j < m_numvectors ; j++)
            {
                ComputeSumTerms( i ,j, ExchField, ElectricField );

                if ( abs (m_kxl) < 1e-5 && abs (m_kyl) < 1e-5 && abs (m_kzl) < 1e-5 && abs (m_Gvx.at(j)) < 0.001 / m_Crystal.a && abs (m_Gvy.at(j)) < 0.001 / m_Crystal.a )
                {
                    //if ( ExchField == 0 ) ExchangeTerm = 0;
                    //else ExchangeTerm = ComputeExchange ( i, j );

                    m_Mmatrix.topLeftCorner    ( m_numvectors, m_numvectors ) ( i, j ) = 0;

                    m_Mmatrix.topRightCorner   ( m_numvectors, m_numvectors ) ( i, j ) = m_Id[i][j] -
                                                                                         ComputeSlabDemagField ( i, j ) +
                                                                                         m_Crystal.ExchangeTerm;

                    m_Mmatrix.bottomLeftCorner ( m_numvectors, m_numvectors ) ( i, j ) = -m_Mmatrix.topRightCorner ( m_numvectors, m_numvectors ) ( i, j );

                    //if      ( ElectricField == 0 ) DMterm = 0;
                    //else if ( ElectricField == 1 ) DMterm = ComputeDMterm ( i, j, ExchField );

                    m_MmatrixComplex.topLeftCorner     ( m_numvectors, m_numvectors ) ( i, j ) = -m_Crystal.DMTerm;
                    m_MmatrixComplex.topRightCorner    ( m_numvectors, m_numvectors ) ( i, j ) = 0;
                    m_MmatrixComplex.bottomLeftCorner  ( m_numvectors, m_numvectors ) ( i, j ) = 0;
                    m_MmatrixComplex.bottomRightCorner ( m_numvectors, m_numvectors ) ( i, j ) = -m_Crystal.DMTerm;
                }
                else
                {
                    //if ( ExchField == 0 ) ExchangeTerm = 0;
                    //else ExchangeTerm = ComputeExchange ( i, j );

                    m_Mmatrix.topLeftCorner    ( m_numvectors, m_numvectors ) ( i, j ) = ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) *
                                                                                           ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) /
                                                                                           ( pow ( norm ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), ky(m_kyl, m_Crystal) + m_Gvy.at(j) ), 2 ) +
                                                                                             pow ( kz ( m_kzl, m_Crystal ), 2 ) ) ) *
                                                                                         m_Msarray[i][j] * DipolarField;

                    m_Mmatrix.topRightCorner   ( m_numvectors, m_numvectors ) ( i, j ) = m_Id[i][j] -
                                                                                         ComputeSlabDemagField ( i, j ) +
                                                                                         pow ( ky(m_kyl, m_Crystal) + m_Gvy.at(j), 2 ) /
                                                                                         ( pow ( norm ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), ky(m_kyl, m_Crystal) + m_Gvy.at(j) ), 2 ) +
                                                                                         pow ( kz ( m_kzl, m_Crystal ), 2 ) ) *
                                                                                         m_Msarray[i][j] * DipolarField +
                                                                                         m_Crystal.ExchangeTerm;

                    m_Mmatrix.bottomLeftCorner ( m_numvectors, m_numvectors ) ( i, j ) = -m_Id[i][j] +
                                                                                         ComputeSlabDemagField ( i, j ) -
                                                                                         pow ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), 2 ) /
                                                                                         ( pow ( norm ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), ky(m_kyl, m_Crystal) + m_Gvy.at(j) ), 2 ) +
                                                                                         pow ( kz ( m_kzl, m_Crystal ), 2 ) ) *
                                                                                         m_Msarray[i][j] * DipolarField -
                                                                                         m_Crystal.ExchangeTerm;

                    //if      ( ElectricField == 0 ) DMterm = 0;
                    //else if ( ElectricField == 1 ) DMterm = ComputeDMterm ( i, j, ExchField );

                    m_MmatrixComplex.topLeftCorner     ( m_numvectors, m_numvectors ) ( i, j ) = -m_Crystal.DMTerm;
                    m_MmatrixComplex.topRightCorner    ( m_numvectors, m_numvectors ) ( i, j ) = 0;
                    m_MmatrixComplex.bottomLeftCorner  ( m_numvectors, m_numvectors ) ( i, j ) = 0;
                    m_MmatrixComplex.bottomRightCorner ( m_numvectors, m_numvectors ) ( i, j ) = -m_Crystal.DMTerm;
                }
            }
        }
        m_Mmatrix.bottomRightCorner ( m_numvectors, m_numvectors ) = -m_Mmatrix.topLeftCorner ( m_numvectors, m_numvectors );
    }
    else if ( ExchField == 2 )
    {
        for ( int i = 0; i < m_numvectors; i++)
        {
            for ( int j = 0; j < m_numvectors ; j++)
            {
                if ( abs (m_kxl) < 1e-5 && abs (m_kyl) < 1e-5 && abs (m_Gvx.at(j)) < 0.001 / m_Crystal.a && abs (m_Gvy.at(j)) < 0.001 / m_Crystal.a )
                {
                    m_Mmatrix.topLeftCorner    ( m_numvectors, m_numvectors ) ( i, j ) = 0;

                    m_Mmatrix.topRightCorner   ( m_numvectors, m_numvectors ) ( i, j ) = m_Id[i][j] -
                                                                                         ComputeSlabDemagField ( i, j ) +
                                                                                         m_Qarray[i][j] *
                                                                                         ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(i) ) * ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) +
                                                                                           ( ky(m_kyl, m_Crystal) + m_Gvy.at(i) ) * ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) );
                                                                                        /*- BoundaryTermArray [i][j]*/

                    m_Mmatrix.bottomLeftCorner ( m_numvectors, m_numvectors ) ( i, j ) = -m_Mmatrix.topRightCorner ( m_numvectors, m_numvectors ) ( i, j );
                }
                else
                {
                    m_Mmatrix.topLeftCorner    ( m_numvectors, m_numvectors ) ( i, j ) = ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) *
                                                                                           ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) /
                                                                                           ( pow ( norm ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), ky(m_kyl, m_Crystal) + m_Gvy.at(j) ), 2 ) ) ) *
                                                                                         m_Msarray[i][j] * DipolarField;

                    m_Mmatrix.topRightCorner   ( m_numvectors, m_numvectors ) ( i, j ) = m_Id[i][j] -
                                                                                         ComputeSlabDemagField ( i, j ) +
                                                                                         pow ( ky(m_kyl, m_Crystal) + m_Gvy.at(j), 2 ) /
                                                                                         ( pow ( norm ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), ky(m_kyl, m_Crystal) + m_Gvy.at(j) ), 2 ) ) *
                                                                                         m_Msarray[i][j] * DipolarField +
                                                                                         m_Qarray[i][j] *
                                                                                         ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(i) ) * ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) +
                                                                                           ( ky(m_kyl, m_Crystal) + m_Gvy.at(i) ) * ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) );
                                                                                         /*- BoundaryTermArray [i][j]*/

                    m_Mmatrix.bottomLeftCorner ( m_numvectors, m_numvectors ) ( i, j ) = -m_Id[i][j] +
                                                                                         ComputeSlabDemagField ( i, j ) -
                                                                                         pow ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), 2 ) /
                                                                                         ( pow ( norm ( kx(m_kxl, m_Crystal) + m_Gvx.at(j), ky(m_kyl, m_Crystal) + m_Gvy.at(j) ), 2 ) ) *
                                                                                         m_Msarray[i][j] * DipolarField -
                                                                                         m_Qarray[i][j] *
                                                                                         ( ( kx(m_kxl, m_Crystal) + m_Gvx.at(i) ) * ( kx(m_kxl, m_Crystal) + m_Gvx.at(j) ) +
                                                                                           ( ky(m_kyl, m_Crystal) + m_Gvy.at(i) ) * ( ky(m_kyl, m_Crystal) + m_Gvy.at(j) ) );
                                                                                         /*+ BoundaryTermArray [i][j]*/
                }

                //if      ( ElectricField == 0 ) DMterm = 0;
                //else if ( ElectricField == 1 ) DMterm = ComputeDMterm ( i, j, ExchField );

                m_MmatrixComplex.topLeftCorner     ( m_numvectors, m_numvectors ) ( i, j ) = -m_Crystal.DMTerm;
                m_MmatrixComplex.topRightCorner    ( m_numvectors, m_numvectors ) ( i, j ) = 0;
                m_MmatrixComplex.bottomLeftCorner  ( m_numvectors, m_numvectors ) ( i, j ) = 0;
                m_MmatrixComplex.bottomRightCorner ( m_numvectors, m_numvectors ) ( i, j ) = -m_Crystal.DMTerm;
            }
        }
        m_Mmatrix.bottomRightCorner ( m_numvectors, m_numvectors ) = -m_Mmatrix.topLeftCorner ( m_numvectors, m_numvectors );
    }
}

void LandauLifshitzGilbert::Solve ( double kxl, double kyl, double kzl, int ExchField, int DipolarField, int ElectricField )  //solves LLG equation with solutions in evals and evecs
{
    m_kxl = kxl;
    m_kyl = kyl;
    m_kzl = kzl;
    CreateMmatrix ( ExchField, DipolarField, ElectricField );

    if ( ElectricField == 0 )
    {
        Eigen::EigenSolver <Eigen::MatrixXd> eigensolver;
        eigensolver.compute ( m_alphaMatrixInverse * m_Mmatrix, true );
        if (eigensolver.info() != Eigen::Success)
        {
            std::cerr << "Eigensolver failed.";
            exit ( EXIT_FAILURE );
        }
        m_evals = eigensolver.eigenvalues();
        m_evecs = eigensolver.eigenvectors();
        BubbleSort( ElectricField );
    }

    else if ( ElectricField == 1 )
    {
        Eigen::ComplexEigenSolver <Eigen::MatrixXcd> eigensolver;
        m_alphaMatrixInverseComplex.imag().setZero();
        m_alphaMatrixInverseComplex.real() = m_alphaMatrixInverse;
        m_MmatrixComplex.real() = m_Mmatrix;
        eigensolver.compute ( m_alphaMatrixInverseComplex * m_MmatrixComplex, true );
        if (eigensolver.info() != Eigen::Success)
        {
            std::cerr << "Eigensolver failed.";
            exit ( EXIT_FAILURE );
        }
        m_evals = eigensolver.eigenvalues();
        m_evecs = eigensolver.eigenvectors();
        BubbleSort( ElectricField );
    }
}

void LandauLifshitzGilbert::Solve_NoVectors ( double kxl, double kyl, double kzl, int ExchField, int DipolarField, int ElectricField )  //solves LLG equation with solutions in evals
{
    m_kxl = kxl;
    m_kyl = kyl;
    m_kzl = kzl;
    CreateMmatrix ( ExchField, DipolarField, ElectricField );

    if ( ElectricField == 0 )
    {
        Eigen::EigenSolver <Eigen::MatrixXd> eigensolver;
        m_Mmatrix = m_alphaMatrixInverse * m_Mmatrix;
        eigensolver.compute ( m_Mmatrix, false );
        if (eigensolver.info() != Eigen::Success)
        {
            std::cerr << "Eigensolver failed.";
            exit ( EXIT_FAILURE );
        }
        m_evals = eigensolver.eigenvalues();
        BubbleSort( ElectricField );
    }

    else if ( ElectricField == 1 )
    {
        Eigen::ComplexEigenSolver <Eigen::MatrixXcd> eigensolver;
        m_alphaMatrixInverseComplex.imag().setZero();
        m_alphaMatrixInverseComplex.real() = m_alphaMatrixInverse;
        m_MmatrixComplex.real() = m_Mmatrix;
        eigensolver.compute ( m_alphaMatrixInverseComplex * m_MmatrixComplex, false );
        if (eigensolver.info() != Eigen::Success)
        {
            std::cerr << "Eigensolver failed.";
            exit ( EXIT_FAILURE );
        }
        m_evals = eigensolver.eigenvalues();
        BubbleSort( ElectricField );
    }

//            Eigen::FullPivLU <Eigen::MatrixXd> LU(m_Mmatrix);
//            if ( LU.rank() != 2 * pow ( 2 * m_sPrecision.Nmax + 1, 2 ) )
//            {
//                std::cerr << "Defective Matrix at kx=" << m_kxl << ", ky =" << m_kyl;
//                exit ( EXIT_FAILURE );
//            }
}

complex <double> LandauLifshitzGilbert::ComputeFields ( Fields eFieldName, double x, double y, int mode )
{
    complex <double> sum(0,0), numerator, denominator, z,
                     z1(0,0), z2(0,0),
                     expTerm, CompTerm;
    double ExchL, RealTerm, Q;

    switch ( eFieldName )
    {
        case SATURATION_MAGNETIZATION:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                sum += FT_Variable ( m_Mat_A.Ms, m_Mat_B.Ms, i, -1 ) * exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) );
            }
            return sum;

        case EXCHANGE_LENGTH:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                sum += FT_Variable ( m_Mat_A.ExchLength, m_Mat_B.ExchLength, i, -1 ) * exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) );
            }
            return sum;

        case GILBERT_DAMPING:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                sum += FT_Variable ( m_Mat_A.alpha, m_Mat_B.alpha, i, -1 ) * exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) );
            }
            return sum;

        case MAGNETOSTATIC_POTENTIAL:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                if ( m_Gvx.at(i) == 0 && m_Gvy.at(i) == 0 ) continue;
                numerator = m_evecs ( i, m_evalsSortedIndex.at(mode) ) *                ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) +
                            m_evecs ( i + m_numvectors, m_evalsSortedIndex.at(mode) ) * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) );
                denominator = iii * ( ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) * ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) +
                                      ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) );
                z = ( numerator / denominator ) *
                    exp ( iii * ( ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) * x +
                                  ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) * y ) );
                sum += z;
            }
            return sum;

        case EXCHANGE_x:  //split into two terms, with derivative acting on either M_0 or m.
            for ( int i = 0; i < m_numvectors; i++ )
            {
                for ( int j = 0; j < m_numvectors; j++ )
                {
                    RealTerm = pow ( norm ( m_Gvx.at(j) + kx(m_kxl, m_Crystal), m_Gvy.at(j) + ky(m_kyl, m_Crystal) ), 2 );

                    ExchL    = FT_Variable ( 2 * m_Mat_A.A / ( m_Crystal.mu0 * m_Mat_A.Ms * m_Mat_A.Ms ), 2 * m_Mat_B.A / ( m_Crystal.mu0 * m_Mat_B.Ms * m_Mat_B.Ms ), i, -1 );
                    expTerm  = exp ( iii * ( ( m_Gvx.at(i) + m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) * x +
                                             ( m_Gvy.at(i) + m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) * y ) );
                    CompTerm = expTerm * m_evecs ( j, m_evalsSortedIndex.at(mode) );
                    z        = -CompTerm * RealTerm * ExchL;
                    sum     += z;
                }
            }

            for ( int i = 0; i < m_numvectors; i++ )
            {
                for ( int j = 0; j < m_numvectors; j++ )
                {
                    RealTerm = m_Gvx.at(i) * ( m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) + m_Gvy.at(i) * ( m_Gvy.at(j) + ky(m_kyl, m_Crystal) );

                    ExchL    = FT_Variable ( 2 * m_Mat_A.A / ( m_Crystal.mu0 * m_Mat_A.Ms * m_Mat_A.Ms ), 2 * m_Mat_B.A / ( m_Crystal.mu0 * m_Mat_B.Ms * m_Mat_B.Ms ), i, -1 );
                    expTerm  = exp ( iii * ( ( m_Gvx.at(i) + m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) * x +
                                             ( m_Gvy.at(i) + m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) * y ) );
                    CompTerm = expTerm * m_evecs ( j, m_evalsSortedIndex.at(mode) );
                    z        = -CompTerm * RealTerm * ExchL;
                    sum     += z;
                }
            }
            return sum;

        case EXCHANGE_y:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                for ( int j = 0; j < m_numvectors; j++ )
                {
                    RealTerm = pow ( norm ( m_Gvx.at(j) + kx(m_kxl, m_Crystal), m_Gvy.at(j) + ky(m_kyl, m_Crystal) ), 2 );

                    ExchL    = FT_Variable ( 2 * m_Mat_A.A / ( m_Crystal.mu0 * m_Mat_A.Ms * m_Mat_A.Ms ), 2 * m_Mat_B.A / ( m_Crystal.mu0 * m_Mat_B.Ms * m_Mat_B.Ms ), i, -1 );
                    expTerm  = exp ( iii * ( ( m_Gvx.at(i) + m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) * x +
                                             ( m_Gvy.at(i) + m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) * y ) );
                    CompTerm = expTerm * m_evecs ( j + m_numvectors, m_evalsSortedIndex.at(mode) );
                    z        = -CompTerm * RealTerm * ExchL;
                    sum     += z;
                }
            }

            for ( int i = 0; i < m_numvectors; i++ )
            {
                for ( int j = 0; j < m_numvectors; j++ )
                {
                    RealTerm = m_Gvx.at(i) * ( m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) + m_Gvy.at(i) * ( m_Gvy.at(j) + ky(m_kyl, m_Crystal) );

                    ExchL    = FT_Variable ( 2 * m_Mat_A.A / ( m_Crystal.mu0 * m_Mat_A.Ms * m_Mat_A.Ms ), 2 * m_Mat_B.A / ( m_Crystal.mu0 * m_Mat_B.Ms * m_Mat_B.Ms ), i, -1 );
                    expTerm  = exp ( iii * ( ( m_Gvx.at(i) + m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) * x +
                                             ( m_Gvy.at(i) + m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) * y ) );
                    CompTerm = expTerm * m_evecs ( j + m_numvectors, m_evalsSortedIndex.at(mode) );
                    z        = -CompTerm * RealTerm * ExchL;
                    sum     += z;
                }
            }
            return sum;

        case EXCHANGE_x_ALT:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                for ( int j = 0; j < m_numvectors; j++ )
                {
                    RealTerm = m_Gvx.at(i) * ( m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) + m_Gvy.at(i) * ( m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) + pow (  norm ( m_Gvx.at(j) + kx(m_kxl, m_Crystal), m_Gvy.at(j) + ky(m_kyl, m_Crystal) ), 2 );
                    Q        = FT_Variable ( m_Mat_A.Q, m_Mat_B.Q, i, -1 );
                    expTerm  = exp ( iii * ( ( m_Gvx.at(i) + m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) * x +
                                             ( m_Gvy.at(i) + m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) * y ) );
                    CompTerm = expTerm * m_evecs ( j, m_evalsSortedIndex.at(mode) );
                    z        = -CompTerm * RealTerm * Q;
                    sum     += z;
                }
            }
            return sum;

        case EXCHANGE_y_ALT:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                for ( int j = 0; j < m_numvectors; j++ )
                {
                    RealTerm = m_Gvx.at(i) * ( m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) + m_Gvy.at(i) * ( m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) + pow (  norm ( m_Gvx.at(j) + kx(m_kxl, m_Crystal), m_Gvy.at(j) + ky(m_kyl, m_Crystal) ), 2 );
                    Q        = FT_Variable ( m_Mat_A.Q, m_Mat_B.Q, i, -1 );
                    expTerm  = exp ( iii * ( ( m_Gvx.at(i) + m_Gvx.at(j) + kx(m_kxl, m_Crystal) ) * x +
                                             ( m_Gvy.at(i) + m_Gvy.at(j) + ky(m_kyl, m_Crystal) ) * y ) );
                    CompTerm = expTerm * m_evecs ( j + m_numvectors, m_evalsSortedIndex.at(mode) );
                    z        = -CompTerm * RealTerm * Q;
                    sum     += z;
                }
            }
            return sum;

        case DIPOLAR_x:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                numerator   = (  -m_Gvx.at(i) - kx( m_kxl, m_Crystal ) ) *
                              ( ( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * m_evecs ( i, m_evalsSortedIndex.at(mode) ) +
                                ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * m_evecs ( i + m_numvectors, m_evalsSortedIndex.at(mode) ) );
                denominator  = -( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) -
                                ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) );
                z = exp ( iii * ( ( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * x +
                                  ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * y ) ) *
                    numerator / denominator;
                sum += z;
            }
            return sum;

        case DIPOLAR_y:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                numerator   = (  -m_Gvy.at(i) - ky( m_kyl, m_Crystal ) ) *
                              ( ( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * m_evecs ( i, m_evalsSortedIndex.at(mode) ) +
                                ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * m_evecs ( i + m_numvectors, m_evalsSortedIndex.at(mode) ) );
                denominator  = -( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) -
                                ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) );
                z = exp ( iii * ( ( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * x +
                                  ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * y ) ) *
                    numerator / denominator;
                sum += z;
            }
            return sum;

        case CURL_OF_DIPOLAR_xWRTy:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                numerator = ( ( -m_Gvx.at(i) - kx(m_kxl, m_Crystal) ) * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) ) *
                            ( (  m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) * m_evecs ( i, m_evalsSortedIndex.at(mode) ) +
                              (  m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) * m_evecs ( i + m_numvectors, m_evalsSortedIndex.at(mode) ) );
                denominator = -( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) * ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) - ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) );
                z = exp ( iii * ( ( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * x +
                                  ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * y ) ) *
                    iii * numerator / denominator;
                sum += z;
            }
            return sum;

        case CURL_OF_DIPOLAR_yWRTx:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                numerator = ( ( -m_Gvy.at(i) - ky(m_kyl, m_Crystal) ) * ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) ) *
                            ( (  m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) * m_evecs ( i, m_evalsSortedIndex.at(mode) ) +
                              (  m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) * m_evecs ( i + m_numvectors, m_evalsSortedIndex.at(mode) ) );
                denominator = -( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) * ( m_Gvx.at(i) + kx(m_kxl, m_Crystal) ) - ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) ) * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal) );
                z = exp ( iii * ( ( m_Gvx.at(i) + kx( m_kxl, m_Crystal ) ) * x +
                                  ( m_Gvy.at(i) + ky( m_kyl, m_Crystal ) ) * y ) ) *
                    iii * numerator / denominator;
                sum += z;
            }
            return sum;

        case MAGNETIZATION_x:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                z = m_evecs ( i, m_evalsSortedIndex.at(mode) ) * exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) );
                sum += z;
            }
            return sum;

        case MAGNETIZATION_y:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                z = m_evecs ( i + m_numvectors, m_evalsSortedIndex.at(mode) ) * exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) );
                sum += z;
            }
            return sum;

        case CURL_OF_MAGNETIZATION_xWRTy:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                z2 = m_evecs ( i, m_evalsSortedIndex.at(mode) ) * exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) );
                z = z2 * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal ) );
                sum += z;
            }
            return sum;

        case CURL_OF_MAGNETIZATION_yWRTx:
            for ( int i = 0; i < m_numvectors; i++ )
            {
                z2 = m_evecs ( i + m_numvectors, m_evalsSortedIndex.at(mode) ) * exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) );
                z = z2 * ( m_Gvy.at(i) + ky(m_kyl, m_Crystal ) );
                sum += z;
            }
            return sum;

        case DEMAG_FIELD:
            for ( int i = 0; i< m_numvectors; i++ )
            {
                sum += exp ( iii * ( m_Gvx.at(i) * x + m_Gvy.at(i) * y ) ) *
                       -4. * m_Crystal.Pi * FT_Variable ( m_Mat_A.Ms, m_Mat_B.Ms, i, -1 ) /
                       ( cosh ( m_Crystal.SlabWidth * norm ( ( m_Gvx.at(i) ), ( m_Gvy.at(i) ) ) ) + sinh ( m_Crystal.SlabWidth * norm ( ( m_Gvx.at(i) ), ( m_Gvy.at(i) ) ) ) );
            }
            return sum;

        default:
            cerr << "Invalid FieldName in function ComputeFields.";
            exit ( EXIT_FAILURE );
    }
}

complex <double> LandauLifshitzGilbert::ComputeFields_kdotx (  Fields eFieldName, double x, double y, int mode )
{
    return ComputeFields (  eFieldName, x, y, mode ) *
           exp ( iii * ( kx(m_kxl, m_Crystal) * x + ky(m_kyl, m_Crystal) * y ) );
}
