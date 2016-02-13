#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <string>
#include <sstream>

//gsl headers
//#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_eigen.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_complex_math.h>

//local headers
#include "constants.h"
#include "declarations.h"
#include "LLG.h"
#include "FieldPlots.h"
#include "Faddeeva.h"

//choose the shape of the rod's cross-section
#define CIRCLE
//#define ELLIPSE

//choose interface type ( diffuse currently only works with circular cross-sections )
#define IDEAL
//#define DIFFUSE

//uncomment for slab geometry
//#define SLAB

//uncomment to output field plot data to FieldPlots.csv instead of calculating GF/response
#define DISPERSIONS       //output file:returns 23
//#define MAG_VEC_SYMMETRY  //output file:returns 24
//#define FIELD_PLOTS       //output file: returns 25
//uncomment to calculate integrated h and H_ex as a function of some parameter
//#define INTEGRATED_FIELDS  //FIELD_PLOTS must also be defined to use this

//use one of these to calculate either GF images or pulse movies
//#define GREEN_FUNCTION
//#define PULSE
//#define DENSITYOFSTATES
//#define ADD_DENSITYOFSTATES

//use one or the other to choose the lattice type
#define SQUARE_LATTICE    //******* Change nummodes when changing lattice type!! *********
//#define HEX_LATTICE    //currently unknown if sigma factor and lattice vector limiting of Ms, etc. works for hex

//list of files in data folder:
//FieldPlots.csv - plots for m, psi, h, H_ex, curl of m, and curl of h
//Integrated_Fields.csv - plots for integrated h and H_ex vs some parameter
//Imaginary_Frequencies.csv - Re(omega)
//Real_Frequencies.csv - Im(oemga)
//FOM.csv - Re(omega) / Im(omega)
//MagVecSymm.csv - magnetization vector and corresponding symmetry plots

using namespace std;

//Required Folders:  Contours, GreenFxn, Info, misc, Progress, Data, Data/DispersionData, Integrated_Fields. For summing: Pyxplot, Sum, VTK

#ifdef SQUARE_LATTICE
int const lattice = 1;
#endif

#ifdef HEX_LATTICE
int const lattice = 2;
#endif

complex <double> const iii ( 0, 1 );

int main()
{
    #include "variables.h"

    double kxlstart, kxlend, kylstart, kylend;
    //Check for lattice type (1=square, 2=hexagonal):
    if ( lattice == 1 )
    {
        kxlstart = -1 + sPrecision.stepsize; kylstart = -1 + sPrecision.stepsize;
        kxlend   =  1.0; kylend   =  1.0;
    }
    else if ( lattice == 2 )
    {
        kxlstart = -2. / sqrt(3.); kylstart = -1.;
        kxlend   =  2. / sqrt(3.); kylend   =  1.;
    }

    int PartNum, kTotal = ceil ( ( ( ( kxlend - kxlstart ) / sPrecision.stepsize + 1 ) * ( ( kylend - kylstart ) / sPrecision.stepsize + 1 ) ) / ( processes + 1 ) ), kMin, kMax;

//Check for existence of info files to determine which part number to run:
    for ( PartNum = 0; PartNum <= floor ( pow ( 2 / sPrecision.stepsize + 1, 2 ) / kTotal ) ; PartNum++ )
    {
        FILE *InfoFile = fopen ( NameFILE ( "Info/Info", sPrecision.Nmax, sPrecision.stepsize, PartNum, ".dat" ).c_str(), "r" ); //replaced sPrecision.Nmax with 0 temporarily
        if ( InfoFile == NULL ) break;
    }
    if ( PartNum > floor ( pow ( 2 / sPrecision.stepsize + 1, 2 ) / kTotal ) ) return 0;

    kMin = PartNum * kTotal;
    kMax = ( PartNum + 1 ) * kTotal;

    ofstream Info ( NameFILE ( "Info/Info", sPrecision.Nmax, sPrecision.stepsize, PartNum, ".dat" ).c_str(), ios::out | ios::trunc );  //replaced sPrecision.Nmax with 0 temporarily
    Info.close();

//  testing the read/write funcions for Eigen
//    Eigen::MatrixXcd WriteTest, ReadTest;
//    WriteTest.resize ( 2, 2 ); ReadTest.resize ( 2, 2 );
//    WriteTest ( 0, 0 ) = 1.+iii;
//    WriteTest ( 0, 1 ) = 1.-iii;
//    WriteTest ( 1, 0 ) = -1.+iii;
//    WriteTest ( 1, 1 ) = -1.-iii;
//
//    WriteMatrix ( "Test.bin", WriteTest, 2, 2 );
//    ReadMatrix  ( "Test.bin", ReadTest, 2, 2 );
//    cout << WriteTest << endl;
//    cout << ReadTest << endl;

//    Response cResponse ( Mat_A, Mat_B, Crystal, sPrecision, sPlotData );
//    cResponse.SetPositions();
//    cResponse.CalculateResponse ( PartNum, kxlstart, kylstart, kxlend, kylend, kMin, kMax, ExchField, DipolarField, ElectricField, FreqGF, kTotal );
//    cResponse.WriteToFile ( PartNum );


//****************************************Beginning of Calculations****************************************

//    Basic_Setup ( Mat_A, Mat_B, Crystal, N_Interface, numvectors, Nmax );  //creates Gvectors, matrices independent of k, and alpha matrix inverse

    #ifdef DISPERSIONS
//    int increment = floor ( PartNum / 10 );
//    sPrecision.Nmax = 8 + ( PartNum - 10 * increment );
//    Mat_B.alpha = 0.0006 * pow ( 10, -( 1 + increment ) );

    LandauLifshitzGilbert cLLG ( Mat_A, Mat_B, Crystal, sPrecision );
    Dispersions ( cLLG, ExchField, DipolarField, ElectricField, PartNum );

    //testing convergence of a single k point:
//    cLLG.Solve( 1.0, 0.4, ExchField, DipolarField );
//
//    std::ofstream ConvFile  ( NameFILE ( "Data/ConvergenceTest_2_", 0, 0, 0, ".csv"   ).c_str(), std::ios::out | std::ios::app );
//    ConvFile.setf  ( std::ios::scientific );
//    ConvFile.precision  ( 15 );
//
//    if ( PartNum == 0 )
//    {
//        ConvFile << "Col. 1:   PartNum" << endl;
//        ConvFile << "Col. 2:   Re (Omega)" << endl;
//        ConvFile << "Col. 3:   Im (Omega)" << endl;
//        ConvFile << "Col. 4:   Re (Omega) / Im (Omega)" << endl;
//        ConvFile << endl << endl << endl;  //creates a new index for pyxplot
//    }
//
//    ConvFile << sPrecision.Nmax << " ";
//    for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
//    {
//        ConvFile << GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) ) << " ";
//        ConvFile << GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) ) << " ";
//        ConvFile << GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) ) / GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) ) << " ";
//    }
//    ConvFile << std::endl;
//    ConvFile.close();

    return 23;
    #endif

    #ifdef MAG_VEC_SYMMETRY
    LandauLifshitzGilbert cLLG ( Mat_A, Mat_B, Crystal, sPrecision );
    MagVecSymm ( cLLG, ExchField, DipolarField );
    return 24;
    #endif

    #ifdef FIELD_PLOTS
    double kxl = 1e-50, kyl = 1e-50, kzl = 0.0;
    FieldPlots cFieldPlots ( Mat_A, Mat_B, Crystal, sPrecision );
    cFieldPlots.PlotFT_Variables();
//    cFieldPlots.Solve ( kxl, kyl, kzl, ExchField, DipolarField );
//    cFieldPlots.Plot ();
    return 25;
    #endif

    #ifdef INTEGRATED_FIELDS
    //sPrecision.Nmax = 6 + PartNum * 1;
    sPrecision.spatial_sizePrime = 100; //+ PartNum * 25;

    double kxl = 1e-50, kyl = 1e-50, kzl = 0.0;
    FieldPlots cIntegratedFields ( Mat_A, Mat_B, Crystal, sPrecision );
    cIntegratedFields.Solve ( kxl, kyl, kzl, ExchField, DipolarField );
    cIntegratedFields.PlotIntegrated ( PartNum );
    return 26;
    #endif

    #ifdef DENSITYOFSTATES
    DensityOfStates cDensityOfStates ( Mat_A, Mat_B, Crystal, sPrecision );
    cDensityOfStates.CalculateDoS( PartNum, kxlstart, kylstart, kxlend, kylend, kMin, kMax, ExchField, DipolarField, ElectricField, kTotal );
    #endif

    #ifdef ADD_DENSITYOFSTATES
    DensityOfStates cDensityOfStates ( Mat_A, Mat_B, Crystal, sPrecision );
    cDensityOfStates.AddDoSMatrix();
    #endif
}
