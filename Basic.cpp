#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

//gsl headers
//#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_complex_math.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_eigen.h>

//local headers
#include "constants.h"
#include "declarations.h"
#include "Faddeeva.h"

using namespace std;
using namespace Eigen;

string NameFILE ( string FileName, int Nmax, double stepsize, int PartNum, string Extension )
{
    stringstream ssFILENAME;
    if ( PartNum < 10 ) ssFILENAME <<
            FileName << Nmax << "_" << stepsize << "_Part00" << PartNum << Extension;
    else if ( PartNum < 100 ) ssFILENAME <<
            FileName << Nmax << "_" << stepsize << "_Part0" << PartNum << Extension;
    else ssFILENAME <<
            FileName << Nmax << "_" << stepsize << "_Part" << PartNum << Extension;

    string FILENAMEString = ssFILENAME.str();
    return FILENAMEString;
}

double kx ( double length, Structure Crystal )
{
    if ( lattice == 1 )
    {
        return length * Crystal.Pi / Crystal.a;
    }
    else if ( lattice == 2 )
    {
        return length * 2 * Crystal.Pi / ( Crystal.a * sqrt(3) );
    }
    else return 0;
}

double ky ( double length, Structure Crystal )
{
    if ( lattice == 1 )
    {
        return length * Crystal.Pi / Crystal.a;
    }
    else if ( lattice == 2 )
    {
        return length * 2 * Crystal.Pi / ( Crystal.a * sqrt(3) );
    }
    else return 0;
}

double kz ( double length, Structure Crystal )
{
    if ( lattice == 1 )
    {
//        return length * Crystal.Pi / Crystal.a;
        return length;
    }
    else if ( lattice == 2 )
    {
//        return length * 2 * Crystal.Pi / ( Crystal.a * sqrt(3) );
        return length;
    }
    else return 0;
}

double rotatexcomp ( double xcomp, double ycomp, double angle )
{
    return xcomp * cos ( angle ) - ycomp * sin ( angle );
}

double rotateycomp ( double xcomp, double ycomp, double angle )
{
    return ycomp * cos ( angle ) + xcomp * sin ( angle );
}

double norm ( double xcomp, double ycomp )
{
    return sqrt ( xcomp * xcomp + ycomp * ycomp );
}

void ChoosePath ( LandauLifshitzGilbert cLLG, int ExchField, int DipolarField, int ElectricField, int PartNum, std::vector < std::vector <double> >& dispersiondata, std::vector <double>& posdata, int &temp, int nPath, double PrevPosition )
{
    std::ofstream DispersionProgressFile ( NameFILE ( "Info/Progress_", 0, 0, PartNum, ".dat" ).c_str(), std::ios::out | std::ios::trunc );
    std::vector <double> DispersionSubset ( 2 * cLLG.Get_Precision().nummodes, 0 );

    switch ( nPath )
    {
        case 1:  //Gamma to M along diagonal:  ( 0, 0 ) to ( 1, 1 )
            for ( double i = cLLG.Get_Precision().stepsize; i < sqrt(2); i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( i / sqrt(2), i / sqrt(2), 0.0 * cLLG.Get_Crystal().Pi / ( 2 * cLLG.Get_Crystal().SlabWidth ), ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                std::cout << posdata.at ( temp - 1 ) << std::endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << std::endl;
            }
            break;
        case 2:  //M to X: ( 1, 1 ) to ( 1, 0 )
            for ( double i = 0; i < 1; i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 1, 1 - i, 0.0 * cLLG.Get_Crystal().Pi / ( 2 * cLLG.Get_Crystal().SlabWidth ), ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                std::cout << posdata.at ( temp - 1 ) << std::endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << std::endl;
            }
            break;
        case 3:  //X to Gamma: ( 1, 0 ) to ( 0, 0 )
            for ( double i = 0; i < 1; i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 1 - i, 1e-5, 0.0 * cLLG.Get_Crystal().Pi / ( 2 * cLLG.Get_Crystal().SlabWidth ), ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                std::cout << posdata.at ( temp - 1 ) << std::endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << std::endl;
            }
            break;
        case 4: //Gamma to X':  ( 0, 0 ) to ( 0, 1 )
            for ( double i = cLLG.Get_Precision().stepsize; i < 1; i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 1e-5, i, 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                std::cout << posdata.at ( temp - 1 ) << std::endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << std::endl;
            }
            break;
        case 5: // M to X: ( 0, 1 ) to ( 1, 1 )
            for ( double i = 0; i < 1; i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( i, 1, 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                std::cout << posdata.at ( temp - 1 ) << std::endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << std::endl;
            }
            break;
        case 6: // Hex: Gamma to M: ( 0, 0 ) to ( 0, 1 )
            for ( double i = cLLG.Get_Precision().stepsize; i < 1; i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 0, i, 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                cout << posdata.at ( temp - 1 ) << endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << endl;
            }
            break;
        case 7: // Hex: M to K: ( 0, 1 ) to ( 1 / sqrt(3), 1 )
            for ( double i = 0; i < 1 / sqrt(3); i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( i, 1, 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                cout << posdata.at ( temp - 1 ) << endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << endl;
            }
            break;
        case 8: // Hex: K to Gamma: ( 1 / sqrt(3), 1 ) to ( 0, 0 )
            for ( double i = 0; i < 2 / sqrt(3); i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 0.5 * ( 2 / sqrt(3) - i ), 0.5 * sqrt(3) * ( 2 / sqrt(3) - i ), 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                cout << posdata.at ( temp - 1 ) << endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << endl;
            }
            break;
        case 9: // Hex: K to K': ( 1 / sqrt*(3), 1 ) to ( 2 / sqrt(3), 0 )
            for ( double i = 0; i < 2 / sqrt(3); i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 0.5 * ( 2 / sqrt(3) - i ), 0.5 * sqrt(3) * ( 2 / sqrt(3) - i ), 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                cout << posdata.at ( temp - 1 ) << endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << endl;
            }
            break;
        case 10: // Hex: K to Gamma: ( 2 / sqrt*(3), 0 ) to ( 0, 0 )
            for ( double i = 0; i < 2 / sqrt(3); i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 0.5 * ( 2 / sqrt(3) - i ), 0.5 * sqrt(3) * ( 2 / sqrt(3) - i ), 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                cout << posdata.at ( temp - 1 ) << endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << endl;
            }
            break;
        case 11: // Hex: Gamma to M': ( 0, 0 ) to (  )
            for ( double i = 0; i < 2 / sqrt(3); i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 0.5 * ( 2 / sqrt(3) - i ), 0.5 * sqrt(3) * ( 2 / sqrt(3) - i ), 0.0, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                cout << posdata.at ( temp - 1 ) << endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << endl;
            }
            break;
        case 12:  //z-direction:  ( 0, 0, 0 ) to ( 0, 0, 1 )
            for ( double i = cLLG.Get_Precision().stepsize; i < 3; i += cLLG.Get_Precision().stepsize )
            {
                dispersiondata.push_back ( DispersionSubset );

                cLLG.Solve_NoVectors( 1.0, 1e-5, i, ExchField, DipolarField, ElectricField );
                for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
                {
                    dispersiondata [temp] [ 2 * j ]     = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).real();
                    dispersiondata [temp] [ 1 + 2 * j ] = cLLG.m_evals ( cLLG.m_evalsSortedIndex.at(j) ).imag();
                }
                posdata.push_back ( i + PrevPosition );
                temp++;
                std::cout << posdata.at ( temp - 1 ) << std::endl;
                DispersionProgressFile << posdata.at ( temp - 1 ) << std::endl;
            }
            break;
        default:
            cerr << "Invalid path number in function ChoosePath.";
            exit ( EXIT_FAILURE );
    }
    DispersionProgressFile.close();
}


void Dispersions ( LandauLifshitzGilbert cLLG, int ExchField, int DipolarField, int ElectricField, int PartNum )
{
    //Calculate and output dispersion data:
//    int numpoints = 0;
//
//    if (  lattice == 1 )
//    {
//        numpoints = ( 2 * ceil ( 1. / cLLG.Get_Precision().stepsize ) + ceil ( sqrt(2) / cLLG.Get_Precision().stepsize ) );  //numpoints is 4 times as large as necessary for circular cross-sections to accomadate ellipses
//    }
//
//    if ( lattice == 2 )
//    {
//        numpoints = (int) ( 1. / cLLG.Get_Precision().stepsize ) + (int) ( 1 / ( sqrt(3) * cLLG.Get_Precision().stepsize ) ) + (int) ( 2 / ( sqrt(3) * cLLG.Get_Precision().stepsize ) ) + 1;
//    }

    int temp = 0;
    std::vector < std::vector <double> > dispersiondata;
//    dispersiondata.resize ( 0, vector <double> ( 2 * cLLG.Get_Precision().nummodes, 0 ) );
    std::vector <double> posdata;
    //2014-08-04: redid this section to accomadate for ellipses
//    double dispersiondata2 [ numpoints ] [ 2 * cLLG.Get_Precision().nummodes ];
//    double posdata2 [ numpoints ];

//    std::ofstream DispersionProgressFile ( NameFILE ( "Info/Progress_", 0, 0, PartNum, ".dat" ).c_str(), std::ios::out | std::ios::trunc );

    std::ofstream ImagFreq  ( NameFILE ( "Data/DispersionData/Imaginary_Frequencies_", 0, 0, 0, ".csv"   ).c_str(), std::ios::out | std::ios::trunc );
    std::ofstream RealFreq  ( NameFILE ( "Data/DispersionData/Real_Frequencies_",      0, 0, 0, ".csv"   ).c_str(), std::ios::out | std::ios::trunc );
    std::ofstream FOMfile   ( NameFILE ( "Data/DispersionData/FOM_",                   0, 0, 0, ".csv"   ).c_str(), std::ios::out | std::ios::trunc );
    std::ofstream DispersionFile   ( NameFILE ( "Data/DispersionData/DispersionData_", 0, 0, 0, ".csv"   ).c_str(), std::ios::out | std::ios::trunc );
    std::ofstream VisitFile ( NameFILE ( "Data/DispersionData/Visit_Data_",            0, 0, PartNum, ".curve" ).c_str(), std::ios::out | std::ios::trunc );
    ImagFreq.setf  ( std::ios::scientific );
    ImagFreq.precision  ( 15 );
    RealFreq.setf  ( std::ios::scientific );
    RealFreq.precision  ( 15 );
    FOMfile.setf   ( std::ios::scientific );
    FOMfile.precision   ( 15 );
    DispersionFile.setf   ( std::ios::scientific );
    DispersionFile.precision   ( 15 );
    VisitFile.setf ( std::ios::scientific );
    VisitFile.precision ( 15 );

    if ( lattice == 1 )
    {
//        for ( double i = cLLG.Get_Precision().stepsize; i < sqrt(2); i += cLLG.Get_Precision().stepsize )
//        {
//            cLLG.Solve( i / sqrt(2), i / sqrt(2), ExchField, DipolarField );
//            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
//            {
//                dispersiondata [temp] [ 2 * j ]     = GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//                dispersiondata [temp] [ 1 + 2 * j ] = GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//            }
//            posdata [temp] = i;
//            temp++;
//            std::cout << posdata[temp-1] << std::endl;
//            DispersionProgressFile << posdata[temp-1] << std::endl;
//        }
//        for ( double i = 0; i < 1; i += cLLG.Get_Precision().stepsize )
//        {
//            cLLG.Solve( 1, 1 - i, ExchField, DipolarField );
//            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
//            {
//                dispersiondata [temp] [ 2 * j ]     = GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//                dispersiondata [temp] [ 1 + 2 * j ] = GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//            }
//            posdata [temp] = i + sqrt(2);
//            temp++;
//            std::cout << posdata[temp-1] << std::endl;
//            DispersionProgressFile << posdata[temp-1] << std::endl;
//        }
//        for ( double i = 0; i < 1; i += cLLG.Get_Precision().stepsize )
//        {
//            cLLG.Solve( 1 - i, 0, ExchField, DipolarField );
//            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
//            {
//                dispersiondata [temp] [ 2 * j ]     = GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//                dispersiondata [temp] [ 1 + 2 * j ] = GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//            }
//            posdata [temp] = i + sqrt(2) + 1;
//            temp++;
//            std::cout << posdata[temp-1] << std::endl;
//            DispersionProgressFile << posdata[temp-1] << std::endl;
//        }

        switch ( cLLG.Get_Crystal().eShape )
        {
            case SHAPE_CIRCLE:
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 1, 0           );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 2, sqrt(2)     );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 3, sqrt(2) + 1 );
                //ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 12, sqrt(2) + 1 );
                break;

            case SHAPE_ELLIPSE:
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 1, 0           );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 2, sqrt(2)     );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 3, sqrt(2) + 1 );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 4, sqrt(2) + 2 );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 5, sqrt(2) + 3 );
                break;
            default:
                cerr << "Invalid shape in function Dispersions.";
                exit ( EXIT_FAILURE );
        }

        for ( int i = 0; i < temp; i++ )
        {
            ImagFreq << posdata[i] << " ";
            RealFreq << posdata[i] << " ";
            FOMfile  << posdata[i] << " ";
            DispersionFile << posdata[i] << " ";
            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
            {
                ImagFreq << dispersiondata [i][2 * j + 1] << " ";
                RealFreq << dispersiondata [i][2 * j] << " ";
                FOMfile  << dispersiondata [i][2 * j + 1] / dispersiondata [i][2 * j] << " ";
                DispersionFile << dispersiondata [i][2 * j + 1] << " " <<
                                  dispersiondata [i][2 * j]     << " ";
            }
            ImagFreq << std::endl;
            RealFreq << std::endl;
            FOMfile  << std::endl;
            DispersionFile << std::endl;
        }
        ImagFreq.close();
        RealFreq.close();
        FOMfile.close();
        DispersionFile.close();

        //output files for creating animations with visit:
        for ( int j = 0; j < 20; j++ )
        {
            if ( j < 9 ) VisitFile << "#Real/0" << j + 1 << endl;
            else VisitFile << "#Real/" << j + 1 << endl;
            for ( int i = 0; i < temp; i++ )
            {
                VisitFile << posdata[i] << " ";
                VisitFile << dispersiondata [i][2 * j + 1] << std::endl;
            }
            VisitFile << std::endl << std::endl;

            if ( j < 9 ) VisitFile << "#Imag/0" << j + 1 << endl;
            else VisitFile << "#Imag/" << j + 1 << endl;
            for ( int i = 0; i < temp; i++ )
            {
                VisitFile << posdata[i] << " ";
                VisitFile << dispersiondata [i][2 * j] << std::endl;
            }
            VisitFile << std::endl << std::endl;

            if ( j < 9 ) VisitFile << "#FOM/0" << j + 1 << endl;
            else VisitFile << "#FOM/" << j + 1 << endl;
            for ( int i = 0; i < temp; i++ )
            {
                VisitFile << posdata[i] << " ";
                VisitFile << dispersiondata [i][2 * j + 1] / dispersiondata [i][2 * j] << std::endl;
            }
            VisitFile << std::endl << std::endl;
        }
        VisitFile.close();
    }

    else if ( lattice == 2 )
    {
//        for ( double i = cLLG.Get_Precision().stepsize; i < 1; i += cLLG.Get_Precision().stepsize )
//        {
//            cLLG.Solve( 0, i, ExchField, DipolarField );
//            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
//            {
//                dispersiondata [temp] [2 * j] = GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//                dispersiondata [temp] [1 + 2 * j] = GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//            }
//            posdata [temp] = i;
//            temp++;
//            cout << posdata[temp-1] << endl;
//            DispersionProgressFile << posdata[temp-1] << endl;
//        }
//        for ( double i = 0; i < 1 / sqrt(3); i += cLLG.Get_Precision().stepsize )
//        {
//            cLLG.Solve( i, 1, ExchField, DipolarField );
//            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
//            {
//                dispersiondata [temp] [ 2 * j] = GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//                dispersiondata [temp] [1 + 2 * j] = GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//            }
//            posdata [temp] = i + 1;
//            temp++;
//            cout << posdata[temp-1] << endl;
//            DispersionProgressFile << posdata[temp-1] << endl;
//        }
//        for ( double i = 0; i < 2 / sqrt(3); i += cLLG.Get_Precision().stepsize )
//        {
//            cLLG.Solve( 0.5 * ( 2 / sqrt(3) - i ), 0.5 * sqrt(3) * ( 2 / sqrt(3) - i ), ExchField, DipolarField );
//            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
//            {
//                dispersiondata [temp] [ 2 * j] = GSL_REAL ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//                dispersiondata [temp] [1 + 2 * j] = GSL_IMAG ( gsl_vector_complex_get ( cLLG.m_evals, 2 * j ) );
//            }
//            posdata [temp] = i + 1 / sqrt(3) + 1;
//            temp++;
//            cout << posdata[temp-1] << endl;
//            DispersionProgressFile << posdata[temp-1] << endl;
//        }

        switch ( cLLG.Get_Crystal().eShape )
        {
            case SHAPE_CIRCLE:
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 6, 0               );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 7, 1               );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 8, 1 + 1 / sqrt(3) );
                break;

            case SHAPE_ELLIPSE:
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 8, 0               );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 6, 2 / sqrt(3)     );
                ChoosePath ( cLLG, ExchField, DipolarField, ElectricField, PartNum, dispersiondata, posdata, temp, 7, 2 / sqrt(3) + 1 );
                break;
            default:
                cerr << "Invalid shape in function Dispersions.";
                exit ( EXIT_FAILURE );
        }

        for ( int i =0; i < temp; i++ )
        {
            ImagFreq << posdata[i] << " ";
            RealFreq << posdata[i] << " ";
            FOMfile << posdata[i] << " ";
            DispersionFile << posdata[i] << " ";
            for ( int j = 0; j < cLLG.Get_Precision().nummodes; j++ )
            {
                ImagFreq << dispersiondata [i][2 * j + 1] << " ";
                RealFreq << dispersiondata [i][2 * j] << " ";
                FOMfile << dispersiondata [i][2 * j + 1] / dispersiondata [i][2 * j] << " ";
                DispersionFile << dispersiondata [i][2 * j + 1] << " " <<
                                  dispersiondata [i][2 * j]     << " ";
            }
            ImagFreq << endl;
            RealFreq << endl;
            FOMfile << endl;
            DispersionFile << endl;
        }
        ImagFreq.close();
        RealFreq.close();
        FOMfile.close();
        DispersionFile.close();

        //output files for creating animations with visit:
        for ( int j = 0; j < 20; j++ )
        {
            if ( j < 9 ) VisitFile << "#Real/0" << j + 1 << endl;
            else VisitFile << "#Real/" << j + 1 << endl;
            for ( int i = 0; i < temp; i++ )
            {
                VisitFile << posdata[i] << " ";
                VisitFile << dispersiondata [i][2 * j + 1] << std::endl;
            }
            VisitFile << std::endl << std::endl;

            if ( j < 9 ) VisitFile << "#Imag/0" << j + 1 << endl;
            else VisitFile << "#Imag/" << j + 1 << endl;
            for ( int i = 0; i < temp; i++ )
            {
                VisitFile << posdata[i] << " ";
                VisitFile << dispersiondata [i][2 * j] << std::endl;
            }
            VisitFile << std::endl << std::endl;

            if ( j < 9 ) VisitFile << "#FOM/0" << j + 1 << endl;
            else VisitFile << "#FOM/" << j + 1 << endl;
            for ( int i = 0; i < temp; i++ )
            {
                VisitFile << posdata[i] << " ";
                VisitFile << dispersiondata [i][2 * j + 1] / dispersiondata [i][2 * j] << std::endl;
            }
            VisitFile << std::endl << std::endl;
        }
        VisitFile.close();
    }
}

//This function describes the distribution of the source ( ex. Gaussian ):
double DistributionFunction ( double xprime, double yprime, Structure Crystal )
{
//    Use this part if using a delta function in space and time:
//    if ( abs (xprime) < a * 1e-5 && abs (yprime) < a * 1e-5 && abs (timeprime) < 1e-5 ) return 1;
//    else return 0;


    return exp ( -( xprime * xprime + yprime * yprime ) / ( 2 * Crystal.sigmaSpace * Crystal.sigmaSpace ) );
}

// This is the term obtained from the Fourier transform of the freq dependent Green's fxn and integral over tprime:
complex <double> FourierIntegral ( complex <double> omega, double time, int n, double EndTime, VectorXcd& OMEGA0, PlotData sPlotData, double Pi )
{
//    Use this part if computing the frequency dependendt Green's Fucntion:
//    Note (2012.09.05): Frequencies in denominator may have mismatched dimensions.  Added alternate.  Denominator may also contain sign error.
//    gsl_complex denominator;
//    denominator = gsl_complex_mul ( gsl_complex_rect ( gyroRatio * mu0 * H0, 0.0 ), gsl_complex_sub ( omega, gsl_complex_rect ( 5, 0.0 ) ) );
//    ALTERNATE:  denominator = gsl_complex_sub ( omega, OMEGA0 );
//    return gsl_complex_div ( gsl_complex_rect ( 1.0, 0.0 ), denominator );

    //Fourier transform:
    //Note (2012.09.04) - This may be incorrect as well.
//    double tFinal;
//
//    if ( time < EndTime ) tFinal = time;
//    else tFinal = EndTime;
//    return gsl_complex_mul ( gsl_complex_mul_real ( gsl_complex_inverse ( omega ), 2 * Crystal.Pi ), gsl_complex_sub ( gsl_complex_exp ( gsl_complex_mul_imag ( omega, -time + tFinal ) ), gsl_complex_exp ( gsl_complex_mul_imag ( omega, -time ) ) ) );


//  Use this for a source at a specific frequency (OMEGA0):
//  Note (2012.09.04) - This seems to be incorrect.  Trying something different below.  (Maybe just sign errors?)
//    double tFinal;
//
//    if ( time < EndTime ) tFinal = time;
//    else tFinal = EndTime;
//    return gsl_complex_mul ( gsl_complex_mul_real ( gsl_complex_inverse ( gsl_complex_sub ( omega, OMEGA0 ) ), 2 * Crystal.Pi ),
//           gsl_complex_sub ( gsl_complex_exp ( gsl_complex_add ( gsl_complex_mul_imag ( omega, -time + tFinal ), gsl_complex_mul_imag ( OMEGA0, -tFinal ) ) ),
//                             gsl_complex_exp ( gsl_complex_mul_imag ( omega, -time ) ) ) );


//  Use this for a source at a specific frequency (OMEGA0):
//    double tFinal;
//    complex <double> result;
//
//    if ( time < EndTime ) tFinal = time;
//    else tFinal = EndTime;
//    result = ( exp ( iii * omega * ( -time + tFinal ) - iii * OMEGA0(n) * tFinal ) -
//             exp ( -iii * omega * time ) ) /
//           ( OMEGA0(n) - omega ); // Note (20130131): changed signs on all time variables
//
//    return result;

//Gaussian distribution centered around tCenter:
    complex <double> arg_t, arg_0, preFactor_t, preFactor_0, result;

//    arg_t = ( pow ( sPlotData.sigmaTime, 2 ) * ( OMEGA0(n) - omega ) - iii * ( time - sPlotData.tCenter ) ) /
//            ( sPlotData.sigmaTime * sqrt(2) );
//    arg_0 = ( pow ( sPlotData.sigmaTime, 2 ) * ( OMEGA0(n) - omega ) + iii * sPlotData.tCenter ) /
//            ( sPlotData.sigmaTime * sqrt(2) );
//    preFactor = iii * sPlotData.sigmaTime * sqrt ( Pi / 2 ) *
//                exp ( - pow ( sPlotData.sigmaTime * ( OMEGA0(n) - omega ), 2 ) / 2. ) *
//                exp ( -iii * sPlotData.tCenter * ( OMEGA0(n) - omega ) ) *
//                exp ( -iii * omega * time );

    arg_t = ( pow ( sPlotData.sigmaTime, 2 ) * ( OMEGA0(n) - omega ) - iii * ( time - sPlotData.tCenter ) ) /
            ( sPlotData.sigmaTime * sqrt(2) );
    arg_0 = ( pow ( sPlotData.sigmaTime, 2 ) * ( OMEGA0(n) - omega ) + iii * sPlotData.tCenter ) /
            ( sPlotData.sigmaTime * sqrt(2) );
    preFactor_t = iii * sPlotData.sigmaTime * sqrt ( 2 ) *
                  exp ( - pow ( ( time - sPlotData.tCenter ) / sPlotData.sigmaTime, 2 ) / 2. ) *
                  exp ( -iii * OMEGA0(n) * time );
    preFactor_0 = iii * sPlotData.sigmaTime * sqrt ( 2 ) *
                  exp ( - pow ( ( sPlotData.tCenter ) / sPlotData.sigmaTime, 2 ) / 2. ) *
                  exp ( -iii * omega * time );

    result = preFactor_t * Faddeeva::Dawson ( arg_t ) - preFactor_0 * Faddeeva::Dawson ( arg_0 );
//
//    if ( isnan (result.real()) || isnan (result.imag()) )
//    {
//        cout << arg_t << endl << arg_0 << endl << preFactor_t << endl << preFactor_0 << endl;
//        exit ( EXIT_FAILURE );
//    }

    return result;

}

//Functions to construct real space magnetization vectors from Bloch's Theorem
//gsl_complex spatial_mx ( double x, double y, int mode )
//{
//    gsl_complex sum = gsl_complex_rect (0.0, 0.0 ), z;
//
//    for ( int i = 0; i < numvectors; i++ )
//    {
//        z = gsl_complex_mul ( gsl_matrix_complex_get ( evecs, i, mode ), gsl_complex_exp ( gsl_complex_rect ( 0.0, Gvx[i] * x + Gvy[i] * y ) ) );
//        sum = gsl_complex_add ( sum, z );
//    }
//    //kdotx = kx(kxl, Crystal) * x + ky(kyl, Crystal) * y;
//    //sum = gsl_complex_mul ( gsl_complex_exp ( gsl_complex_rect ( 0.0, kdotx ) ), sum );
//    return sum;
//}
//
//gsl_complex spatial_my ( double x, double y, int mode )
//{
//    gsl_complex sum = gsl_complex_rect (0.0, 0.0 ), z;
//
//    for ( int i = 0; i < numvectors; i++ )
//    {
//        z = gsl_complex_mul ( gsl_matrix_complex_get ( evecs, i + numvectors, mode ), gsl_complex_exp ( gsl_complex_rect ( 0.0, Gvx[i] * x + Gvy[i] * y ) ) );
//        sum = gsl_complex_add ( sum, z );
//    }
//    //kdotx = kx(kxl, Crystal) * x + ky(kyl, Crystal) * y;
//    //sum = gsl_complex_mul ( gsl_complex_exp ( gsl_complex_rect ( 0.0, kdotx ) ), sum );
//    return sum;
//}
//
//void SymmCalc ( double kxl, double kyl, FILE *evalsfile, FILE *evecsfile, Structure Crystal, int ExchField, int DipolarField )
//{
//    double kxl2, kyl2;
//
//    kxl2 = kxl; kyl2 = kyl;
//    CreateMmatrix ( kxl2, kyl2, Crystal, ExchField, DipolarField );
//    gsl_blas_dgemm ( CblasNoTrans, CblasNoTrans, 1.0, alphaMatrixInverse, Mmatrix, 0.0, MatProd );
//    gsl_eigen_nonsymmv ( MatProd, evals, evecs, w );
//    gsl_eigen_nonsymmv_sort ( evals, evecs, GSL_EIGEN_SORT_ABS_ASC );
//    gsl_vector_complex_fwrite ( evalsfile, evals );
//    gsl_matrix_complex_fwrite ( evecsfile, evecs );
//    cout << kxl2 << "  " << kyl2 << endl;
//}

//calculates magnetization vector in unit cell and examines symmetry (currently set only for square lattice)
void MagVecSymm ( LandauLifshitzGilbert cLLG, int ExchField, int DipolarField, int ElectricField )
{
    double kxl = 1e-50, kyl = 1e-50, kzl = 0.0, mode = 0;
    double xpos = 0, ypos = 0, xpos2 = 0, ypos2 = 0, spatial_stepsize = cLLG.Get_a() / 100, phase = 0;
    std::complex<double> mx, my, mx0, my0;
    ofstream MagVecSymm ( "Data/MagVecSymm.csv", ios::out | ios::trunc );

    MagVecSymm << "Col. 1:  xpos / a" << endl <<
                  "Col. 2:  ypos / a" << endl <<
                  "Col. 3:  ( |m| - |m_r| ) / |m|" << endl <<
                  "Col. 4:  Re ( ( mx - mx_r ) / mx )" << endl <<
                  "Col. 5:  Im ( ( mx - mx_r ) / mx )" << endl <<
                  "Col. 6:  Re ( ( my - my_r ) / my )" << endl <<
                  "Col. 7:  Im ( ( my - my_r ) / my )" << endl <<
                  "Col. 8:  Re (mx)" << endl <<
                  "Col. 9:  Im (mx)" << endl <<
                  "Col. 10: Re (my)" << endl <<
                  "Col. 11: Im (my)" << endl << endl;
    cLLG.ExportData ( MagVecSymm );

    cLLG.Solve ( kxl, kyl, kzl, ExchField, DipolarField, ElectricField );

    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
    phase = atan2 ( -mx0.imag(), mx0.real() );  //finds phase such that Im (mx) = 0
    exp ( iii * phase ) * cLLG.m_evecs;


    for ( xpos = -cLLG.Get_a() / 2; xpos < 0; xpos += spatial_stepsize )
        {
            for ( ypos = -cLLG.Get_a() / 2; ypos < xpos + spatial_stepsize / 2; ypos += spatial_stepsize )
            {
                cout << xpos << " " << ypos << endl;

                xpos2 = xpos; ypos2 = ypos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    MagVecSymm << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " ";
                    MagVecSymm << mx0.real() << " " << mx0.imag() << " " << my0.real() << " " << my0.imag() << " ";
                    //MagVecSymm << sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " " << sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << endl;
                }
                MagVecSymm << endl;

                xpos2 = xpos; ypos2 = -ypos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos, ypos, mode );
                    mx  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    if ( abs ( mx0 ) < 1e-10 )
                    {
                        mx = mx0;
                        my = my0;
                    }
    //                MagVecSymm << abs ( sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) - sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) ) / sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " "
    //                           << abs ( sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) - sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) ) / sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << " ";
    //                MagVecSymm << sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) << " " << sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) << endl;
                    MagVecSymm << abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                        sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                        sqrt ( norm ( mx0 ) + norm ( my0 ) ) << " ";
                    MagVecSymm << ( mx0.real() - mx.real() ) / ( mx0.real() ) << " " << ( mx0.imag() - mx.imag() ) / ( mx0.imag() ) << " "
                               << ( my0.real() - my.real() ) / ( my0.real() ) << " " << ( my0.imag() - my.imag() ) / ( my0.imag() ) << " ";
                    MagVecSymm << mx.real() << " " << mx.imag() << " " << my.real() << " " << my.imag() << " ";
                }
                MagVecSymm << endl;

                xpos2 = -xpos; ypos2 = ypos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos, ypos, mode );
                    mx  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    if ( abs ( mx0 ) < 1e-10 )
                    {
                        mx = mx0;
                        my = my0;
                    }
    //                MagVecSymm << abs ( sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) - sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) ) / sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " "
    //                           << abs ( sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) - sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) ) / sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << " ";
    //                MagVecSymm << sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) << " " << sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) << endl;
                    MagVecSymm << abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                        sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                        sqrt ( norm ( mx0 ) + norm ( my0 ) ) << " ";
                    MagVecSymm << ( mx0.real() - mx.real() ) / ( mx0.real() ) << " " << ( mx0.imag() - mx.imag() ) / ( mx0.imag() ) << " "
                               << ( my0.real() - my.real() ) / ( my0.real() ) << " " << ( my0.imag() - my.imag() ) / ( my0.imag() ) << " ";
                    MagVecSymm << mx.real() << " " << mx.imag() << " " << my.real() << " " << my.imag() << " ";
                }
                MagVecSymm << endl;

                xpos2 = -xpos; ypos2 = -ypos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos, ypos, mode );
                    mx  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    if ( abs ( mx0 ) < 1e-10 )
                    {
                        mx = mx0;
                        my = my0;
                    }
    //                MagVecSymm << abs ( sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) - sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) ) / sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " "
    //                           << abs ( sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) - sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) ) / sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << " ";
    //                MagVecSymm << sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) << " " << sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) << endl;
                    MagVecSymm << abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                        sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                        sqrt ( norm ( mx0 ) + norm ( my0 ) ) << " ";
                    MagVecSymm << ( mx0.real() - mx.real() ) / ( mx0.real() ) << " " << ( mx0.imag() - mx.imag() ) / ( mx0.imag() ) << " "
                               << ( my0.real() - my.real() ) / ( my0.real() ) << " " << ( my0.imag() - my.imag() ) / ( my0.imag() ) << " ";
                    MagVecSymm << mx.real() << " " << mx.imag() << " " << my.real() << " " << my.imag() << " ";
                }
                MagVecSymm << endl;

                xpos2 = ypos; ypos2 = xpos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos, ypos, mode );
                    mx  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    if ( abs ( mx0 ) < 1e-10 )
                    {
                        mx = mx0;
                        my = my0;
                    }
    //                MagVecSymm << abs ( sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) - sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) ) / sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " "
    //                           << abs ( sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) - sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) ) / sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << " ";
    //                MagVecSymm << sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) << " " << sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) << endl;
                    MagVecSymm << abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                        sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                        sqrt ( norm ( mx0 ) + norm ( my0 ) ) << " ";
                    MagVecSymm << ( mx0.real() - mx.real() ) / ( mx0.real() ) << " " << ( mx0.imag() - mx.imag() ) / ( mx0.imag() ) << " "
                               << ( my0.real() - my.real() ) / ( my0.real() ) << " " << ( my0.imag() - my.imag() ) / ( my0.imag() ) << " ";
                    MagVecSymm << mx.real() << " " << mx.imag() << " " << my.real() << " " << my.imag() << " ";
                }
                MagVecSymm << endl;

                xpos2 = ypos; ypos2 = -xpos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos, ypos, mode );
                    mx  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    if ( abs ( mx0 ) < 1e-10 )
                    {
                        mx = mx0;
                        my = my0;
                    }
        //                MagVecSymm << abs ( sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) - sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) ) / sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " "
        //                           << abs ( sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) - sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) ) / sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << " ";
        //                MagVecSymm << sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) << " " << sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) << endl;
                    MagVecSymm << abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                        sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                        sqrt ( norm ( mx0 ) + norm ( my0 ) ) << " ";
                    MagVecSymm << ( mx0.real() - mx.real() ) / ( mx0.real() ) << " " << ( mx0.imag() - mx.imag() ) / ( mx0.imag() ) << " "
                               << ( my0.real() - my.real() ) / ( my0.real() ) << " " << ( my0.imag() - my.imag() ) / ( my0.imag() ) << " ";
                    MagVecSymm << mx.real() << " " << mx.imag() << " " << my.real() << " " << my.imag() << " ";
                }
                MagVecSymm << endl;

                xpos2 = -ypos; ypos2 = xpos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos, ypos, mode );
                    mx  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    if ( abs ( mx0 ) < 1e-10 )
                    {
                        mx = mx0;
                        my = my0;
                    }
        //                MagVecSymm << abs ( sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) - sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) ) / sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " "
        //                           << abs ( sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) - sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) ) / sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << " ";
        //                MagVecSymm << sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) << " " << sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) << endl;
                    MagVecSymm << abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                        sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                        sqrt ( norm ( mx0 ) + norm ( my0 ) ) << " ";
                    MagVecSymm << ( mx0.real() - mx.real() ) / ( mx0.real() ) << " " << ( mx0.imag() - mx.imag() ) / ( mx0.imag() ) << " "
                               << ( my0.real() - my.real() ) / ( my0.real() ) << " " << ( my0.imag() - my.imag() ) / ( my0.imag() ) << " ";
                    MagVecSymm << mx.real() << " " << mx.imag() << " " << my.real() << " " << my.imag() << " ";
                }
                MagVecSymm << endl;

                xpos2 = -ypos; ypos2 = -xpos;
                MagVecSymm << xpos2 / cLLG.Get_a() << " " << ypos2 / cLLG.Get_a() << " ";
                for ( mode = 0; mode < 20; mode++ )
                {
                    mx0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos, ypos, mode );
                    my0 = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos, ypos, mode );
                    mx  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_x, xpos2, ypos2, mode );
                    my  = cLLG.ComputeFields_kdotx ( MAGNETIZATION_y, xpos2, ypos2, mode );
                    if ( abs ( mx0 ) < 1e-10 )
                    {
                        mx = mx0;
                        my = my0;
                    }
        //                MagVecSymm << abs ( sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) - sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) ) / sqrt ( pow ( mx0.real(), 2 ) + pow ( mx0.imag(), 2 ) ) << " "
        //                           << abs ( sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) - sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) ) / sqrt ( pow ( my0.real(), 2 ) + pow ( my0.imag(), 2 ) ) << " ";
        //                MagVecSymm << sqrt ( pow ( mx.real(), 2 ) + pow ( mx.imag(), 2 ) ) << " " << sqrt ( pow ( my.real(), 2 ) + pow ( my.imag(), 2 ) ) << endl;
                    MagVecSymm << abs ( sqrt ( norm ( mx0 ) + norm ( my0 ) ) -
                                        sqrt ( norm ( mx  ) + norm ( my  ) ) ) /
                                        sqrt ( norm ( mx0 ) + norm ( my0 ) ) << " ";
                    MagVecSymm << ( mx0.real() - mx.real() ) / ( mx0.real() ) << " " << ( mx0.imag() - mx.imag() ) / ( mx0.imag() ) << " "
                               << ( my0.real() - my.real() ) / ( my0.real() ) << " " << ( my0.imag() - my.imag() ) / ( my0.imag() ) << " ";
                    MagVecSymm << mx.real() << " " << mx.imag() << " " << my.real() << " " << my.imag() << " ";
                }
                MagVecSymm << endl;
            }
        }
}

//These funcions obtained from https://forum.kde.org/viewtopic.php?f=74&t=107565
//template<typename MatrixType>
//void WriteMatrix(string filename, const MatrixType& m)
//{
//    ofstream f(filename.c_str(), ios::binary);
//    f.write((char *)&m.rows(), sizeof(m.rows()));
//    f.write((char *)&m.cols(), sizeof(m.cols()));
//    f.write((char *)&m.data(), sizeof(typename MatrixType::Scalar)*m.cols()*m.cols());
//    f.close();
//}
//
//template<typename MatrixType>
//void ReadMatrix(string filename, MatrixType& m)
//{
//    typename MatrixType::Index rows, cols;
//    ifstream f(filename.c_str(), ios::binary);
//    f.read((char *)&rows, sizeof(rows));
//    f.read((char *)&cols, sizeof(cols));
//    m.resize(rows, cols);
//    f.read((char *)&m.data(), sizeof(typename MatrixType::Scalar)*rows*cols);
//    if (f.bad())
//    throw std::exception("Error reading matrix");
//    f.close();
//}

void WriteMatrix(string filename, MatrixXcd& m, int rows, int cols)
{
    ofstream f(filename.c_str(), ios::binary | ios::trunc);
    f.write((char *)m.data(), 16 * rows * cols);
    f.close();
}

void ReadMatrix(string filename, MatrixXcd& m, int rows, int cols)
{
    ifstream f(filename.c_str(), ios::binary);
    f.read((char *)m.data(), 16 * rows * cols);
//    if (f.bad())
//        throw std::exception("Error reading matrix");
    f.close();
}

void WriteMatrixReal (string filename, MatrixXd& m, int rows, int cols)
{
    ofstream f(filename.c_str(), ios::binary | ios::trunc);
    f.write((char *)m.data(), 8 * rows * cols);
    f.close();
}

void ReadMatrixReal (string filename, MatrixXd& m, int rows, int cols)
{
    ifstream f(filename.c_str(), ios::binary);
    f.read((char *)m.data(), 8 * rows * cols);
//    if (f.bad())
//        throw std::exception("Error reading matrix");
    f.close();
}

