#include "constants.h"
#include "declarations.h"

using namespace std;
using namespace Eigen;


DensityOfStates::DensityOfStates ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision )
    : LandauLifshitzGilbert ( Mat_A, Mat_B, Crystal, sPrecision )
    {

    }

void DensityOfStates::CalculateDoS ( int PartNum, double kxlstart, double kylstart, double kxlend, double kylend,
                                     double kMin, double kMax, int ExchField, int DipolarField, int ElectricField, int kTotal )  //calculates dispersion contours and denstity of states
{
    //Declare variables
    int mode, kCount = 0, OmegaIndex = 0;
    double kxl, kyl, kzl = 0, EigenOmega;
    MatrixXd DoSMatrix = MatrixXd::Zero ( m_sPrecision.MaxOmega / m_sPrecision.DeltaOmega, 2 );
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

    //initialize first column of DoSMatrix:
    for ( int i = 0; i < DoSMatrix.rows(); i++ )
    {
        DoSMatrix ( i, 0 ) = i * m_sPrecision.DeltaOmega;
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

            Solve ( kxl, kyl, kzl, ExchField, DipolarField, ElectricField );

//          Output k-space contour data:
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
                OmegaIndex = 0;
                EigenOmega = m_evals ( m_evalsSortedIndex.at(mode) ).imag() * -m_Crystal.gyroRatio * m_Crystal.mu0 * m_Crystal.H0;
                if ( EigenOmega < 0 ) continue;  //ignore negative frequencies

                while (1)
                {
                    if ( EigenOmega > m_sPrecision.MaxOmega ) break;
                    //if ( EigenOmega < ( OmegaIndex + 1 ) * m_sPrecision.DeltaOmega )
                    //{
                        DoSMatrix ( OmegaIndex, 1 ) += exp ( -pow ( DoSMatrix ( OmegaIndex, 0 ) + m_sPrecision.DeltaOmega / 2 - EigenOmega, 2 ) /
                                                             ( 2 * pow ( 5 * m_sPrecision.DeltaOmega, 2 ) ) ) /
                                                       ( sqrt ( 2 * m_Crystal.Pi ) * 5 * m_sPrecision.DeltaOmega );
                        //break;
                    //}
                    OmegaIndex += 1;
                    if ( OmegaIndex >= DoSMatrix.rows() ) break;
                }
            }
            time ( &m_finish );
            m_timeEst = difftime ( m_finish, m_start ) * ( kTotal * ( PartNum + 1 )  - kCount ) / ( kCount - kTotal * PartNum ) / 3600.0;
            if ( m_timeEst > 24 ) ProgressFile << kxl << " " << kyl << " ETA: " << m_timeEst / 24 << " days" << endl;
            else ProgressFile << endl << kxl << " " << kyl << " ETA: " << m_timeEst << " hours" << endl;
        }
    }
    if      ( lattice == 1 ) DoSMatrix.col (1) = DoSMatrix.col (1) * m_sPrecision.stepsize * m_sPrecision.stepsize * ( m_Crystal.Pi / m_Crystal.a ) * ( m_Crystal.Pi / m_Crystal.a );
    else if ( lattice == 2 ) DoSMatrix.col (1) = DoSMatrix.col (1) * m_sPrecision.stepsize * m_sPrecision.stepsize * ( 2 * m_Crystal.Pi / ( m_Crystal.a * sqrt(3) ) ) * ( 2 * m_Crystal.Pi / ( m_Crystal.a * sqrt(3) ) );
    ProgressFile.close(); ERRORfile.close(); DispersionFile.close(); LinewidthFile.close();

    WriteMatrixReal ( NameFILE ( "DoS/DensityOfStates", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".bin" ).c_str(), DoSMatrix, DoSMatrix.rows(), DoSMatrix.cols() );
    ofstream DoSFile ( "DoS/DensityOfStates.dat", ios::out | ios::trunc );
    DoSFile << DoSMatrix;
}

void DensityOfStates::AddDoSMatrix () //adds DensityOfStates files from calculations
{
    int PartNum = 0;
    MatrixXd DoSMatrixTemp = MatrixXd::Zero ( m_sPrecision.MaxOmega / m_sPrecision.DeltaOmega, 2 ),
             DoSMatrix     = MatrixXd::Zero ( m_sPrecision.MaxOmega / m_sPrecision.DeltaOmega, 2 );
    ofstream DoSFile ( "DoS/DensityOfStates.dat", ios::out | ios::trunc );

    while (1)
    {
        FILE *DoSFileCheck = fopen ( NameFILE ( "DoS/DensityOfStates", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".bin" ).c_str(), "r" );
        if ( DoSFileCheck == NULL ) break;

        ReadMatrixReal ( NameFILE ( "DoS/DensityOfStates", m_sPrecision.Nmax, m_sPrecision.stepsize, PartNum, ".bin" ).c_str(), DoSMatrixTemp, DoSMatrixTemp.rows(), DoSMatrixTemp.cols() );
        DoSMatrix.col(1) += DoSMatrixTemp.col(1);
        PartNum++;
    }
    DoSMatrix.col(0) = DoSMatrixTemp.col(0);

    DoSFile << DoSMatrix;
}
