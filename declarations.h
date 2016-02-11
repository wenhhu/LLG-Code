#ifndef DECLARATIONS_H_INCLUDED
#define DECLARATIONS_H_INCLUDED

#include <vector>
#include <string>
//#include <gsl/gsl_matrix.h>

#include <Eigen/Dense>

struct Material
{
    double A;
    double Ms;
    double alpha;
    double Q;
    double ExchLength;
    double DMvector;
};

enum Shape
{
    SHAPE_CIRCLE,
    SHAPE_ELLIPSE,
};

enum Interface
{
    INTERFACE_IDEAL,
    INTERFACE_DIFFUSE,
};

enum SlabType
{
    FINITE,
    INFINITE,
};

struct Structure
{
    double Pi;
    double mu0;
    double gyroRatio;
    double a;
    double f;
    double t_Interface;
    double H0;
    double R;
    double Ratio;  //ratio of x-axis to y-axis for ellipse
    double SlabWidth;
    double sigmaSpace;
    Shape eShape;
    Interface eInterface;
    SlabType eSlabType;
    double eCharge;
    double eField;
    double E_so;
    double ExchangeTerm;
    std::complex <double> DMTerm;
};

struct Precision
{
    int    Nmax;
    int    NmaxLimiter;
    int    Nmax_Sum;
    int    N_Interface;
    int    nummodes;
    int    spatial_sizePrime;
    int    numPositionsPrime;
    double spatial_stepsizePrime;
    double stepsize;
    double DeltaOmega;
    double MaxOmega;
};

struct PlotData
{
    int spatial_size;
    int numPositions;
    int temporalSteps;
    double xmax;
    double ymax;
    double temporal_stepsize;
    double EndTime;
    int NumFreq;
    Eigen::VectorXd OMEGA0_nodim;
    Eigen::VectorXcd OMEGA0;
    double sigmaTime;
    double tCenter;


//    int NumFreq = 4;
//    double OMEGA0_nodim [ ] = { 5, 15, 20, 50 };
//    gsl_vector_complex *OMEGA0 = gsl_vector_complex_calloc ( NumFreq );
};

enum Fields
{
    SATURATION_MAGNETIZATION,
    EXCHANGE_LENGTH,
    GILBERT_DAMPING,
    MAGNETOSTATIC_POTENTIAL,
    EXCHANGE_x,
    EXCHANGE_y,
    EXCHANGE_x_ALT,
    EXCHANGE_y_ALT,
    DIPOLAR_x,
    DIPOLAR_y,
    CURL_OF_DIPOLAR_xWRTy,
    CURL_OF_DIPOLAR_yWRTx,
    MAGNETIZATION_x,
    MAGNETIZATION_y,
    CURL_OF_MAGNETIZATION_xWRTy,
    CURL_OF_MAGNETIZATION_yWRTx,
    DEMAG_FIELD,
};

//class definition headers
#include "LLG.h"
#include "FieldPlots.h"
#include "Response.h"
#include "DensityOfStates.h"


std::string NameFILE ( std::string FileName, int Nmax, double stepsize, int PartNum, std::string Extension );
double rotatexcomp ( double xcomp, double ycomp, double angle );
double rotateycomp ( double xcomp, double ycomp, double angle );
double norm ( double xcomp, double ycomp );
void Dispersions ( LandauLifshitzGilbert cLLG, int ExchField, int DipolarField, int ElectricField, int PartNum );
double kx ( double length, Structure Crystal );
double ky ( double length, Structure Crystal );
double kz ( double length, Structure Crystal );
double DistributionFunction ( double xprime, double yprime, Structure Crystal );
std::complex <double> FourierIntegral ( std::complex <double> omega, double time, int n, double EndTime, Eigen::VectorXcd& OMEGA0, PlotData sPlotData, double Pi );
//gsl_complex spatial_mx ( double x, double y, int mode );   //Uses global variable evecs
//gsl_complex spatial_my ( double x, double y, int mode );   //Uses global variable evecs
void SymmCalc ( double kxl, double kyl, FILE *evalsfile, FILE *evecsfile, Structure Crystal, int ExchField, int DipolarField, int ElectricField );
void MagVecSymm ( LandauLifshitzGilbert cLLG, int ExchField, int DipolarField, int ElectricField );
void WriteMatrix ( std::string filename, Eigen::MatrixXcd& m, int rows, int cols);  //writes/reads Gxx/Gxy/Gyx/Gyy matrices to/from a binary file
void ReadMatrix  ( std::string filename, Eigen::MatrixXcd& m, int rows, int cols);
void WriteMatrixReal ( std::string filename, Eigen::MatrixXd& m, int rows, int cols);
void ReadMatrixReal  ( std::string filename, Eigen::MatrixXd& m, int rows, int cols);

#endif // DECLARATIONS_H_INCLUDED
