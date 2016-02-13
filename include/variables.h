#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED

//constants and crystal parameters
Structure Crystal;
Crystal.Pi          = 3.14159265358979;
Crystal.mu0         = 4e-7 * Crystal.Pi;
Crystal.gyroRatio   = 1.76084e11;
Crystal.a           = 400e-9;
Crystal.f           = 0.4;
Crystal.t_Interface = 2.0e-9;
Crystal.H0          = 7.9577e4 * 0.1;
Crystal.Ratio       = 1.5;  //set to >1 for major axis along x-axis and < 1 for major axis along y-axis
Crystal.SlabWidth   = 500e-9 / 2; //Note: actual thickness is two times this value
Crystal.eCharge     = 1.60217657e-19;
Crystal.eField      = 1.5e7;
Crystal.E_so        = 6.10427e-19;

#ifdef IDEAL
Crystal.t_Interface = 1.0e-50;
Crystal.eInterface = INTERFACE_IDEAL;
#endif

#ifdef DIFFUSE
Crystal.eInterface = INTERFACE_DIFFUSE;
#endif

#ifdef SLAB
Crystal.eSlabType = FINITE;
#endif

#ifndef SLAB
Crystal.eSlabType = INFINITE;
#endif

#ifdef CIRCLE
Crystal.eShape = SHAPE_CIRCLE;
Crystal.Ratio  = 1.0;
if ( lattice == 1 )      Crystal.R = sqrt ( Crystal.Ratio * (Crystal.f * Crystal.a * Crystal.a) / Crystal.Pi ) - 0.5 * Crystal.t_Interface;
else if ( lattice == 2 ) Crystal.R = Crystal.a * sqrt ( Crystal.Ratio * ( Crystal.f / Crystal.Pi ) * 3 / ( 2 * sqrt (3) ) ) - 0.5 * Crystal.t_Interface;
Crystal.sigmaSpace = Crystal.R / 4;
#endif

#ifdef ELLIPSE
Crystal.eShape = SHAPE_ELLIPSE;
if ( lattice == 1 )      Crystal.R = sqrt ( Crystal.Ratio * (Crystal.f * Crystal.a * Crystal.a) / Crystal.Pi ) - 0.5 * Crystal.t_Interface;
else if ( lattice == 2 ) Crystal.R = Crystal.a * sqrt ( Crystal.Ratio * ( Crystal.f / Crystal.Pi ) * 3 / ( 2 * sqrt (3) ) ) - 0.5 * Crystal.t_Interface;
Crystal.sigmaSpace = Crystal.R / 4;
#endif

//magnetic properties for various materials
double A = 1.38, B = 0.15, C = 0.1826, phi = 1 - A * C - B * C * C;   // constants for FeGa

//               {      A,           Ms,   alpha,     a   }
Material Fe    = { 8.3e-12,    1.71092e6,   0.0019 },
         Co    = { 10.3e-12,   1.40056e6,   0.011  },
         Ni    = { 3.4e-12,    0.485423e6,  0.064  },
         SmCo  = { 22.0e-12,   8.514789e5,  0.0    },
         BaFeO = { 4.1e-12,    3.740141e5,  0.0    },
         YIG   = { 4.15e-12,   1.5915494e5, 0.0006 }, //check alpha for YIG
         FeGa  = { Fe.A * phi, 1.38426e6,   0.0    }, //for ~18% Ga
         Vacuum = { 4.15e-12 * 0.9,  1.5915494e5 * 0.9, 0.0006 },
         Mat_A = YIG,  //cylinders
         Mat_B = Vacuum;  //host
//         Mat_A.alpha = 0;
//         Mat_B.alpha = 0;

if ( Mat_A.Ms < 1e-8 )
{
    Mat_A.Q = 0;
    Mat_A.ExchLength = 0;
}
else
{
    Mat_A.Q = ( 2 * Mat_A.A ) / ( Mat_A.Ms * Crystal.mu0 * Crystal.H0 );
    Mat_A.ExchLength = 2 * Mat_A.A / ( Crystal.mu0 * Mat_A.Ms * Mat_A.Ms );
}

if ( Mat_B.Ms < 1e-8 )
{
    Mat_B.Q = 0;
    Mat_B.ExchLength = 0;
}
else
{
    Mat_B.Q = ( 2 * Mat_B.A ) / ( Mat_B.Ms * Crystal.mu0 * Crystal.H0 );
    Mat_B.ExchLength = 2 * Mat_B.A / ( Crystal.mu0 * Mat_B.Ms * Mat_B.Ms );
}

//Parameters that only affect plotting region, detail, etc.
PlotData sPlotData;
sPlotData.spatial_size      = 100;
sPlotData.numPositions      = ( 2 * sPlotData.spatial_size + 1 ) * ( 2 * sPlotData.spatial_size + 1 );
sPlotData.xmax              = 20.0;
sPlotData.ymax              = sPlotData.xmax;
sPlotData.temporal_stepsize = 0.02e-9;
sPlotData.EndTime           = 0.01e-9;  //does nothing when using Gaussian distribution in FourierIntegral function
sPlotData.NumFreq           = 32;    //must match number of frequencies defined below
sPlotData.OMEGA0_nodim.resize ( sPlotData.NumFreq );
sPlotData.OMEGA0.resize       ( sPlotData.NumFreq );
sPlotData.sigmaTime         = 0.1e-9;
sPlotData.tCenter           = 4 * sPlotData.sigmaTime;

for ( int i = 0; i < sPlotData.NumFreq; i++ )
{
    sPlotData.OMEGA0 ( i ) = 8.5e9 + i * 0.0625e9;
}
//sPlotData.OMEGA0_nodim(0) = 5.0;
//sPlotData.OMEGA0_nodim(1) = 10.0;
//sPlotData.OMEGA0_nodim(2) = 15.0;
//sPlotData.OMEGA0_nodim(3) = 20.0;

#ifdef PULSE
sPlotData.temporalSteps = 20;
#endif

#ifdef GREEN_FUNCTION
sPlotData.temporalSteps = 1;
#endif

//int const spatial_size  = 50,
//          numPositions  = ( 2 * spatial_size + 1 ) * ( 2 * spatial_size + 1 );             //spatial_size_prime needs to be a multiple of spatial_size / ( 2 * (xmax || ymax) )
//#ifdef PULSE
//int const temporalSteps = 20;
//#endif
//
//#ifdef GREEN_FUNCTION
//int const temporalSteps = 1;
//#endif
//
//double    xmax              = 20.0,
//          ymax              = xmax,    //xmax, ymax are in units of a
//          temporal_stepsize = 0.01e-9,
//          EndTime           = 0.01e-9;
//
//int NumFreq = 4;
//double OMEGA0_nodim [ ] = { 5, 15, 20, 50 };
//gsl_vector_complex *OMEGA0 = gsl_vector_complex_calloc ( NumFreq );

//Variables affecting precision of calculations
Precision sPrecision;
sPrecision.Nmax                  = 6; //Nmax >= NmaxLimiter
sPrecision.NmaxLimiter           = 6; //Limits Nmax in Ms, alpha and exchange length
sPrecision.Nmax_Sum              = sPrecision.Nmax + sPrecision.NmaxLimiter;
sPrecision.N_Interface           = 20;
if ( lattice == 1 )      sPrecision.nummodes = 121;   //Square -> nummodes = 121(Nmax=5), Nmax=12; Hex -> nummodes = 127(Nmax=6), Nmax=14
else if ( lattice == 2 ) sPrecision.nummodes = 127;
sPrecision.spatial_sizePrime     = 50;
sPrecision.numPositionsPrime     = ( 2 * sPrecision.spatial_sizePrime + 1 ) * ( 2 * sPrecision.spatial_sizePrime + 1 );
sPrecision.spatial_stepsizePrime = 1. / sPrecision.spatial_sizePrime;
sPrecision.stepsize              = 0.1; //k-space step size
sPrecision.DeltaOmega            = 0.01e9; //freq. step size for calculating density of states
sPrecision.MaxOmega              = 50e9;  //max frequency used in DoS calculation

//Toggles and misc options
int const ExchField     = 1,    //Set to 1 for derived, 2 for alternate, 0 for no exchange field
          DipolarField  = 1,    //set to zero to artificially remove dipolar field
          ElectricField = 1,    //set to one to activate E-field, 0 to remove it
          processes    = 156;

#ifdef PULSE
int const FreqGF       = 0;    //set temporalSteps to one if FreqGF=1
#endif

#ifdef GREEN_FUNCTION
int const FreqGF       = 1;    //set temporalSteps to one if FreqGF=1
#endif

#endif // VARIABLES_H_INCLUDED
