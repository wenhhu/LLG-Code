#ifndef LLG_H_INCLUDED
#define LLG_H_INCLUDED

#include <vector>
#include <fstream>

//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_eigen.h>

#include <Eigen/Dense>

#include "declarations.h"

class LandauLifshitzGilbert
{
    protected:
        //protected variables
        Material  m_Mat_A;
        Material  m_Mat_B;
        Structure m_Crystal;
        Precision m_sPrecision;
        int       m_numvectors, m_numvectors_Sum;
        double    m_kxl;
        double    m_kyl;
        double    m_kzl;
        std::vector <double> m_Gvx, m_Gvy, m_Gvx_Sum, m_Gvy_Sum;
        std::vector < std::vector <double> > m_Id, m_Qarray, m_Msarray, m_Msarray_Sum, m_DMarray_Sum,
                                             m_ExchLengthArray, m_ExchLengthArray_Sum;
//        gsl_matrix                   *m_Mmatrix,
//                                     *m_alphaMatrix,
//                                     *m_alphaMatrixInverse,
//                                     *m_MatProd;
//        gsl_eigen_nonsymmv_workspace *m_w;

        Eigen::MatrixXd m_Mmatrix,
                        m_alphaMatrix,
                        m_alphaMatrixInverse,
                        m_MatProd;

        Eigen::MatrixXcd m_MmatrixComplex,
                         m_alphaMatrixInverseComplex;

        Eigen::Matrix<int, 6, 2> m_Hex_Gv_indexer;
        Eigen::VectorXi m_Gv1_index, m_Gv2_index;

        //protected default constructor
        LandauLifshitzGilbert () {}

        //protected functions
        void CreateGvectors ();  //defines elements of Gvx, Gvy: dependent on
        double Structure_Factor ( int l, int j, int n );
        double C_n ( double C_a, double C_b, int n );
        double LanczosSigmaFactor ( int l, int j );
        double FT_Variable ( double C_a, double C_b, int l, int j );
        double Structure_Factor_Sum ( int l, int j, int n );
        double C_n_Sum ( double C_a, double C_b, int n );
        double LanczosSigmaFactor_Sum ( int l, int j );
        double FT_Variable_Sum ( double C_a, double C_b, int l, int j );
        void BuildMatrices ();  //Builds all matrices independent of kx, ky
        double ComputeExchange ( int i, int j ); //obsolete - ComputeSumTerms is used in place of this function
        std::complex <double> ComputeDMterm   ( int i, int j, int ExchField ); //obsolete - ComputeSumTerms is used in place of this function
        void ComputeSumTerms ( int i, int j, int ExchField, int ElectricField );  //computes both exchange and DM terms at once to save computation time
        double ComputeSlabDemagField ( int i, int j );
        void CreateMmatrix ( int ExchField, int DipolarField, int ElectricField );
    public:
        //constructor
        LandauLifshitzGilbert ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision );

        //public variables
        Eigen::VectorXcd m_evals;
        Eigen::MatrixXcd m_evecs;
        std::vector <int> m_evalsSortedIndex;

        //public functions
        double Get_a () { return m_Crystal.a; }
        Precision Get_Precision () { return m_sPrecision; }
        Structure Get_Crystal () { return m_Crystal; }
        void SetParameters ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision );  //sets MC parameters and initializes matrices independent of k
        void ExportData ( std::ofstream& filename ); //exports variables to beginning of chosen data file. data will then start at index 1 for pyxplot
        void BubbleSort ( int ElectricField );
        void Solve           ( double kxl, double kyl, double kzl, int ExchField, int DipolarField, int ElectricField );  //solves LLG equation with solutions in evals and evecs
        void Solve_NoVectors ( double kxl, double kyl, double kzl, int ExchField, int DipolarField, int ElectricField );  //solves LLG equation with solutions in evals
        std::complex <double> ComputeFields       ( Fields eFieldName, double x, double y, int mode );
        std::complex <double> ComputeFields_kdotx ( Fields eFieldName, double x, double y, int mode );  //fields multiplied by e^(ikr)
};


#endif // LLG_H_INCLUDED
