#ifndef DENSITYOFSTATES_H_INCLUDED
#define DENSITYOFSTATES_H_INCLUDED

//calculates dispersion contours in BZ and obtains density of states
class DensityOfStates : public LandauLifshitzGilbert
{
    private:
        //protected default constructor
        DensityOfStates () {}

        //private variables
        time_t m_start, m_finish;
        double m_timeEst;

    public:
        //constructor
        DensityOfStates ( Material Mat_A, Material Mat_B, Structure Crystal, Precision sPrecision );

        //public functions
        void CalculateDoS ( int PartNum, double kxlstart, double kylstart, double kxlend, double kylend,
                            double kMin, double kMax, int ExchField, int DipolarField, int ElectricField, int kTotal );  //calculates dispersion contours and denstity of states
        void AddDoSMatrix ();  //adds DensityOfStates files from calculations
};

#endif // DENSITYOFSTATES_H_INCLUDED
