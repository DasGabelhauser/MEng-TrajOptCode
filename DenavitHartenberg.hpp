#ifndef DENAVIT_HARTENBERG
#define DENAVIT_HARTENBERG

#include <cmath>
#include <vector>

namespace EricsTrajOpt {
    using namespace EricsTrajOpt;
    /* class to represent a Denavit-Hartenberg matrix. */
    class DenavitHartenberg{
        public:
        //initialise a DH class with all elements of rotation matrix as 0
        DenavitHartenberg() {}
        
        /*  Initialise a DH matrix with DH parameters d, theta, r and alpha.
            d and r are in arbitrary distance units.
            theta and alpha are in degrees. */
        DenavitHartenberg(const double d, const double theta, const double r, const double alpha);

        /* overloaded * operator. multiplies together two DH matrices */
        DenavitHartenberg operator* (DenavitHartenberg b);

        /*  set an element of the DH matrix to a value.
            args:
            val - the value to set the element to.
            col - the column index of the element.
            row - the row index of the element.
            returns:
            false if row or column indices are out of range for the DH matrix.
            true if row and column indices are valid.
        */
        bool setElement (const double val, const int col, const int row);

        std::vector<std::vector<double>> getDHMat() {
            return m_DHMat;
        }
        
        /*  Getter for the rotation sub-matrix of the DH matrix*/
        std::vector<std::vector<double>> getRotMat() {
            std::vector<std::vector<double>> rotMat = {
                {m_DHMat[0][0], m_DHMat[0][1], m_DHMat[0][2]}, 
                {m_DHMat[1][0], m_DHMat[1][1], m_DHMat[1][2]}, 
                {m_DHMat[2][0], m_DHMat[2][1], m_DHMat[2][2]}};
            return rotMat;
        }

        /*  Getter for the translation sub-vector of the DH matrix*/
        std::vector<double> getTransVec () {
            std::vector<double> transVec = {m_DHMat[0][3], m_DHMat[1][3], m_DHMat[2][3]};
            return transVec;
        }

        private:

        /* 2D vector representing a DH matrix */
        std::vector<std::vector<double>> m_DHMat = 
            {  
                {0,0,0,0},
                {0,0,0,0},
                {0,0,0,0},
                {0,0,0,0}   
            };
    };
}

#endif //DENAVIT_HARTENBERG