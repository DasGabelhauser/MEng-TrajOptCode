#include "DenavitHartenberg.hpp"
#include <cmath>

namespace EricsTrajOpt {
    using namespace std;

    /*  Initialise a DH matrix with DH parameters d, theta, r and alpha.
        d and r are in arbitrary distance units.
        theta and alpha are in degrees. */
    DenavitHartenberg::DenavitHartenberg(const double d, const double theta, const double r, const double alpha) {
        //convert angle units
        double radTheta = (theta * M_PI) / 180;
        double radAlpha = (alpha * M_PI) / 180;
        //calculate DH matrix
        double rotMat[4][4] = { {cos(radTheta), -sin(radTheta)*cos(radAlpha), sin(radTheta)*sin(radAlpha) , r*cos(radTheta)},
                                {sin(radTheta), cos(radTheta)*cos(radAlpha) , -cos(radTheta)*sin(radAlpha), r*sin(radTheta)},
                                {0            , sin(radAlpha)               , cos(radAlpha)               , d              },
                                {0            , 0                           , 0                           , 1              }};
        //copy to DHMat member. 
        for (int row = 0; row < 4; row++) {
            for (int col = 0; col < 4; col++) {
                m_DHMat[row][col] = rotMat[row][col];
            }
        }
    }
    /* overloaded * operator. multiplies together two DH matrices */
    DenavitHartenberg DenavitHartenberg::operator* (DenavitHartenberg b) {
        DenavitHartenberg returnDH = *new DenavitHartenberg();
        double val;
        for (int row = 0; row < 4; row++) {
            for (int col = 0; col < 4; col++) {
                //accumulate value of current element of new DH object
                val = 0;
                for (int i = 0; i < 4; i++) {
                    val += m_DHMat[row][i]*b.getDHMat()[i][col];
                }
                //set value of current element of new DH object
                returnDH.setElement(val, col, row);
            }
        }
        return returnDH;
    }

    /*  set an element of the DH matrix to a value.
        args:
        val - the value to set the element to.
        col - the column index of the element.
        row - the row index of the element.
        returns:
        false if row or column indices are out of range for the DH matrix.
        true if row and column indices are valid.
    */
    bool DenavitHartenberg::setElement (const double val, const int col, const int row) {
        //check indices are within accpetable range
        if (col < 0 or col > 3 or row < 0 or row > 3) {
            return false;
        }
        //set value
        m_DHMat[row][col] = val;
        return true;
    }


}
