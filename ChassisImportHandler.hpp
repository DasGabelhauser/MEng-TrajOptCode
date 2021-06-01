#ifndef CHASSIS_IMPORT_HANDLER
#define CHASSIS_IMPORT_HANDLER



#include "MeshOps.hpp"

namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    // structure storing the data of one line of an ARE project blueprint file
    typedef struct {
        int skeleton;
        int organType;
        double posX;
        double posY;
        double posZ;
        double rotX;
        double rotY;
        double rotZ;
    } blueprintLine;

    /*  Read in the data from a blueprint file.
        args:
        filePathName - the absolute filepath and filename of the blueprint file.
        returns:
        an std::vector of blueprintLine objects containing each line of the blueprint
        file in order
    */
    std::vector<blueprintLine> readChassisBlueprint(std::string filePathName);

    /*  Import a chassis file and convex decompose it into convex sub-hulls.
        args: 
        filePathName - the absolute filepath and filename of the chassis file. The chassis file should not have the organ clips added.
        blueprint - an std::vector of blueprintLine objects desribing where to add organ clips to the model
        returns:
        an std::vector of Hull objects that collectively represent the chassis model
    */
    std::vector<Hull> importConvexDecomposedChassis(std::string filePathName, std::vector<blueprintLine> blueprint);

    /*  Import a chassis file and convex decompose it into convex sub-hulls.
        args: 
        filePathName - the absolute filepath and filename of the chassis file. The chassis file should have the organ clips added.
        returns:
        a Hull object representing the convex hull of the chassis model.
    */
    Hull importConvexHullChassis(std::string filePathName);
}

#endif //CHASSIS_IMPORT_HANDLER