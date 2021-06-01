#ifndef CHASSIS_IMPORT_HANDLER
#define CHASSIS_IMPORT_HANDLER

#include "coppeliaGeometricRoutines-master/geom.h"
#include <string>

namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

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

    /*  Read an STL file into a coppelia CObbStruct object.
        args:
        filePathName - the absolute filepath and filename of the STL file.
        returns:
        a CObbStruct* containing the mesh data.
    */
    CObbStruct* readSTL2CoppMesh(std::string filePathName);

    /*  Read multiple STL files into a single coppelia CObbStruct object.
        args:
        filepathnames - a std::vector of the absolute filepaths and filenames of the STL files.
        returns:
        a CObbStruct* containing the combined mesh data.
    */
    CObbStruct* readMultipleSTLs2CoppMesh(std::vector<std::string> filepathnames);

}

#endif //CHASSIS_IMPORT_HANDLER