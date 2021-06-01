#ifndef ROBOT_BUILDER
#define ROBOT_BUILDER

#include <vector>
#include "ChassisImportHandler.hpp"
#include "coppeliaGeometricRoutines-master/geom.h"

namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    class RobotBuilder {
        public:
        RobotBuilder(std::string robotCode);

        private:
        std::vector<CObbStruct*> staticObjects;
        std::vector<CObbStruct*> chassisConvexDecomp;
        std::vector<CObbStruct*> armSegments;
        CObbStruct* chassisConvexHull;
        std::vector<blueprintLine> blueprint;
        std::vector<std::vector<double>> DHParamsAt0 = {
            {60+(29.4/2),   -90,    0   ,   90},
            {66         ,   -90,    -425,   180},
            {9          ,   0  ,    -392,   180},
            {53         ,   90 ,    0   ,  -90},
            {95         ,   0  ,    0   ,   90},
            {76         ,   0  ,    0   ,   0},
        };
    };

}


#endif //ROBOT_BUILDER