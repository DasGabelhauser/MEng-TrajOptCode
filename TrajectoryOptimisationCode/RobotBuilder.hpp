#ifndef ROBOT_BUILDER
#define ROBOT_BUILDER

#include <vector>
#include "ChassisImportHandler.hpp"
#include "MeshOps.hpp"
#include "TrajectoryOptimisationNLP.hpp"
#include "DenavitHartenberg.hpp"

namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    /*  Class that optimises the trajectories of the steps to build a specific robot identified
        by the robt's number. Unifinished.
    */
    class RobotBuilder {
        public:
        RobotBuilder(std::string robotCode);

        private:
        std::vector<Hull> staticObjects;
        std::vector<Hull> chassisConvexDecomp;
        std::vector<Hull> armSegments;
        Hull chassisConvexHull;
        std::vector<blueprintLine> blueprint;
        //DH paramaters of the UR5 robot arm at 0 joint positions.
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