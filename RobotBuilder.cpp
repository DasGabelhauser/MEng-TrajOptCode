#include "/usr/local/include/coin-or/IpIpoptApplication.hpp"
#include "RobotBuilder.hpp"
#include "ChassisImportHandler.hpp"
#include "MeshOps.hpp"
#include "TrajectoryOptimisationNLP.hpp"
#include "DenavitHartenberg.hpp"
#include <string>


namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    RobotBuilder::RobotBuilder(std::string robotCode) {

        // read in all necessary .stl files
        std::string modelsFilePath = "/home/eric/Documents/MENG-Robot-Fab-Project/Optimisation Code/Models/";
        std::string chassisFilePath = modelsFilePath + "Chassis/";
        blueprint = readChassisBlueprint(chassisFilePath + "blueprints/BP" + robotCode + ".csv");
        chassisConvexDecomp = importConvexDecomposedChassis(chassisFilePath + "meshesNoClipsJustSpace/mesh" + robotCode + ".stl", blueprint);
        chassisConvexHull = importConvexHullChassis(chassisFilePath + "completeMeshes/mesh" + robotCode + ".stl");

        std::string staticPartsPath = modelsFilePath + "Static Parts Assembly/";
        std::ifstream staticPartsListFile(staticPartsPath+"staticPartsList.txt", std::ios::in);
        std::string staticPartFilename;

        std::getline(staticPartsListFile, staticPartFilename);
        while(!staticPartsListFile.eof()) {
            staticObjects.push_back(readBinaryStl(staticPartsPath + staticPartFilename));
            std::getline(staticPartsListFile, staticPartFilename);
        }
        staticPartsListFile.close();

        //load in higher quality arm segment models
        /*std::string armSegmentsFilePath = modelsFilePath + "ArmSegments/convexHulls/";
        for (int armSegment = 1; armSegment < 7; armSegment++) {
            //armSegments.push_back(readBinaryStl(armSegmentsFilePath + "UR5_link_" + std::to_string(armSegment) + "_CH.stl"));
            armSegments.push_back(readBinaryStl(staticPartsPath + "Static Parts Assembly_RobotFabFrame_"+std::to_string(armSegment)+"_"+std::to_string(armSegment)+".stl"));
        }*/

        //load in lower quality arm segment models
        std::string armSegmentsFilePath = modelsFilePath + "ArmSegments/intelligentConvexHulls/";
        for (int armSegment = 1; armSegment < 7; armSegment++) {
            armSegments.push_back(readBinaryStl(armSegmentsFilePath + "Link"+std::to_string(armSegment)+"Cube.stl"));
        }

        //export hulls for verification purposes
        /*Hull totalHull = *new Hull();
        for (std::vector<Hull>::size_type subHull = 0; subHull < chassisConvexDecomp.size(); subHull++) {
            totalHull = totalHull + chassisConvexDecomp[subHull];
        }

        Hull totalStatic = *new Hull();
        for (std::vector<Hull>::size_type subHull = 0; subHull < staticObjects.size(); subHull++) {
            totalStatic = totalStatic + staticObjects[subHull];
        }

        totalHull.exportBinaryStl(modelsFilePath + "tempTotalHull.stl");
        chassisConvexHull.exportBinaryStl(modelsFilePath + "tempConvexHull.stl");
        totalStatic.exportBinaryStl(modelsFilePath + "tempStatics.stl");*/

        //create IPOPT app
        // Create a new instance of IpoptApplication
        //  (use a SmartPtr, not raw)
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        
        // Change some options
        // Note: The following choices are only examples, they might not be
        //       suitable for your optimization problem.
        app->Options()->SetNumericValue("tol", 1e-7);
        //app->Options()->SetNumericValue("acceptable_tol", 1e-6);
        app->Options()->SetIntegerValue("print_level", 3);
        app->Options()->SetIntegerValue("print_frequency_iter", 1000);
        app->Options()->SetIntegerValue("max_iter", 60000);
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("output_file", "ipopt.out");
        app->Options()->SetStringValue("linear_solver", "mumps");
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");

        //enable derivative checker
        app->Options()->SetStringValue("derivative_test", "first-order");
        app->Options()->SetStringValue("derivative_test_print_all", "no");
        // The following overwrites the default name (ipopt.opt) of the options file
        // app->Options()->SetStringValue("option_file_name", "hs071.opt");
        
        // Initialize the IpoptApplication and process the options
        Ipopt::ApplicationReturnStatus status;
        status = app->Initialize();

        /* create a test trajectory optimisation NLP and solve it. */
        double startPos[6] = {0.0,0.0,0.0,0.0,0.0,-360.0};
        double endPos[6] = {300, 400, 250, -60, 100, 45};
        int numSegs = 5;
        std::vector<std::vector<double>> returnX;
        Ipopt::SmartPtr<Ipopt::TNLP> trajOptNLP = new  TrajectoryOptimisationNLP(
            numSegs, 
            startPos, 
            endPos,
            staticObjects,
            *new std::vector<Hull>,
            armSegments,
            DHParamsAt0,
            returnX,
            360,
            180,
            40
        );
        status = app->OptimizeTNLP(trajOptNLP);
    }

}