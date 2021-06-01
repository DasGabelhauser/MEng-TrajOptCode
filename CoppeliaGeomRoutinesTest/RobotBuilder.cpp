#include "/usr/local/include/coin-or/IpIpoptApplication.hpp"
#include "RobotBuilder.hpp"
#include "ChassisImportHandler.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>


namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    RobotBuilder::RobotBuilder(std::string robotCode) {

        // read in all necessary .stl files
        std::string modelsFilePath = "/home/eric/Documents/MENG-Robot-Fab-Project/Optimisation Code/Models/";
        std::string chassisFilePath = modelsFilePath + "Chassis/";
        //blueprint = readChassisBlueprint(chassisFilePath + "blueprints/BP" + robotCode + ".csv");
        //chassisConvexDecomp = importConvexDecomposedChassis(chassisFilePath + "meshesNoClipsJustSpace/mesh" + robotCode + ".stl", blueprint);
        //chassisConvexHull = importConvexHullChassis(chassisFilePath + "completeMeshes/mesh" + robotCode + ".stl");

        std::string staticPartsPath = modelsFilePath + "Static Parts Assembly/";
        std::ifstream staticPartsListFile(staticPartsPath+"staticPartsList.txt", std::ios::in);
        std::string staticPartFilename;
        std::vector<std::string> staticPolyhedronFilepathnames;
        
        std::getline(staticPartsListFile, staticPartFilename);
        while(!staticPartsListFile.eof()) {
            std::cout << "reading static part\n";
            staticObjects.push_back(readSTL2CoppMesh(staticPartsPath + staticPartFilename));
            std::getline(staticPartsListFile, staticPartFilename);
        }
        staticPartsListFile.close();


        //staticObjects.push_back(readSTL2CoppMesh(staticPartsPath + "Static Parts Assembly.stl"));

        std::string armSegmentsFilePath = modelsFilePath + "ArmSegments/intelligentConvexHulls/";
        for (int armSegment = 1; armSegment < 6; armSegment++) {
            armSegments.push_back(readSTL2CoppMesh(armSegmentsFilePath + "Link" + std::to_string(armSegment) + "Cube.stl"));
            std::cout << "reading arm segment\n";
            //armSegments.push_back(readSTL2CoppMesh(staticPartsPath + "Static Parts Assembly_RobotFabFrame_"+std::to_string(armSegment)+"_"+std::to_string(armSegment)+".stl"));
        }

        std::cout << "done reading stls\n";

        //set up time measurement variables
        auto start = std::chrono::steady_clock::now();
        auto end = std::chrono::steady_clock::now();
        auto time = end-start;

        //mesh distance caching variables to speed up calculation
        int mesh1Caching [6][40];
        int mesh2Caching [6][40];

        //find distances between all collision pairs at various positions.
        //time how long it takes to calculate distances
        for (double t = 0; t < 1; t = t + 0.001) {
            C7Vector mesh1Transformation;
            mesh1Transformation.X=C3Vector(0.5*cos(0.1*t),0.0,0.0); // sinusoidal x-axis movement
            mesh1Transformation.Q.setIdentity();
            C7Vector mesh2Transformation;
            mesh2Transformation.X=C3Vector(0.0,0.4*cos(0.2*t),0.3*cos(0.1*t)); // sinusoidal y- and z-axis movement
            C3Vector eulerAngles(0.0,0.0,t*1.2); // rotation around vertical axis
            mesh2Transformation.Q.setEulerAngles(eulerAngles);
            
            start = std::chrono::steady_clock::now();
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < (int)staticObjects.size(); j++) {
                    //std::cout << "arm segment: " << i+1 << ", staticObj: " << j << "\n";
                    simReal dist=REAL_MAX;
                    geom_getMeshMeshDistanceIfSmaller(
                        armSegments[i],mesh1Transformation,staticObjects[j],mesh2Transformation,dist,
                        nullptr, nullptr, &mesh1Caching[i][j], &mesh2Caching[i][j]);
            
                    //printf("mesh-mesh minimum distance: %f\n",dist);
                }
            }
            //calculate how long distance calculations took
            end = std::chrono::steady_clock::now();
            time = end-start;
            std::cout << "time: " << time.count() << "\n";
            //calculate how long it takes to calculate time (about 30ns), as a control measure
            start = std::chrono::steady_clock::now();
            end = std::chrono::steady_clock::now();
            time = end-start;
            std::cout << "time: " << time.count() << "\n";
        }
    }
}