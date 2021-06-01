#include "MeshOps.hpp"
#include "/usr/local/include/coin-or/IpIpoptApplication.hpp"
#include "MeshDifferenceQP.hpp"
#include "TrajectoryOptimisationNLP.hpp"
#include "RobotBuilder.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <array>
using namespace Ipopt;

std::vector<double> rotTransMat2PosNorm(
        std::vector<std::vector<double>> rotMat,
        std::vector<double> transMat
    ) {
        std::vector<double> posNormVec = *new std::vector<double>(transMat);

        //multiplying the rotation vector by [1;0;0] and [0;1;0] is just selecting the 0th and 1st columns
        for(int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                posNormVec.push_back(rotMat[j][i]);
            }
        }
        return posNormVec;
    }

std::vector<double> pose2PosNorm(
        double posX, 
        double posY, 
        double posZ, 
        double rotX, 
        double rotY, 
        double rotZ
    ) {
        std::vector<std::vector<double>> rotMat = EricsTrajOpt::makeRotMat(rotX, rotY, rotZ);
        std::vector<double> transVec = {posX, posY, posZ};
        return rotTransMat2PosNorm(rotMat, transVec);
    }

int main(int, char**) {

    
    std::vector<double> posNorm = pose2PosNorm(500,400,250,-60,100,45);
    for(int i = 0; i < 9; i++) {
        std::cout << posNorm[i] << "\n";
    }
    //launch a RobotBuilder to build a robot. unfinished.
    EricsTrajOpt::RobotBuilder("3134");

    //old development contents of launcher, kept around in case they are needed to test something
///////////////////////

    /*// Create a new instance of your nlp
    //  (use a SmartPtr, not raw)
    double startPos[6] = {0.0,0.0,0.0,0.0,0.0,-360.0};
    double endPos[6] = {60.0, 120.0, 180.0, 240.0, 300.0, 360.0};
    int numSegs;
    std::cin >> numSegs;
    SmartPtr<TNLP> mynlp = new  TrajectoryOptimisationNLP(
                                    numSegs, 
                                    startPos, 
                                    endPos,
                                    360,
                                    180,
                                    40
                                );
    
    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    
    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    app->Options()->SetNumericValue("tol", 1e-7);
    app->Options()->SetIntegerValue("print_level", 3);
    app->Options()->SetIntegerValue("max_iter", 6000);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("linear_solver", "ma27");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    //enable derivative checker
    app->Options()->SetStringValue("derivative_test", "first-order");
    app->Options()->SetStringValue("derivative_test_print_all", "no");
    // The following overwrites the default name (ipopt.opt) of the options file
    // app->Options()->SetStringValue("option_file_name", "hs071.opt");
    
    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded )
    {
       std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
       return (int) status;
    }
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if( status == Solve_Succeeded )
    {
       std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else
    {
       std::cout << std::endl << std::endl << "*** The problem FAILED! Status code: " << status << std::endl;
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.

    return (int) status;*/

///////////////////////

    /*std::string filePath = "/home/eric/Documents/MENG-Robot-Fab-Project/Models/OpenSCAD/MESH_DISTANCE_TESTS/";
    std::string fileName1 = "cube_10_35.stl";
    std::string fileName2 = "cube_10_neg35.stl";*/

    /*std::string filePath = "/home/eric/Documents/MENG-Robot-Fab-Project/Optimisation Code/Models/Chassis/convexClips/";
    std::string fileName1 = "slide_clip_CH.stl";
    std::string fileName2 = "slide_clip_CH.stl";*/
    
    /*EricsTrajOpt::Hull hull1 = EricsTrajOpt::readBinaryStl(filePath+fileName1);
    EricsTrajOpt::Hull hull2 = EricsTrajOpt::readBinaryStl(filePath+fileName2);
    hull1.setHullPermPosRot(0,0,-100,0,0,0);
    // Create a new instance of your nlp
    //  (use a SmartPtr, not raw)

    SmartPtr<TNLP> mynlp = new Mesh_Difference_QP(&hull1, &hull2);

    
    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    
    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    app->Options()->SetNumericValue("tol", 1e-7);
    app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("jac_d_constant", "yes");
    app->Options()->SetStringValue("jac_c_constant", "yes");
    app->Options()->SetStringValue("linear_solver", "mumps");
    // The following overwrites the default name (ipopt.opt) of the options file
    // app->Options()->SetStringValue("option_file_name", "hs071.opt");
    
    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded )
    {
       std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
       return (int) status;
    }
    
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if( status == Solve_Succeeded )
    {
       std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else
    {
       std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.

    return (int) status;*/
    
}