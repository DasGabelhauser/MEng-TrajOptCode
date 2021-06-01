// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
//modified by Eric Walker

#include "TrajectoryOptimisationNLP.hpp"
#include "/usr/local/include/coin-or/IpSolveStatistics.hpp"

#include <cassert>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "matplotlibcpp.h"
#include <chrono>


#define numJoints 6
#define x_VALUE(joint,controlPoint) x[(controlPoint+3) * numJoints + joint]
#define x_INDEX(joint,controlPoint) ((controlPoint+3) * numJoints + joint)
#define duration x[n-1]
#define diffDiff 1e-8 //fixed difference to change variables by when approximating derivatives
#define obstacleClearance 3 //minimum distance in mm between obstacles

using namespace Ipopt;

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

enum splineCoeffIndices {A = 0, B = 1, C = 2, D = 3, deltaT = 4};
enum controlPointIndices {C0 = 3, Cneg1 = 2, Cneg2 = 1, Cneg3 = 0};

// constructor
TrajectoryOptimisationNLP::TrajectoryOptimisationNLP(
    int numIntervals, 
    double * startPose, 
    double * endPose, 
    std::vector<EricsTrajOpt::Hull> staticObstacles,
    std::vector<EricsTrajOpt::Hull> gripperObjects,
    std::vector<EricsTrajOpt::Hull> armLinks,
    std::vector<std::vector<double>> DHParamsAt0,
    std::vector<std::vector<double>> &returnX,
    double maxPos/* = 360*/, 
    double maxVel/* = 180*/, 
    double maxAcc/* = 40*/
    ) 
    :   m_numIntervals{numIntervals},  
        m_numControlPoints{numIntervals + 3}, 
        m_startPose{startPose[0], startPose[1], startPose[2], startPose[3], startPose[4], startPose[5]}, 
        m_endPose{endPose[0], endPose[1], endPose[2], endPose[3], endPose[4], endPose[5]}, 
        m_staticObstacles{staticObstacles},
        m_gripperObjects{gripperObjects},
        m_armLinks{armLinks},
        m_DHParamsAt0{DHParamsAt0},
        m_returnX{returnX},
        m_maxPos{maxPos}, 
        m_maxVel{maxVel}, 
        m_maxAcc{maxAcc}
    {

    //populate staticObj_2_ArmLinks with Mesh_Difference_QPs to find
    //distance between pairs of static objects and arm links
    for (int armLink = 1; armLink < 6; armLink++) {
        for (int staticObject = 0; staticObject < (int)m_staticObstacles.size(); staticObject++) {
            staticObj_2_ArmLinks[armLink-1].push_back(
                new Mesh_Difference_QP( &(m_armLinks[armLink]),
                                        &(m_staticObstacles[staticObject]))              
                );
        }
    }

    //populate link1_2_ArmLinks with Mesh_Difference_QPs to find
    //distance between arm link 1 and arm links 3-6
    for (int armLink = 0; armLink < 4; armLink++) {
        link1_2_ArmLinks [armLink] = 
            new Mesh_Difference_QP( &(m_armLinks[0]),
                                    &(m_armLinks[armLink+2]));              
    }
    //populate link2_2_ArmLinks with Mesh_Difference_QPs to find
    //distance between arm link 2 and arm links 4-6
    for (int armLink = 0; armLink < 3; armLink++) {
        link2_2_ArmLinks [armLink] = 
            new Mesh_Difference_QP( &(m_armLinks[1]),
                                    &(m_armLinks[armLink+3]));              
    }
    //populate link3_2_ArmLinks with Mesh_Difference_QPs to find
    //distance between arm link 3 and arm links 5-6
    for (int armLink = 0; armLink < 2; armLink++) {
        link3_2_ArmLinks [armLink] = 
            new Mesh_Difference_QP( &(m_armLinks[2]),
                                    &(m_armLinks[armLink+4]));              
    }

    //populate gripperObj_2_ArmLinks with Mesh_Difference_QPs to find
    //distance between pairs of gripper objects and arm links
    for (int armLink = 0; armLink < 5; armLink++) {
        for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
            gripperObj_2_ArmLinks[armLink].push_back(
                new Mesh_Difference_QP( &(m_armLinks[armLink]),
                                        &(m_gripperObjects[gripperObj]))  
            );
        }
    }

    //populate gripperObj_2_ArmLinks with Mesh_Difference_QPs to find
    //distance between pairs of gripper objects and static objects
    for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
        gripperObj_2_staticObj.push_back(*new std::vector<Ipopt::SmartPtr<Ipopt::TNLP>>);
        for (int staticObject = 0; staticObject < (int)m_staticObstacles.size(); staticObject++) {
            gripperObj_2_staticObj[gripperObj].push_back(
                new Mesh_Difference_QP( &(m_gripperObjects[gripperObj]),
                                        &(m_staticObstacles[staticObject]))              
                );
        }
    }
    
    // Create a new instance of IpoptApplication
    //  to solve obstacle distance calculations with
    obstDistSolver = IpoptApplicationFactory();
    
    // Change some options
    obstDistSolver->Options()->SetNumericValue("tol", 1e-7);
    obstDistSolver->Options()->SetIntegerValue("print_level", 0);
    obstDistSolver->Options()->SetStringValue("mu_strategy", "adaptive");
    obstDistSolver->Options()->SetStringValue("output_file", "ipopt.out");
    obstDistSolver->Options()->SetStringValue("jac_d_constant", "yes");
    obstDistSolver->Options()->SetStringValue("jac_c_constant", "yes");
    obstDistSolver->Options()->SetStringValue("linear_solver", "mumps");
    
    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = obstDistSolver->Initialize();
    if( status != Solve_Succeeded )
    {
       std::cout << std::endl << std::endl << "*** Error during obstDistSolver initialization!" << std::endl;
    }
}

// destructor
TrajectoryOptimisationNLP::~TrajectoryOptimisationNLP()
{ }

// [TNLP_get_nlp_info]
// returns the size of the problem
bool TrajectoryOptimisationNLP::get_nlp_info(
    Index&          n,
    Index&          m,
    Index&          nnz_jac_g,
    Index&          nnz_h_lag,
    IndexStyleEnum& index_style
)
{
    //std::cout << "get_nlp_info\n";
    // The problem has number of design variables equal to 
    // one spline coefficient per joint, per control point plus a duration
    // [A1_-3,      A2_-3       ... A6_-3     ],
    // [A1_-2,      A2_-2       ... A6_-2     ],
    // [A1_-1,      A2_-1       ... A6_-1     ],
    // [A1_0,       A2_0        ... A6_0      ],
    //          ...
    // [A1_numSegs, A2_numSegs  ... A6_numSegs],
    // [T]
    // (in vector form, not a matrix)

    n = m_numControlPoints * numJoints + 1;
 
    /*  constraints table
        num                                     |   purpose
        1                                       |   keeping durations > 0 
        numJoints * numIntervals                |   max pos[joint][interval] < limit
        numJoints * numIntervals                |   min pos[joint][interval] > limit
        numJoints * numIntervals                |   max vel[joint][interval] < limit
        numJoints * numIntervals                |   min vel[joint][interval] > limit
        numJoints * numIntervals                |   max acc[joint][interval] < limit
        numJoints * numIntervals                |   min acc[joint][interval] > limit
        6                                       |   ensure startPos are correct
        6                                       |   ensure startVel are 0 
        9                                       |   ensure endPos are correct
        6                                       |   ensure endVel are 0 
        5 * numIntervals                        |   collision constraints for static obstacles->link2-6
        4 * numIntervals                        |   collision constraints for link1->link3-6
        3 * numIntervals                        |   collision constraints for link2->link4-6
        2 * numIntervals                        |   collision constraints for link3->link5-6
        6 * numIntervals*size(m_gripperObjects) |   collision constraints for each gripper object->link1-5, static objects     
    */
            m = (1) + 
        (numJoints * m_numIntervals) +
        (numJoints * m_numIntervals) +
        (numJoints * m_numIntervals) +
        (numJoints * m_numIntervals) +
        (numJoints * m_numIntervals) +
        (numJoints * m_numIntervals) +
        (6) +
        (6) +
        (9) +
        (6) +
        ((5+4+3+2) * m_numIntervals) + 
        (6*m_numIntervals*m_gripperObjects.size());
 
    /*  nnz in jacobian is:
        1 from duration constraint
        numIntervals * numJoints * 4 
            for each of the min/max pos constraints (so *2)
            (as each segment's constraints are only affected by the preceeding 4 control points)
        numIntervals * numJoints * 5 
            for each of the min/max vel/acc constraints (so *4)
            (as each segment's constraints are only affected by the preceeding 4 control points and the duration)
        6*3 
            as first 3 control points can affect startPoss (the 4th, C0, does not)
        6*3 
            as C-3 and C-1 + duration can affect startVels
        9*3*6 
            as last 3 control points can affect endPoss (d/dC-3 = 0)
        6*3 
            as C-2 and C0 + duration can affect endVels
        (2+3+4+5+6)*numIntervals*4 
            for the static obstacles->link2-6 constraints (4 control points for each joint that
            can affect each link, for each interval)
        (2+3+4+5)*numIntervals*4
            for the link1->link3-6 constraints (4 control points for each joint that
            can affect each link, for each interval. The first joint cannot affect the distance)
        (2+3+4)*numIntervals*4
            for the link2->link4-6 constraints (4 control points for each joint that
            can affect each link, for each interval. The first two joints cannot affect the distance)
        (2+3)*numIntervals*4
            for the link3->link5-6 constraints (4 control points for each joint that
            can affect each link, for each interval. The first three joints cannot affect the distance)
        (6+5+4+3+2+1)*numIntervals*4*size(m_gripperObjects)
            for the gripper object->static objects,link1-5 constraints.
    */
    nnz_jac_g = (1) +
                (2 * m_numIntervals * numJoints * 4) +
                (4 * m_numIntervals * numJoints * 5) +
                (6*3) +
                (6*3) + 
                (9*3*6) +
                (6*3) + 
                ((2+3+4+5+6)*m_numIntervals*4) + 
                ((2+3+4+5)*m_numIntervals*4) + 
                ((2+3+4)*m_numIntervals*4) + 
                ((2+3)*m_numIntervals*4) + 
                ((6+5+4+3+2+1)*m_numIntervals*4*m_gripperObjects.size());
 
    // the Hessian is gonna be very dense. We only need the lower-left entries, however.
    //nnz_h_lag = (30*m_numIntervals)+(9*6)+1;
    nnz_h_lag = (n*(n+1))/2;
 
    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;
    //std::cout << "get_nlp_info done\n";
    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool TrajectoryOptimisationNLP::get_bounds_info(
    Index   n,
    Number* x_l,
    Number* x_u,
    Index   m,
    Number* g_l,
    Number* g_u
) {
    //std::cout << "get_bounds_info\n";
    for (int i = 0; i < n; i++) {
        //no constraints on variables
        x_l[i] = -2e19;
        x_u[i] = 2e19;
    }

    // constraint on duration
    g_l[0] = 0;
    g_u[0] = 2e19;

    int offset = numJoints * m_numIntervals;
    for (int i = 0; i < numJoints * m_numIntervals; i++) {
        //set limits for max pos
        g_l[i + (offset * 0) + 1] = -2e19;
        g_u[i + (offset * 0) + 1] = m_maxPos;
        //set limits for min pos
        g_l[i + (offset*1) + 1] = -m_maxPos;
        g_u[i + (offset*1) + 1] = 2e19;
        //set limits for max vel
        g_l[i + (offset*2) + 1] = -2e19;
        g_u[i + (offset*2) + 1] = m_maxVel;
        //set limits for min vel
        g_l[i + (offset*3) + 1] = -m_maxVel;
        g_u[i + (offset*3) + 1] = 2e19;
        //set limits for max acc
        g_l[i + (offset*4) + 1] = -2e19;
        g_u[i + (offset*4) + 1] = m_maxAcc;
        //set limits for min acc
        g_l[i + (offset*5) + 1] = -m_maxAcc;
        g_u[i + (offset*5) + 1] = 2e19;

    }

    offset = 1 + (6 * numJoints * m_numIntervals);
    // set limits for start position and start/end velocity
    for (int i = 0; i < 6; i++) {
        g_l[i + offset + 0 ] = m_startPose[i];
        g_u[i + offset + 0 ] = m_startPose[i];
        g_l[i + offset + 6 ] = 0;
        g_u[i + offset + 6 ] = 0;

        g_l[i + offset + 21] = 0;
        g_u[i + offset + 21] = 0;
    }

    //convert end position from pose to PosNorm form
    std::vector<double> endPosNorm = pose2PosNorm(
        m_endPose[0],
        m_endPose[1],
        m_endPose[2],
        m_endPose[3],
        m_endPose[4],
        m_endPose[5]
    );

    //set limits for end poistion in PosNorm form
    for (int i = 0; i < 9; i++) {
        g_l[i + offset + 12] = endPosNorm[i];
        g_u[i + offset + 12] = endPosNorm[i];
    }

    offset = 1 + (6 * numJoints * m_numIntervals) + 27;
    // set limits for obstacle constraints
    for (   int i = 0; 
            i < ((5+4+3+2) * m_numIntervals) + 
                (6*m_numIntervals*(int)m_gripperObjects.size());
            i++) {
        //all obstacle constraints can be no closer than obstacleClearance,
        //and infinitely distant
        g_l[offset + i] = obstacleClearance;
        g_u[offset + i] = 2e19;
    }
    //std::cout << "get_bounds_info done\n";
    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool TrajectoryOptimisationNLP::get_starting_point(
    Index   n,
    bool    init_x,
    Number* x,
    bool    init_z,
    Number* z_L,
    Number* z_U,
    Index   m,
    bool    init_lambda,
    Number* lambda) {
    //std::cout << "get_starting_point\n";
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
 
    // find the largest angle from the start position to the end position
    // TODO do something to get endPose in joint space from endPose in tool space
    /*double endPose[6] = {-56.6551, 29.8349, 73.8019, 83.0549, 131.832, 79.096};
    double largestAngle = 0;
    double angleChanges[6];
    for (int i = 0; i < numJoints; i++) {
        angleChanges[i] = endPose[i] - m_startPose[i];
        if (std::abs(angleChanges[i]) > largestAngle) {
            largestAngle = std::abs(angleChanges[i]);
        }
    }*/
    // find time taken to move largest angle with no acceleration and maximum velocity.
    // x[n-1] is the last entry of x, and is the variable representing the duration.
    //double totalDuration = largestAngle/m_maxVel;
    // find initial velocities
    /*double initialVelocities[6];
    for (int i = 0; i < numJoints; i++) {
        initialVelocities[i] = angleChanges[i]/totalDuration;
    }*/
    // calculate the spline coefficients to give the desired initial trajectories
    // each joint follows the pattern 
    // [startPos-vel, startPos, startPos+vel, startPos+2*vel ... startPos+(numIntervals-1)*vel]
    for (int i = -3; i < m_numIntervals; i++) {
        for (int joint = 0; joint < numJoints; joint++) {
            x_VALUE(joint,i) = 0;//m_startPose[joint] + (i+2)*initialVelocities[joint];
        }
    }
    x[n-1] = 10;//totalDuration*2;

    
    /*if (getStartPositionCounter == 1) {
        graphSolution(x, n);
    }
    getStartPositionCounter++;*/
    //std::cout << "get_starting_point done\n";
    return true;

}
// [TNLP_get_starting_point]

// [TNLP_eval_f]
// returns the value of the objective function
bool TrajectoryOptimisationNLP::eval_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number&       obj_value) {
    //std::cout << "eval_f\n";
    //objective function is just the duration, stored in x[n-1].
    obj_value = duration;
    return true;
    //std::cout << "eval_f done\n";
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool TrajectoryOptimisationNLP::eval_grad_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number*       grad_f) {
    //std::cout << "eval_grad_f\n";
    // the coefficients of the spline do nothing to do with f
    for (int i = 0; i < numJoints*m_numControlPoints; i++) {
        // set df/d(splineCoeff) to 0
        grad_f[i] = 0;
    }
    //f = x[n-1], so df/dx[n-1] = 1
    grad_f[n-1] = 1;
    //std::cout << "eval_grad_f done\n";
    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool TrajectoryOptimisationNLP::eval_g(
    Index         n,
    const Number* x,
    bool          new_x,
    Index         m,
    Number*       g
)
{
    std::cout << "eval_g\n";
    //variables to hold spline parameters
    double splineParameters[5];
    //calculate duration constraint
    g[0] = duration;

    //calculate min/max pos/vel/acc for each joint and interval
    int offset = 1;
    double maxPos, minPos, maxVel, minVel, maxAcc, minAcc;
    for (int interval = 0; interval < m_numIntervals; interval++) {
        for (int joint = 0; joint < numJoints; joint++) {
            findIntervalCoefficients(x, n, joint, interval, splineParameters);
            findMaxMinPos(splineParameters, &maxPos, &minPos);
            findMaxMinVel(splineParameters, &maxVel, &minVel);
            findMaxMinAcc(splineParameters, &maxAcc, &minAcc);
            g[(0 * m_numIntervals + interval)*numJoints + joint + offset] = maxPos;
            g[(1 * m_numIntervals + interval)*numJoints + joint + offset] = minPos;
            g[(2 * m_numIntervals + interval)*numJoints + joint + offset] = maxVel;
            g[(3 * m_numIntervals + interval)*numJoints + joint + offset] = minVel;
            g[(4 * m_numIntervals + interval)*numJoints + joint + offset] = maxAcc;
            g[(5 * m_numIntervals + interval)*numJoints + joint + offset] = minAcc; 
        }
    }

    //calculate start/end pos/vel
    offset = 1 + (6 * numJoints * m_numIntervals);
    double endPosJSpace[6] = {0,0,0,0,0,0};
    for (int joint = 0; joint < numJoints; joint++) {
        //caculate start pos/vel
        findIntervalCoefficients(x, n, joint, 0, splineParameters);
        g[offset + joint + (numJoints * 0)] = calculatePos(splineParameters, 0);
        g[offset + joint + (numJoints * 1)] = calculateVel(splineParameters, 0);
        //calculate end pos/vel
        findIntervalCoefficients(x, n, joint, m_numIntervals-1, splineParameters);
        endPosJSpace[joint] = calculatePos(splineParameters, 1);
        g[offset + joint + (numJoints * 2) + 9] = calculateVel(splineParameters, 1);
    }
    std::vector<EricsTrajOpt::DenavitHartenberg> endPosDH = getDHAtJointPose(m_DHParamsAt0, endPosJSpace);
    std::vector<double> endPosPosNorm = DH2PosNorm(endPosDH[5]);
    for (int i = 0; i < 9; i++) {
        g[offset + (numJoints * 2) + i] = endPosPosNorm[i];
    }

    //get offset index of first obstacle distance constraint
    offset = 1 + (6 * numJoints * m_numIntervals) + 6+6+9+6;
    //calculate obstacle distance constraints
    calcAllCollisions(x,n,g+offset);

    std::cout << "eval_g done\n";
    return true;
}
// [TNLP_eval_g]

// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool TrajectoryOptimisationNLP::eval_jac_g(
    Index         n,
    const Number* x,
    bool          new_x,
    Index         m,
    Index         nele_jac,
    Index*        iRow,
    Index*        jCol,
    Number*       values
) {
    //variables to store offsets
    int jacOffset = 0;
    int gOffset = 0;

    //variable to store spline control points to differentiate W.R.T.
    double controlPoints[4];
    //variables to store differentiated max/min pos/vel/acc
    double dMaxPos[4];
    double dMinPos[4];
    double dMaxVel[5];
    double dMinVel[5];
    double dMaxAcc[5];
    double dMinAcc[5];

    if( values == NULL )
    {
        std::cout << "eval_jac_g first time\n";
        // copy the structure of the Jacobian

        //constraint on time
        iRow[0] = 0;
        jCol[0] = n-1; //duration is last constraint

        //constraints on max/min pos
        //set offsets
        jacOffset = 1;
        gOffset = 1;
        for (int i = 0; i < 2; i++) { // 2 loops for min/max pos constraints
            for (int interval = 0; interval < m_numIntervals; interval++) { //m_numIntervals for each interval's constraints
                for (int joint = 0; joint < numJoints; joint++) { //numJoints for each joint's constraints
                    for (int j = 0; j < 4; j++) { //4 variables can affect each constraint
                        //((i*m_numIntervals + interval)*numJoints + joint)*4 + jacOffset + j
                        iRow[(i*m_numIntervals*numJoints*4) + (interval*numJoints*4) + (joint*4) + j + jacOffset] = 
                            (i*m_numIntervals*numJoints) + (interval * numJoints) + joint + gOffset;
                        jCol[(i*m_numIntervals*numJoints*4) + (interval*numJoints*4) + (joint*4) + j + jacOffset] = 
                            x_INDEX(joint, interval - 3 + j);
                    }
                }
            }
        }

        //constraints on max/min vel/acc
        //increase the offset
        jacOffset += 2*m_numIntervals*numJoints*4;
        gOffset += 2*m_numIntervals*numJoints;
        for (int i = 0; i < 4; i++) { // 4 loops for min/max vel/acc constraints
            for (int interval = 0; interval < m_numIntervals; interval++) { //m_numIntervals for each interval's constraints
                for (int joint = 0; joint < numJoints; joint++) { //numJoints for each joint's constraints
                    for (int j = 0; j < 4; j++) { //5 variables can affect each constraint, but the 5th is not iterable (duration)
                        iRow[(i*m_numIntervals*numJoints*5) + (interval*numJoints*5) + (joint*5) + j + jacOffset] = 
                            (i*m_numIntervals*numJoints) + (interval * numJoints) + joint + gOffset;
                        jCol[(i*m_numIntervals*numJoints*5) + (interval*numJoints*5) + (joint*5) + j + jacOffset] = 
                            x_INDEX(joint, interval - 3 + j);
                    }
                    //min/max vel/acc constraints are affected by duration, add to sparsity structure
                    iRow[(i*m_numIntervals*numJoints*5) + (interval*numJoints*5) + (joint*5) + 4 + jacOffset] = 
                        (i*m_numIntervals*numJoints) + (interval * numJoints) + joint + gOffset;
                    jCol[(i*m_numIntervals*numJoints*5) + (interval*numJoints*5) + (joint*5) + 4 + jacOffset] = n-1;
                }
            }
        }

        //increase the offset
        jacOffset += 4*m_numIntervals*numJoints*5;
        gOffset += 4*m_numIntervals*numJoints;
        //start/end Pos/Vel constraints
        //derivatives are ordered C-3, C-2, C-1, C0, duration, with elements missing as appropriate
        for (int joint = 0; joint < numJoints; joint++) { //numJoints for each joint's startPos, startVel, endPos and endVel
            //startPos constraint derivatives
            iRow[jacOffset + (numJoints*3)*0 + joint*3 + 0] = gOffset + joint;
            jCol[jacOffset + (numJoints*3)*0 + joint*3 + 0] = x_INDEX(joint, -3);
            iRow[jacOffset + (numJoints*3)*0 + joint*3 + 1] = gOffset + joint;
            jCol[jacOffset + (numJoints*3)*0 + joint*3 + 1] = x_INDEX(joint, -2);
            iRow[jacOffset + (numJoints*3)*0 + joint*3 + 2] = gOffset + joint;
            jCol[jacOffset + (numJoints*3)*0 + joint*3 + 2] = x_INDEX(joint, -1);

            //startVel
            iRow[jacOffset + (numJoints*3)*1 + joint*3 + 0] = gOffset + 6 + joint;
            jCol[jacOffset + (numJoints*3)*1 + joint*3 + 0] = x_INDEX(joint, -3);
            iRow[jacOffset + (numJoints*3)*1 + joint*3 + 1] = gOffset + 6 + joint;
            jCol[jacOffset + (numJoints*3)*1 + joint*3 + 1] = x_INDEX(joint, -1);
            iRow[jacOffset + (numJoints*3)*1 + joint*3 + 2] = gOffset + 6 + joint;
            jCol[jacOffset + (numJoints*3)*1 + joint*3 + 2] = n-1; //duration is at x[n-1]

            //endPos
            //endPos is defined in PosNorm form in tool space, so all 9 constraints are affected by all
            //numJoints*3 variables
            for (int i = 0; i < 9; i++) {
                iRow[jacOffset + (numJoints*3)*(2+i) + joint*3 + 0] = gOffset + 12 + i;
                jCol[jacOffset + (numJoints*3)*(2+i) + joint*3 + 0] = x_INDEX(joint, (m_numIntervals-1)-2);
                iRow[jacOffset + (numJoints*3)*(2+i) + joint*3 + 1] = gOffset + 12 + i;
                jCol[jacOffset + (numJoints*3)*(2+i) + joint*3 + 1] = x_INDEX(joint, (m_numIntervals-1)-1);
                iRow[jacOffset + (numJoints*3)*(2+i) + joint*3 + 2] = gOffset + 12 + i;
                jCol[jacOffset + (numJoints*3)*(2+i) + joint*3 + 2] = x_INDEX(joint, (m_numIntervals-1)-0);
            }

            //endVel
            iRow[jacOffset + (numJoints*3)*11 + joint*3 + 0] = gOffset + 21 + joint;
            jCol[jacOffset + (numJoints*3)*11 + joint*3 + 0] = x_INDEX(joint, (m_numIntervals-1)-2);
            iRow[jacOffset + (numJoints*3)*11 + joint*3 + 1] = gOffset + 21 + joint;
            jCol[jacOffset + (numJoints*3)*11 + joint*3 + 1] = x_INDEX(joint, (m_numIntervals-1)-0);
            iRow[jacOffset + (numJoints*3)*11 + joint*3 + 2] = gOffset + 21 + joint;
            jCol[jacOffset + (numJoints*3)*11 + joint*3 + 2] = n-1; //duration is at x[n-1]
        }

        //increase the offset
        jacOffset += (6*3) + (6*3) + (9*3*6) + (6*3);
        gOffset += 6+6+9+6;
        for (int interval = 0; interval < m_numIntervals; interval++) {
            for (int controlPoint = -3; controlPoint < 1; controlPoint++) {
                //there is no easy pattern to these constraints, so values are somewhat hardcoded
                iRow[jacOffset + 0] = gOffset+0;//static->L2
                jCol[jacOffset + 0] = x_INDEX(0,interval+controlPoint);
                iRow[jacOffset + 1] = gOffset+1;//static->L3
                jCol[jacOffset + 1] = x_INDEX(0,interval+controlPoint);
                iRow[jacOffset + 2] = gOffset+2;//static->L4
                jCol[jacOffset + 2] = x_INDEX(0,interval+controlPoint);
                iRow[jacOffset + 3] = gOffset+3;//static->L5
                jCol[jacOffset + 3] = x_INDEX(0,interval+controlPoint);
                iRow[jacOffset + 4] = gOffset+4;//static->L6
                jCol[jacOffset + 4] = x_INDEX(0,interval+controlPoint);
                
                iRow[jacOffset + 5] = gOffset+0;//static->L2
                jCol[jacOffset + 5] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 6] = gOffset+1;//static->L3
                jCol[jacOffset + 6] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 7] = gOffset+2;//static->L4
                jCol[jacOffset + 7] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 8] = gOffset+3;//static->L5
                jCol[jacOffset + 8] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 9] = gOffset+4;//static->L6
                jCol[jacOffset + 9] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 10] = gOffset+5;//L1->L3
                jCol[jacOffset + 10] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 11] = gOffset+6;//L1->L4
                jCol[jacOffset + 11] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 12] = gOffset+7;//L1->L5
                jCol[jacOffset + 12] = x_INDEX(1,interval+controlPoint);
                iRow[jacOffset + 13] = gOffset+8;//L1->L6
                jCol[jacOffset + 13] = x_INDEX(1,interval+controlPoint);

                iRow[jacOffset + 14] = gOffset+1;//static->L3
                jCol[jacOffset + 14] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 15] = gOffset+2;//static->L4
                jCol[jacOffset + 15] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 16] = gOffset+3;//static->L5
                jCol[jacOffset + 16] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 17] = gOffset+4;//static->L6
                jCol[jacOffset + 17] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 18] = gOffset+5;//L1->L3
                jCol[jacOffset + 18] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 19] = gOffset+6;//L1->L4
                jCol[jacOffset + 19] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 20] = gOffset+7;//L1->L5
                jCol[jacOffset + 20] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 21] = gOffset+8;//L1->L6
                jCol[jacOffset + 21] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 22] = gOffset+9;//L2->L4
                jCol[jacOffset + 22] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 23] = gOffset+10;//L2->L5
                jCol[jacOffset + 23] = x_INDEX(2,interval+controlPoint);
                iRow[jacOffset + 24] = gOffset+11;//L2->L6
                jCol[jacOffset + 24] = x_INDEX(2,interval+controlPoint);

                iRow[jacOffset + 25] = gOffset+2;//static->L4
                jCol[jacOffset + 25] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 26] = gOffset+3;//static->L5
                jCol[jacOffset + 26] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 27] = gOffset+4;//static->L6
                jCol[jacOffset + 27] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 28] = gOffset+6;//L1->L4
                jCol[jacOffset + 28] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 29] = gOffset+7;//L1->L5
                jCol[jacOffset + 29] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 30] = gOffset+8;//L1->L6
                jCol[jacOffset + 30] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 31] = gOffset+9;//L2->L4
                jCol[jacOffset + 31] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 32] = gOffset+10;//L2->L5
                jCol[jacOffset + 32] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 33] = gOffset+11;//L2->L6
                jCol[jacOffset + 33] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 34] = gOffset+12;//L3->L5
                jCol[jacOffset + 34] = x_INDEX(3,interval+controlPoint);
                iRow[jacOffset + 35] = gOffset+13;//L3->L6
                jCol[jacOffset + 35] = x_INDEX(3,interval+controlPoint);

                iRow[jacOffset + 36] = gOffset+3;//static->L5
                jCol[jacOffset + 36] = x_INDEX(4,interval+controlPoint);
                iRow[jacOffset + 37] = gOffset+4;//static->L6
                jCol[jacOffset + 37] = x_INDEX(4,interval+controlPoint);
                iRow[jacOffset + 38] = gOffset+7;//L1->L5
                jCol[jacOffset + 38] = x_INDEX(4,interval+controlPoint);
                iRow[jacOffset + 39] = gOffset+8;//L1->L6
                jCol[jacOffset + 39] = x_INDEX(4,interval+controlPoint);
                iRow[jacOffset + 40] = gOffset+10;//L2->L5
                jCol[jacOffset + 40] = x_INDEX(4,interval+controlPoint);
                iRow[jacOffset + 41] = gOffset+11;//L2->L6
                jCol[jacOffset + 41] = x_INDEX(4,interval+controlPoint);
                iRow[jacOffset + 42] = gOffset+12;//L3->L5
                jCol[jacOffset + 42] = x_INDEX(4,interval+controlPoint);
                iRow[jacOffset + 43] = gOffset+13;//L3->L6
                jCol[jacOffset + 43] = x_INDEX(4,interval+controlPoint);

                iRow[jacOffset + 44] = gOffset+4;//static->L6
                jCol[jacOffset + 44] = x_INDEX(5,interval+controlPoint);
                iRow[jacOffset + 45] = gOffset+8;//L1->L6
                jCol[jacOffset + 45] = x_INDEX(5,interval+controlPoint);
                iRow[jacOffset + 46] = gOffset+11;//L2->L6
                jCol[jacOffset + 46] = x_INDEX(5,interval+controlPoint);
                iRow[jacOffset + 47] = gOffset+13;//L3->L6
                jCol[jacOffset + 47] = x_INDEX(5,interval+controlPoint);

                for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                    iRow[jacOffset + 47 + (gripperObj*21) + 1] = 
                        gOffset+13+(gripperObj*6) + 1;//static->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 1] = 
                        x_INDEX(0,interval+controlPoint);

                    iRow[jacOffset + 47 + (gripperObj*21) + 2] = 
                        gOffset+13+(gripperObj*6) + 1;//static->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 2] = 
                        x_INDEX(1,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 3] = 
                        gOffset+13+(gripperObj*6) + 2;//L1->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 3] = 
                        x_INDEX(1,interval+controlPoint);

                    iRow[jacOffset + 47 + (gripperObj*21) + 4] = 
                        gOffset+13+(gripperObj*6) + 1;//static->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 4] = 
                        x_INDEX(2,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 5] = 
                        gOffset+13+(gripperObj*6) + 2;//L1->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 5] = 
                        x_INDEX(2,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 6] = 
                        gOffset+13+(gripperObj*6) + 3;//L2->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 6] = 
                        x_INDEX(2,interval+controlPoint);

                    iRow[jacOffset + 47 + (gripperObj*21) + 7] = 
                        gOffset+13+(gripperObj*6) + 1;//static->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 7] = 
                        x_INDEX(3,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 8] = 
                        gOffset+13+(gripperObj*6) + 2;//L1->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 8] = 
                        x_INDEX(3,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 9] = 
                        gOffset+13+(gripperObj*6) + 3;//L2->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 9] = 
                        x_INDEX(3,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 10] = 
                        gOffset+13+(gripperObj*6) + 4;//L3->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 10] = 
                        x_INDEX(3,interval+controlPoint);

                    iRow[jacOffset + 47 + (gripperObj*21) + 11] = 
                        gOffset+13+(gripperObj*6) + 1;//static->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 11] = 
                        x_INDEX(4,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 12] = 
                        gOffset+13+(gripperObj*6) + 2;//L1->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 12] = 
                        x_INDEX(4,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 13] = 
                        gOffset+13+(gripperObj*6) + 3;//L2->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 13] = 
                        x_INDEX(4,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 14] = 
                        gOffset+13+(gripperObj*6) + 4;//L3->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 14] = 
                        x_INDEX(4,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 15] = 
                        gOffset+13+(gripperObj*6) + 5;//L4->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 15] = 
                        x_INDEX(4,interval+controlPoint);

                    iRow[jacOffset + 47 + (gripperObj*21) + 16] = 
                        gOffset+13+(gripperObj*6) + 1;//static->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 16] = 
                        x_INDEX(5,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 17] = 
                        gOffset+13+(gripperObj*6) + 2;//L1->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 17] = 
                        x_INDEX(5,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 18] = 
                        gOffset+13+(gripperObj*6) + 3;//L2->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 18] = 
                        x_INDEX(5,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 19] = 
                        gOffset+13+(gripperObj*6) + 4;//L3->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 19] = 
                        x_INDEX(5,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 20] = 
                        gOffset+13+(gripperObj*6) + 5;//L4->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 20] = 
                        x_INDEX(5,interval+controlPoint);
                    iRow[jacOffset + 47 + (gripperObj*21) + 21] = 
                        gOffset+13+(gripperObj*6) + 6;//L5->obj[gripperObj]
                    jCol[jacOffset + 47 + (gripperObj*21) + 21] = 
                        x_INDEX(5,interval+controlPoint);
                }



                //adjust jacOffset for next interval and control point
                jacOffset += 48 + (21*m_gripperObjects.size());
            }
            //adjust gOffset for next interval
            gOffset += 14 + (6*m_gripperObjects.size());
        }


        

    }
    else
    {
        std::cout << "eval_jac_g subsequent time\n";
        //calculate the values of the jacobian

        //calculate intervalDuration
        double intervalDuration = duration/m_numIntervals;

        values[0] = 1; // time constraint is linear with gradient 1

        //constraints on max/min pos/vel/acc
        //set offsets
        jacOffset = 1;
        gOffset = 1;
        int jacPosOffsetPer = (m_numIntervals * numJoints * 4);
        int jacVelAccOffsetPer = (m_numIntervals * numJoints * 5);
        for (int interval = 0; interval < m_numIntervals; interval++) { //m_numIntervals for each interval's constraints
            for (int joint = 0; joint < numJoints; joint++) { //numJoints for each joint's constraints
                //calculate differential of max/min pos/vel/acc for this interval and joint
                findIntervalControlPoints(x, n, joint, interval, controlPoints);
                differentiateMaxMinPos(controlPoints, dMaxPos, dMinPos);
                differentiateMaxMinVel(controlPoints, duration, dMaxVel, dMinVel);
                differentiateMaxMinAcc(controlPoints, duration, dMaxAcc, dMinAcc);
                //copy max/min pos derivatives to values
                for (int j = 0; j < 4; j++) { //4 variables can affect each pos constraint
                    values[(interval*numJoints + joint)*4 + j + (jacOffset  + 0*jacPosOffsetPer)] = 
                        dMaxPos[j];  
                    values[(interval*numJoints + joint)*4 + j + (jacOffset  + 1*jacPosOffsetPer)] = 
                        dMinPos[j];  
                }
                //copy max/min vel/acc derivatives to values
                for (int j = 0; j < 5; j++) { //5 variables can affect each vel/acc constraint
                    values[(interval*numJoints + joint)*5 + j + (jacOffset  + 2*jacPosOffsetPer + 0*jacVelAccOffsetPer)] = 
                        dMaxVel[j];
                    values[(interval*numJoints + joint)*5 + j + (jacOffset  + 2*jacPosOffsetPer + 1*jacVelAccOffsetPer)] = 
                        dMinVel[j];
                    values[(interval*numJoints + joint)*5 + j + (jacOffset  + 2*jacPosOffsetPer + 2*jacVelAccOffsetPer)] = 
                        dMaxAcc[j];
                    values[(interval*numJoints + joint)*5 + j + (jacOffset  + 2*jacPosOffsetPer + 3*jacVelAccOffsetPer)] = 
                        dMinAcc[j];
                }
            }
        }
        
        //set the offset
        jacOffset = 1 +
                    2*m_numIntervals*numJoints*4 +
                    4*m_numIntervals*numJoints*5;
        gOffset = 1 + 
                  6*m_numIntervals*numJoints;
        //start/end Pos/Vel constraints
        //derivatives are ordered C-3, C-2, C-1, C0, duration, with elements missing as appropriate
        for (int joint = 0; joint < numJoints; joint++) { //numJoints for each joint's startPos, startVel, endPos and endVel
            //startPos constraint derivatives
            values[jacOffset + (numJoints*3)*0 + joint*3 + 0] = 1.0/(6.0);
            values[jacOffset + (numJoints*3)*0 + joint*3 + 1] = 4.0/(6.0);
            values[jacOffset + (numJoints*3)*0 + joint*3 + 2] = 1.0/(6.0);

            //startVel
            values[jacOffset + (numJoints*3)*1 + joint*3 + 0] = -3.0/(6.0*intervalDuration);
            values[jacOffset + (numJoints*3)*1 + joint*3 + 1] = 3.0/(6.0*intervalDuration);
            values[jacOffset + (numJoints*3)*1 + joint*3 + 2] = -(  x_VALUE(joint, -1) * (3.0/6) + 
                                                                    x_VALUE(joint, -3) * (-3.0/6) ) / 
                                                                 (intervalDuration*duration);

            //endVel
            values[jacOffset + (numJoints*3)*11 + joint*3 + 0] = -3.0/(6.0*intervalDuration);
            values[jacOffset + (numJoints*3)*11 + joint*3 + 1] = 3.0/(6.0*intervalDuration);
            values[jacOffset + (numJoints*3)*11 + joint*3 + 2] = -(  x_VALUE(joint, m_numIntervals-0-1) * (3.0/6) + 
                                                                    x_VALUE(joint, m_numIntervals-2-1) * (-3.0/6) ) / 
                                                                 (intervalDuration*duration);;
        }

        //endPos
        double endPosDerivatives[9*numJoints*3];
        differentiateEndPos(x, n, endPosDerivatives);
        for (int endPosDerivative = 0; endPosDerivative < 9*numJoints*3; endPosDerivative++) {
            values[jacOffset + (numJoints*3)*2 + endPosDerivative] = 
                endPosDerivatives[endPosDerivative];
        }

        //increase the offset
        jacOffset += (6*3) + (6*3) + (9*3*6) + (6*3);
        gOffset += 6+6+9+6;
        double tempIntervalControlPoints[6][4];
        double tempSplineParams[6][5];
        double pos[6];
        std::vector<EricsTrajOpt::DenavitHartenberg> DHParams;
        //variable to hold time and the size of the time step
        double time = 0;
        double dT = 1e-3;

        //variables to temporarily hold distance and index in distances
        double tempDistance;
        int tempIndex;
        std::vector<double> forwardIntervalMinDists(48 + (21*m_gripperObjects.size()), 2e19);
        std::vector<double> backwardIntervalMinDists(48 + (21*m_gripperObjects.size()), 2e19);
        std::vector<double> intervalMinDists[2] = {forwardIntervalMinDists, backwardIntervalMinDists};
        double differences[2] = {diffDiff, -diffDiff};
        for (int interval = 0; interval < m_numIntervals; interval++) {
            for (int joint = 0; joint < 6; joint++) {
                findIntervalControlPoints(x,n,joint,interval,tempIntervalControlPoints[joint]);
            }
            //iterate for each control point to differentiate
            for (int diffControlPoint = 0; diffControlPoint < 4; diffControlPoint++) {
                //iterate for each joint to differentiate
                for (int diffJoint = 0; diffJoint < 6; diffJoint++) {
                    for (int direction = 0; direction < 2; direction++) {
                        //calculate spline parameters for this interval, joint, control point and direction
                        tempIntervalControlPoints[diffJoint][diffControlPoint] = x_VALUE(diffJoint, interval+diffControlPoint-3) + differences[direction];
                        for (int joint = 0; joint < 6; joint++) {
                            calculateIntervalSplineParams(tempIntervalControlPoints[joint], duration, tempSplineParams[joint]);
                        }

                        //run differentiation. each joint fills in a different subsection of intervalMinDists
                        time = 0;
                        while (time < tempSplineParams[0][5]) {
                            //get arm joint positions at this time
                            for (int joint = 0; joint < 6; joint++) {
                                pos[joint] = calculatePos(tempSplineParams[joint], time/tempSplineParams[0][5]);
                            }
                            //rotate arm and gripper hulls to position
                            DHParams = getDHAtJointPose(m_DHParamsAt0, pos);
                            for (int link = 0; link < 6; link++) {
                                m_armLinks[link].setHullDH(DHParams[link]);
                            }
                            for (int gripperObj = 0; gripperObj < 6; gripperObj++) {
                                m_gripperObjects[gripperObj].setHullDH(DHParams[5]);
                            }
                            switch (diffJoint) {
                                //if we are differentiating joint 1:
                                case 0:
                                    //calculate distance between static objects and arm links 2-6
                                    for (int link = 2; link < 7; link++) {
                                        tempDistance = calcStatic_2_linkCollision(link);
                                        tempIndex = 0 + link - 2;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between each gripper object and static objects plus arm links 1-5
                                    for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                                        //calculate distance between gripper object and static objects
                                        tempDistance = calcGripper_2_StaticCollision(gripperObj);
                                        tempIndex = 47 + (gripperObj*21) + 1;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    break;
                                //if we are differentiating joint 2:
                                case 1:
                                    //calculate distance between static objects and arm links 2-6
                                    for (int link = 2; link < 7; link++) {
                                        tempDistance = calcStatic_2_linkCollision(link);
                                        tempIndex = 5 + link - 2;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 1 and arm links 3-6
                                    for (int link = 3; link < 7; link++) {
                                        tempDistance = calcLink1_2_linkCollision(link);
                                        tempIndex = 10 + link - 3;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between each gripper object and static objects plus arm links 1-1
                                    for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                                        //calculate distance between gripper object and static objects
                                        tempDistance = calcGripper_2_StaticCollision(gripperObj);
                                        tempIndex = 47 + (gripperObj*21) + 2;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                        //calculate distance between gripper object and arm links 1-1
                                        for (int link = 1; link < 2; link++) {
                                            tempDistance = calcGripper_2_linkCollision(link, gripperObj);
                                            tempIndex = 47 + (gripperObj*21) + 2 + link;
                                            //update distances if necessary
                                            if (tempDistance < intervalMinDists[direction][tempIndex]){
                                                intervalMinDists[direction][tempIndex] = tempDistance;
                                            }
                                        }
                                    }
                                    break;
                                //if we are differentiating joint 3:
                                case 2:
                                    //calculate distance between static objects and arm links 3-6
                                    for (int link = 3; link < 7; link++) {
                                        tempDistance = calcStatic_2_linkCollision(link);
                                        tempIndex = 14 + link - 3;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 1 and arm links 3-6
                                    for (int link = 3; link < 7; link++) {
                                        tempDistance = calcLink1_2_linkCollision(link);
                                        tempIndex = 18 + link - 3;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 2 and arm links 4-6
                                    for (int link = 4; link < 7; link++) {
                                        tempDistance = calcLink2_2_linkCollision(link);
                                        tempIndex = 22 + link - 4;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between each gripper object and static objects plus arm links 1-2
                                    for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                                        //calculate distance between gripper object and static objects
                                        tempDistance = calcGripper_2_StaticCollision(gripperObj);
                                        tempIndex = 47 + (gripperObj*21) + 4;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                        //calculate distance between gripper object and arm links 1-2
                                        for (int link = 1; link < 3; link++) {
                                            tempDistance = calcGripper_2_linkCollision(link, gripperObj);
                                            tempIndex = 47 + (gripperObj*21) + 4 + link;
                                            //update distances if necessary
                                            if (tempDistance < intervalMinDists[direction][tempIndex]){
                                                intervalMinDists[direction][tempIndex] = tempDistance;
                                            }
                                        }
                                    }
                                    break;
                                //if we are differentiating joint 4:
                                case 3:
                                    //calculate distance between static objects and arm links 4-6
                                    for (int link = 4; link < 7; link++) {
                                        tempDistance = calcStatic_2_linkCollision(link);
                                        tempIndex = 25 + link - 4;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 1 and arm links 4-6
                                    for (int link = 4; link < 7; link++) {
                                        tempDistance = calcLink1_2_linkCollision(link);
                                        tempIndex = 28 + link - 4;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 2 and arm links 4-6
                                    for (int link = 4; link < 7; link++) {
                                        tempDistance = calcLink2_2_linkCollision(link);
                                        tempIndex = 31 + link - 4;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 3 and arm links 5-6
                                    for (int link = 5; link < 7; link++) {
                                        tempDistance = calcLink3_2_linkCollision(link);
                                        tempIndex = 34 + link - 5;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between each gripper object and static objects plus arm links 1-3
                                    for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                                        //calculate distance between gripper object and static objects
                                        tempDistance = calcGripper_2_StaticCollision(gripperObj);
                                        tempIndex = 47 + (gripperObj*21) + 7;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                        //calculate distance between gripper object and arm links 1-3
                                        for (int link = 1; link < 4; link++) {
                                            tempDistance = calcGripper_2_linkCollision(link, gripperObj);
                                            tempIndex = 47 + (gripperObj*21) + 7 + link;
                                            //update distances if necessary
                                            if (tempDistance < intervalMinDists[direction][tempIndex]){
                                                intervalMinDists[direction][tempIndex] = tempDistance;
                                            }
                                        }
                                    }
                                    break;
                                //if we are differentiating joint 5:
                                case 4:
                                    //calculate distance between static objects and arm links 5-6
                                    for (int link = 5; link < 7; link++) {
                                        tempDistance = calcStatic_2_linkCollision(link);
                                        tempIndex = 36 + link - 5;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 1 and arm links 5-6
                                    for (int link = 5; link < 7; link++) {
                                        tempDistance = calcLink1_2_linkCollision(link);
                                        tempIndex = 38 + link - 5;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 2 and arm links 5-6
                                    for (int link = 5; link < 7; link++) {
                                        tempDistance = calcLink2_2_linkCollision(link);
                                        tempIndex = 40 + link - 5;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 3 and arm links 5-6
                                    for (int link = 5; link < 7; link++) {
                                        tempDistance = calcLink3_2_linkCollision(link);
                                        tempIndex = 42 + link - 5;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between each gripper object and static objects plus arm links 1-4
                                    for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                                        //calculate distance between gripper object and static objects
                                        tempDistance = calcGripper_2_StaticCollision(gripperObj);
                                        tempIndex = 47 + (gripperObj*21) + 11;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                        //calculate distance between gripper object and arm links 1-4
                                        for (int link = 1; link < 5; link++) {
                                            tempDistance = calcGripper_2_linkCollision(link, gripperObj);
                                            tempIndex = 47 + (gripperObj*21) + 11 + link;
                                            //update distances if necessary
                                            if (tempDistance < intervalMinDists[direction][tempIndex]){
                                                intervalMinDists[direction][tempIndex] = tempDistance;
                                            }
                                        }
                                    }
                                    break;
                                //if we are differentiating joint 6:
                                case 5:
                                    //calculate distance between static objects and arm links 6-6
                                    for (int link = 6; link < 7; link++) {
                                        tempDistance = calcStatic_2_linkCollision(link);
                                        tempIndex = 36 + link - 6;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 1 and arm links 6-6
                                    for (int link = 6; link < 7; link++) {
                                        tempDistance = calcLink1_2_linkCollision(link);
                                        tempIndex = 38 + link - 6;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 2 and arm links 6-6
                                    for (int link = 6; link < 7; link++) {
                                        tempDistance = calcLink2_2_linkCollision(link);
                                        tempIndex = 40 + link - 6;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between arm link 3 and arm links 6-6
                                    for (int link = 6; link < 7; link++) {
                                        tempDistance = calcLink3_2_linkCollision(link);
                                        tempIndex = 42 + link - 6;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                    }
                                    //calculate distance between each gripper object and static objects plus arm links 1-5
                                    for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                                        //calculate distance between gripper object and static objects
                                        tempDistance = calcGripper_2_StaticCollision(gripperObj);
                                        tempIndex = 47 + (gripperObj*21) + 16;
                                        //update distances if necessary
                                        if (tempDistance < intervalMinDists[direction][tempIndex]){
                                            intervalMinDists[direction][tempIndex] = tempDistance;
                                        }
                                        //calculate distance between gripper object and arm links 1-5
                                        for (int link = 1; link < 6; link++) {
                                            tempDistance = calcGripper_2_linkCollision(link, gripperObj);
                                            tempIndex = 47 + (gripperObj*21) + 16 + link;
                                            //update distances if necessary
                                            if (tempDistance < intervalMinDists[direction][tempIndex]){
                                                intervalMinDists[direction][tempIndex] = tempDistance;
                                            }
                                        }
                                    }
                                    break;
                            }
                            time+=dT;
                        } //end of time iteration loop
                    } //end of direction iterating loop
                } //end of diffJoint loop
                
                //calculate derivatives for this interval and controlPoint and put in values array
                for (int i = 0; i < 48 + (21*(int)m_gripperObjects.size()); i++) {
                    values[jacOffset + i] = (intervalMinDists[0][i] - intervalMinDists[1][i]) / (2*diffDiff);
                }
                //adjust jacOffset for next interval and control point
                jacOffset += 48 + (21*m_gripperObjects.size());
            } //end of diffCOntrolPoint loop
            //adjust gOffset for next interval
            gOffset += 14 + (6*m_gripperObjects.size());
        }
    }
    std::cout << "eval_jac_g done\n";
    return true;
}
// [TNLP_eval_jac_g]

// [TNLP_eval_h]
//return the structure or values of the Hessian
bool TrajectoryOptimisationNLP::eval_h(
    Index         n,
    const Number* x,
    bool          new_x,
    Number        obj_factor,
    Index         m,
    const Number* lambda,
    bool          new_lambda,
    Index         nele_hess,
    Index*        iRow,
    Index*        jCol,
    Number*       values)
{
    if( values == NULL )
    {
        //std::cout << "eval_h first time\n";
        
        //a table of which joints have common constraints, and therefore cause nonzero hessian entries
        /*bool jointInteractions[6][6] = {{true , false, true , true , true , true },
                                        {false, true , false, true , true , true },
                                        {true , false, true , false, true , true },
                                        {true , true , false, true , false, true },
                                        {true , true , true , false, true , true },
                                        {true , true , true , true , true , true }};*/
        bool jointInteractions[6][6] = {{true , false, false, false, false, false},
                                        {false, true , false, false, false, false},
                                        {false, false, true , false, false, false},
                                        {false, false, false, true , false, false},
                                        {false, false, false, false, true , false},
                                        {false, false, false, false, false, true }};
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.

        //
        Index idx = 0;
        //iterate through each column of h, (except duration), keeping control point and joint number separate
        for (Index controlPoint_j = 0; controlPoint_j < m_numControlPoints; controlPoint_j++) {
            for (Index joint_j = 0; joint_j < numJoints; joint_j++) {
                //iterate through each row of h, (except duration), keeping control point and joint number separate
                for (Index controlPoint_i = controlPoint_j; controlPoint_i < m_numControlPoints; controlPoint_i++) {
                    for (Index joint_i = 0; joint_i < numJoints; joint_i++) {
                        //we only need lower left triangle of h, row must be >= column
                        if (x_INDEX(joint_i, controlPoint_i-3) >= x_INDEX(joint_j, controlPoint_j-3)) {
                            //row control point must be used by interval with column control point as C-3
                            if (controlPoint_i - controlPoint_j < 4) {
                                //check if pair of joint has an interaction, as defined in jointInteractions
                                if (jointInteractions[joint_i][joint_j]) {
                                    hessianIndexMap [std::make_pair(x_INDEX(joint_j, controlPoint_j-3), 
                                                                    x_INDEX(joint_i, controlPoint_i-3))
                                                    ] = idx;
                                    jCol[idx] = x_INDEX(joint_j, controlPoint_j-3);
                                    iRow[idx] = x_INDEX(joint_i, controlPoint_i-3);
                                    /*std::cout << idx << 
                                                " | C" << joint_j << "_" << controlPoint_j-3 << 
                                                " | C" << joint_i << "_" << controlPoint_i-3 << "\n";*/
                                    idx++;
                                    
                                }
                            }
                        }
                    }
                }
                // each column also needs an entry for duration
                hessianIndexMap [std::make_pair(x_INDEX(joint_j, controlPoint_j-3), 
                                                n-1) // index of duration
                                ] = idx;
                jCol[idx] = x_INDEX(joint_j, controlPoint_j-3);
                iRow[idx] = n-1;
                /*std::cout << idx << 
                            " | C" << joint_j << "_" << controlPoint_j-3 << 
                            " | duration" << "\n";*/
                idx++;
            }
        }
        // an entry for 2nd-order duration derivative is required
        hessianIndexMap [std::make_pair(n-1, n-1) // index of duration
                        ] = idx;
        jCol[idx] = n-1;
        iRow[idx] = n-1;
        /*std::cout << idx << 
            " | duration" << 
            " | duration" << "\n";*/
        idx++;
    }
    else {
        //std::cout << "eval_h subsequent time\n";
        // return the values. This is a symmetric matrix, fill the lower left
        // triangle only

        //variables to hold sub-hessians of min/max pos/vel/acc
        double  dMaxPos[4][4], dMinPos[4][4], 
                dMaxVel[5][5], dMinVel[5][5], 
                dMaxAcc[5][5], dMinAcc[5][5];
        //variables to hold the lambda of min/max pos/vel/acc for the current joint and interval
        double  lambdaMaxPos, lambdaMinPos, 
                lambdaMaxVel, lambdaMinVel, 
                lambdaMaxAcc, lambdaMinAcc;
        //variables to hold the lambda of start/end vel for the current joint
        double  lambdaStartVel, lambdaEndVel;
        
        //variable to hold control points of spline of interval and joint
        double controlPoints[4];

        //calculate intervalDuration
        double intervalDuration = duration/m_numIntervals;
    
        // fill the objective portion
        // as the objective is linear, all elements of the hessian will be zero. use this to initialise values.
        for (Index idx = 0; idx < nele_hess; idx++) {
            values[idx] = 0;
        }

        // fill the constraints portion

        // duration constraint does not affect hessian, ignore

        // aggregate values from max/min pos/vel/acc constraints
        for (Index interval = 0; interval < m_numIntervals; interval++) {
            for (Index joint = 0; joint < numJoints; joint++) {
                //find lambda for ech constraint for this joint and interval
                lambdaMaxPos = lambda[(m_numIntervals * numJoints * 0) + (interval*numJoints + joint) + 1];
                lambdaMinPos = lambda[(m_numIntervals * numJoints * 1) + (interval*numJoints + joint) + 1];
                lambdaMaxVel = lambda[(m_numIntervals * numJoints * 2) + (interval*numJoints + joint) + 1];
                lambdaMinVel = lambda[(m_numIntervals * numJoints * 3) + (interval*numJoints + joint) + 1];
                lambdaMaxAcc = lambda[(m_numIntervals * numJoints * 4) + (interval*numJoints + joint) + 1];
                lambdaMinAcc = lambda[(m_numIntervals * numJoints * 5) + (interval*numJoints + joint) + 1];
                //calculate sub-hessians of max/min pos/vel/acc for interval and joint
                findIntervalControlPoints(x, n, joint, interval, controlPoints);
                hessianMaxMinPos(controlPoints, dMaxPos, dMinPos);
                hessianMaxMinVel(controlPoints, duration, dMaxVel, dMinVel);
                hessianMaxMinAcc(controlPoints, duration, dMaxAcc, dMinAcc);
                //aggregate sub-hessians to correct place in h/values
                for (Index col = 0; col < 4; col++) {
                    for (Index row = col; row < 4; row++) {
                        //for controlPoint-controlPoint derivatives
                        values  [hessianIndexMap[std::make_pair(x_INDEX(joint,interval-3+col),
                                                                x_INDEX(joint,interval-3+row)
                                                               )
                                                ]
                                ] += (  (lambdaMaxPos*dMaxPos[col][row]) + (lambdaMinPos*dMinPos[col][row]) +
                                        (lambdaMaxVel*dMaxVel[col][row]) + (lambdaMinVel*dMinVel[col][row]) +
                                        (lambdaMaxAcc*dMaxAcc[col][row]) + (lambdaMinAcc*dMinAcc[col][row]));
                        /*std::cout << "index: " << hessianIndexMap[std::make_pair(x_INDEX(joint,interval-3+col),
                                                                    x_INDEX(joint,interval-3+row)
                                                                   )] << " | " <<
                                    x_INDEX(joint,interval-3+col) << " | " <<
                                    x_INDEX(joint,interval-3+row) << " || " << 
                                    joint << " " << interval << " " << col << " " << row << "\n";*/
                    }
                    //for controlPoint-duration derivatives (vel and acc only)
                    values[hessianIndexMap[std::make_pair(x_INDEX(joint,interval-3+col), n-1)]] 
                            += ((lambdaMaxVel*dMaxVel[col][4]) + (lambdaMinVel*dMinVel[col][4]) +
                                (lambdaMaxAcc*dMaxAcc[col][4]) + (lambdaMinAcc*dMinAcc[col][4]));
                    /*std::cout << hessianIndexMap[std::make_pair(x_INDEX(joint,interval-3+col),
                                                                n-1
                                                               )] << " | " <<
                                x_INDEX(joint,interval-3+col) << " | " <<
                                n-1 << "\n";*/
                }
                //for duration-duration derivative (vel and acc only)
                values[hessianIndexMap[std::make_pair(n-1, n-1)]] 
                        += ((lambdaMaxVel*dMaxVel[4][4]) + (lambdaMinVel*dMinVel[4][4]) +
                            (lambdaMaxAcc*dMaxAcc[4][4]) + (lambdaMinAcc*dMinAcc[4][4]));
                /*std::cout << hessianIndexMap[std::make_pair(n-1, n-1)] << " | " <<
                                n-1<< " | " << n-1 << "\n";*/
            }
        }

        //start and end pos constraints have no nonzero second order derivatives, ignore

        //start and end vel constraints have 3 nonzero second order derivatives per joint each
        //these are the derivatives calculated in the jacobian of g, further differentiated by duration
        for (Index joint = 0; joint < numJoints; joint++) {
            lambdaStartVel = lambda[joint + 6 + (6*m_numIntervals*numJoints) + 1];
            lambdaEndVel = lambda[joint + 18 + (6*m_numIntervals*numJoints) + 1];
            //aggregate startVel S.O.Ds
            //d^2/(dC-1 dduration)
            values[hessianIndexMap[std::make_pair(x_INDEX(joint,-1), n-1)]] 
                +=  lambdaStartVel * 
                    (-3.0/(6.0*intervalDuration*duration));
            //d^2/(dC-3 dduration)
            values[hessianIndexMap[std::make_pair(x_INDEX(joint,-3), n-1)]] 
                +=  lambdaStartVel * 
                    (3.0/(6.0*intervalDuration*duration));
            //d^2/(dduration^2)
            values[hessianIndexMap[std::make_pair(n-1, n-1)]] 
                +=  lambdaStartVel * 
                    ((x_VALUE(joint,-3) - x_VALUE(joint,-1)) / 
                    (intervalDuration*duration*duration));

            //aggregate endVel S.O.Ds
            //d^2/(dC0 dduration)
            values[hessianIndexMap[std::make_pair(x_INDEX(joint, (m_numIntervals-1)-0), n-1)]] 
                +=  lambdaEndVel * 
                    (-3.0/(6.0*intervalDuration*duration));
            //d^2/(dC-2 dduration)
            values[hessianIndexMap[std::make_pair(x_INDEX(joint,(m_numIntervals-1)-2), n-1)]] 
                +=  lambdaEndVel * 
                    (3.0/(6.0*intervalDuration*duration));
            //d^2/(dduration^2)
            values[hessianIndexMap[std::make_pair(n-1, n-1)]] 
                +=  lambdaEndVel * 
                    ((x_VALUE(joint,(m_numIntervals-1)-2) - x_VALUE(joint,(m_numIntervals-1)-0)) / 
                    (intervalDuration*duration*duration));
            /*std::cout << 
            hessianIndexMap[std::make_pair(x_INDEX(joint,-1), n-1)] << " | " << 
            hessianIndexMap[std::make_pair(x_INDEX(joint,-3), n-1)] << " | " << 
            hessianIndexMap[std::make_pair(x_INDEX(joint, (m_numIntervals-1)-0), n-1)] << " | " << 
            hessianIndexMap[std::make_pair(x_INDEX(joint,(m_numIntervals-1)-2), n-1)] << " | " << 
            hessianIndexMap[std::make_pair(n-1, n-1)] << "\n";*/

        }
        //double breakTHings = 1/0;
    }
    //std::cout << "eval_h done\n";
    return true;
}
// [TNLP_eval_h]

// [TNLP_finalize_solution]
void TrajectoryOptimisationNLP::finalize_solution(
    SolverReturn               status,
    Index                      n,
    const Number*              x,
    const Number*              z_L,
    const Number*              z_U,
    Index                      m,
    const Number*              g,
    const Number*              lambda,
    Number                     obj_value,
    const IpoptData*           ip_data,
    IpoptCalculatedQuantities* ip_cq
)
{
    /* A variety of actions to take upon completion of the NLP solver. frequently changed during development,
        hence the commented sections */
 
    /*std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
    for( Index i = 0; i < n; i++ )
    {
       std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }*/
 
    /*std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
    for( Index i = 0; i < n; i++ )
    {
       std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
    }
    for( Index i = 0; i < n; i++ )
    {
       std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
    }*/
 
    /*std::cout << std::endl << std::endl << "Objective value" << std::endl;
    std::cout << "f(x*) = " << obj_value << std::endl;*/
 
    std::cout << std::endl << "Final value of the constraints:" << std::endl;
    for( Index i = 0; i < m; i++ )
    {
       std::cout << "g(" << i << ") = " << g[i] << std::endl;
    }

    // write final values of x to the returnX vector
    double splineCoeffs[5];
    double endPose[6];
    for (int joint = 0; joint < numJoints; joint++) {
        m_returnX.push_back(*new std::vector<double>);
        for (int controlPoint = -3; controlPoint < m_numIntervals; controlPoint++) {
            m_returnX[joint].push_back(x_VALUE(joint, controlPoint));
        }
        //also calculate the end pose of the arm so that a model of the arm at the final position can be generated
        findIntervalCoefficients(x, n, joint, m_numIntervals - 1, splineCoeffs);
        std::cout << calculatePos(splineCoeffs, 1) << "    ";
        endPose[joint] = calculatePos(splineCoeffs, 1);
    }
    std::cout << "\n";
    m_returnX.push_back(*new std::vector<double>);
    m_returnX[numJoints].push_back(duration);


    /* create a model of the UR5 arm in the final position of the trajectory and export */
    std::string filePath = "/home/eric/Documents/MENG-Robot-Fab-Project/Optimisation Code/Models/ArmSegments/convexHulls/";
    std::string fileName1 = "UR5_link_1_CH.stl";
    std::string fileName2 = "UR5_link_2_CH.stl";
    std::string fileName3 = "UR5_link_3_CH.stl";
    std::string fileName4 = "UR5_link_4_CH.stl";
    std::string fileName5 = "UR5_link_5_CH.stl";
    std::string fileName6 = "UR5_link_6_CH.stl";
    std::string outFilename = "UR5_whole_arm.stl";
    EricsTrajOpt::Hull hull1 = EricsTrajOpt::readBinaryStl(filePath+fileName1);
    EricsTrajOpt::Hull hull2 = EricsTrajOpt::readBinaryStl(filePath+fileName2);
    EricsTrajOpt::Hull hull3 = EricsTrajOpt::readBinaryStl(filePath+fileName3);
    EricsTrajOpt::Hull hull4 = EricsTrajOpt::readBinaryStl(filePath+fileName4);
    EricsTrajOpt::Hull hull5 = EricsTrajOpt::readBinaryStl(filePath+fileName5);
    EricsTrajOpt::Hull hull6 = EricsTrajOpt::readBinaryStl(filePath+fileName6);
    EricsTrajOpt::Hull wholeArm;

    std::vector<EricsTrajOpt::DenavitHartenberg> endDH = getDHAtJointPose(m_DHParamsAt0, endPose);
    
    hull1.setHullDH(endDH[0]);
    hull2.setHullDH(endDH[1]);
    hull3.setHullDH(endDH[2]);
    hull4.setHullDH(endDH[3]);
    hull5.setHullDH(endDH[4]);
    hull6.setHullDH(endDH[5]);
    wholeArm = hull1 + hull2 + hull3 + hull4 + hull5 + hull6;
    wholeArm.exportBinaryStl(filePath+outFilename);

    /* calculate the end positon of the arm in PosNorm form for verification purposes. */
    std::vector<double> posNorm = DH2PosNorm(endDH[5]);
    for(int i = 0; i < 9; i++) {
        std::cout << posNorm[i] << "\n";
    }

    /* display a graph of the final trajectory. */
    graphSolution(x, n);
}
// [TNLP_finalize_solution]

// function to find the minimum and maximum pos of a spline within an interval,
// given spline coefficients (A-D) 
void TrajectoryOptimisationNLP::findMaxMinPos(
    const double    splineParameters[5],
    double *        maxPos, 
    double *        minPos) {
    
    //temporary double to hold values we calculate
    double tempValue;



    // initial value of max/minPos - either value at tau=0 is > value at tau=1
    // or value at tau=0 <= value at tau=1. Set initial values appropriately.
    if (calculatePos(splineParameters, 0) > calculatePos(splineParameters, 1)) {
        *maxPos = calculatePos(splineParameters, 0);
        *minPos = calculatePos(splineParameters, 1);
    } else {
        *maxPos = calculatePos(splineParameters, 1);
        *minPos = calculatePos(splineParameters, 0);
    }
    //check for real roots of derivative of position using discriminant
    if (4*splineParameters[B]*splineParameters[B] > 12*splineParameters[A]*splineParameters[C]) {
        //if stationary points of pos are real, calculate them
        double stationaryPosPoints[2] = { (-2*splineParameters[B] + 
                                           std::sqrt(4*splineParameters[B]*splineParameters[B] - 
                                                     12*splineParameters[A]*splineParameters[C])) / 
                                          (6*splineParameters[A]) ,
                                          (-2*splineParameters[B] - 
                                           std::sqrt(4*splineParameters[B]*splineParameters[B] - 
                                                     12*splineParameters[A]*splineParameters[C])) /
                                          (6*splineParameters[A]) };
        for (int i = 0; i < 2; i++) {
            if (stationaryPosPoints[i] >= 0 and
                stationaryPosPoints[i] <= 1) {
                //if stationary point is within the domain, calculate pos value and store in tempValue
                tempValue = calculatePos(splineParameters, stationaryPosPoints[i]);
                //compare to current min/maxPos and update if necessary.
                if (tempValue > *maxPos) {
                    *maxPos = tempValue;
                } else if (tempValue < *minPos) {
                    *minPos = tempValue;
                }
            }
        }
    }
}

// function to find the minimum and maximum vel of a spline within an interval,
// given spline coefficients (A-C) and segment duration (deltaT)
void TrajectoryOptimisationNLP::findMaxMinVel(
    const double    splineParameters[5],
    double *        maxVel, 
    double *        minVel) {

    //temporary double to hold values we calculate
    double tempValue;

    // initial value of max/minVel - either value at tau=0 is > value at tau=1
    // or value at tau=0 <= value at tau=1. Set initial values appropriately.
    if (calculateVel(splineParameters, 0) > calculateVel(splineParameters, 1)) {
        *maxVel = calculateVel(splineParameters, 0);
        *minVel = calculateVel(splineParameters, 1);
    } else {
        *maxVel = calculateVel(splineParameters, 1);
        *minVel = calculateVel(splineParameters, 0);
    }
    //check for stationary point of velocity
    if (splineParameters[A] != 0) {
        //calculate stationary velocity point
        double stationaryVelPoint = -splineParameters[B]/(3*splineParameters[A]);
        if (stationaryVelPoint >= 0 and
            stationaryVelPoint <= 1) {
            //if stationary point is within the domain, calculate vel value and store in tempValue
            tempValue = calculateVel(splineParameters, stationaryVelPoint);
            //compare to current min/maxVel and update if necessary.
            if (tempValue > *maxVel) {
                *maxVel = tempValue;
            } else if (tempValue < *minVel) {
                *minVel = tempValue;
            }
        }
    }
}

// function to find the minimum and maximum acc of a spline within an interval,
// given spline coefficients (A-B) and segment duration (deltaT)
void TrajectoryOptimisationNLP::findMaxMinAcc(
    const double    splineParameters[5],
    double *        maxAcc, 
    double *        minAcc) {

    // initial value of max/minAcc - either value at tau=0 is > value at tau=1
    // or value at tau=0 <= value at tau=1. Set initial values appropriately.
    if (calculateAcc(splineParameters, 0) > calculateAcc(splineParameters, 1)) {
        *maxAcc = calculateAcc(splineParameters, 0);
        *minAcc = calculateAcc(splineParameters, 1);
    } else {
        *maxAcc = calculateAcc(splineParameters, 1);
        *minAcc = calculateAcc(splineParameters, 0);
    }
    // acc has no stationary points that need considering. case where acc is constant rather
    // than linear is accounted for when value at tau=0 == value at tau=1.
}

// function to differentiate the minimum and maximum pos of a spline within an interval,
// given spline control points C1-C4 
void TrajectoryOptimisationNLP::differentiateMaxMinPos(
    const double    controlPoints [4],
    double          dMaxPos[4], 
    double          dMinPos[4]) {
    
    //variable to store a master copy of differentiation variables, in the correct order
    double  masterControlPoints[4] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]};
    //variable to store a temporary copy of differentiation variables
    double  tempControlPoints[4] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]};

    //variables to temporarily store forward and backward components of each differentiation approximation
    double  forwardSplineParams[5], backwardSplineParams[5],
            forwardMaxPos, forwardMinPos, backwardMaxPos, backwardMinPos;
    
    for(int i = 0; i < 4; i++) {
        //calculate spline parameters at forward and backward points
        tempControlPoints[i] = masterControlPoints[i] + diffDiff;
        calculateIntervalSplineParams(tempControlPoints, 0, forwardSplineParams);
        tempControlPoints[i] = masterControlPoints[i] - diffDiff;
        calculateIntervalSplineParams(tempControlPoints, 0, backwardSplineParams);
        //reset tempControlPoints
        tempControlPoints[i] = masterControlPoints[i];

        //find max and min pos at forward and backward points
        findMaxMinPos(forwardSplineParams, &forwardMaxPos, &forwardMinPos);
        findMaxMinPos(backwardSplineParams, &backwardMaxPos, &backwardMinPos);

        //calculate approximate derivative W.R.T controlPoints[i]
        dMaxPos[i] = (forwardMaxPos-backwardMaxPos)/(2*diffDiff);
        dMinPos[i] = (forwardMinPos-backwardMinPos)/(2*diffDiff);
    }
}
    
// function to differentiate the minimum and maximum vel of a spline within an interval,
// given spline control points C1-C4 and segment duration (deltaT)
void TrajectoryOptimisationNLP::differentiateMaxMinVel(
    const double    controlPoints [4],
    const double    totalDuration,
    double          dMaxVel[5], 
    double          dMinVel[5]) {

    //variable to store a master copy of controlPoints and totalDuration
    double  masterControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};
    //variable to store a temporary copy of controlPoints and totalDuration
    double  tempControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};

    //variables to temporarily store forward and backward components of each differentiation approximation
    double  forwardSplineParams[5], backwardSplineParams[5],
            forwardMaxVel, forwardMinVel, backwardMaxVel, backwardMinVel;
    
    for(int i = 0; i < 5; i++) {
        //calculate spline parameters at forward and backward points
        tempControlPoints[i] = masterControlPoints[i] + diffDiff;
        calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], forwardSplineParams);
        tempControlPoints[i] = masterControlPoints[i] - diffDiff;
        calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], backwardSplineParams);
        //reset tempControlPoints
        tempControlPoints[i] = masterControlPoints[i];

        //find max and min pos at forward and backward points
        findMaxMinVel(forwardSplineParams, &forwardMaxVel, &forwardMinVel);
        findMaxMinVel(backwardSplineParams, &backwardMaxVel, &backwardMinVel);

        //calculate approximate derivative W.R.T controlPoints[i]
        dMaxVel[i] = (forwardMaxVel-backwardMaxVel)/(2*diffDiff);
        dMinVel[i] = (forwardMinVel-backwardMinVel)/(2*diffDiff);
    }
}

// function to differentiate the minimum and maximum acc of a spline within an interval,
// given spline control points C1-C4 and segment duration (deltaT)
void TrajectoryOptimisationNLP::differentiateMaxMinAcc(
    const double    controlPoints [4],
    const double    totalDuration, 
    double          dMaxAcc[5], 
    double          dMinAcc[5]) {

    //variable to store a master copy of controlPoints and totalDuration
    double  masterControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};
    //variable to store a temporary copy of controlPoints and totalDuration
    double  tempControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};

    //variables to temporarily store forward and backward components of each differentiation approximation
    double  forwardSplineParams[5], backwardSplineParams[5],
            forwardMaxAcc, forwardMinAcc, backwardMaxAcc, backwardMinAcc;
    
    for(int i = 0; i < 5; i++) {
        //calculate spline parameters at forward and backward points
        tempControlPoints[i] = masterControlPoints[i] + diffDiff;
        calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], forwardSplineParams);
        tempControlPoints[i] = masterControlPoints[i] - diffDiff;
        calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], backwardSplineParams);
        //reset tempControlPoints
        tempControlPoints[i] = masterControlPoints[i];

        //find max and min pos at forward and backward points
        findMaxMinAcc(forwardSplineParams, &forwardMaxAcc, &forwardMinAcc);
        findMaxMinAcc(backwardSplineParams, &backwardMaxAcc, &backwardMinAcc);

        //calculate approximate derivative W.R.T controlPoints[i]
        dMaxAcc[i] = (forwardMaxAcc-backwardMaxAcc)/(2*diffDiff);
        dMinAcc[i] = (forwardMinAcc-backwardMinAcc)/(2*diffDiff);
    }
}

/*  function to find the hessian of MaxMinPos for a particular joint and interval
    only the lower left of the return matrices are filled
    order is:
    [dC-3dC-3, dC-2dC-3, dC-1dC-3, dC0dC-3]
    [dC-3dC-2, dC-2dC-2, dC-1dC-2, dC0dC-2]
    [dC-3dC-1, dC-2dC-1, dC-1dC-1, dC0dC-1]
    [dC-3dC0 , dC-2dC0 , dC-1dC0 , dC0dC0 ] */
void TrajectoryOptimisationNLP::hessianMaxMinPos (
    const double    controlPoints [4],
    double          dMaxPos[4][4], 
    double          dMinPos[4][4]) {
    
    //variable to store a master copy of controlPoints and totalDuration
    //order is reversed from because of difference in order between controlPoints and dMin/MaxPos
    double  masterControlPoints[4] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]};
    //variable to store a temporary copy of controlPoints
    double  tempControlPoints[4] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]};

    /*  variables to temporarily store forward and backward components of each differentiation approximation
        index 0 means +diffDiff, index 1 means -diffDiff.
        first index is j(col), second index is i(row)
    */
    double  splineParams[2][2][5], maxPos[2][2], minPos[2][2];
    
    for(int j = 0; j < 4; j++) { //iterate through columns of hessian
        for(int i = j; i < 4; i++) { //start at diagonal, iterate through rows
            //calculate spline parameters at forward and backward points
            tempControlPoints[j] = masterControlPoints[j] + diffDiff;
            tempControlPoints[i] = masterControlPoints[i] + diffDiff;
            calculateIntervalSplineParams(tempControlPoints, 0, splineParams[0][0]);
            tempControlPoints[j] = masterControlPoints[j] - diffDiff;
            tempControlPoints[i] = masterControlPoints[i] + diffDiff;
            calculateIntervalSplineParams(tempControlPoints, 0, splineParams[1][0]);
            tempControlPoints[j] = masterControlPoints[j] + diffDiff;
            tempControlPoints[i] = masterControlPoints[i] - diffDiff;
            calculateIntervalSplineParams(tempControlPoints, 0, splineParams[0][1]);
            tempControlPoints[j] = masterControlPoints[j] - diffDiff;
            tempControlPoints[i] = masterControlPoints[i] - diffDiff;
            calculateIntervalSplineParams(tempControlPoints, 0, splineParams[1][1]);
            //reset tempControlPoints
            tempControlPoints[j] = masterControlPoints[j];
            tempControlPoints[i] = masterControlPoints[i];

            //find max and min pos at forward and backward points
            findMaxMinPos(splineParams[0][0], &maxPos[0][0], &minPos[0][0]);
            findMaxMinPos(splineParams[1][0], &maxPos[1][0], &minPos[1][0]);
            findMaxMinPos(splineParams[0][1], &maxPos[0][1], &minPos[0][1]);
            findMaxMinPos(splineParams[1][1], &maxPos[1][1], &minPos[1][1]);

            //calculate approximate derivative W.R.T controlPoints[j], controlPoints[i]
            dMaxPos[j][i] = (maxPos[0][0]-maxPos[1][0]-maxPos[0][1]+maxPos[1][1])/(4*diffDiff*diffDiff);
            dMinPos[j][i] = (minPos[0][0]-minPos[1][0]-minPos[0][1]+minPos[1][1])/(4*diffDiff*diffDiff);
        }
    }
}

/*  function to find the hessian of MaxMinVel for a particular joint and interval
    only the lower left of the return matrices are filled
    order is:
    [dC-3dC-3, dC-2dC-3, dC-1dC-3, dC0dC-3, dCTdC-3]
    [dC-3dC-2, dC-2dC-2, dC-1dC-2, dC0dC-2, dCTdC-2]
    [dC-3dC-1, dC-2dC-1, dC-1dC-1, dC0dC-1, dCTdC-1]
    [dC-3dC0 , dC-2dC0 , dC-1dC0 , dC0dC0 , dCTdC0 ]
    [dC-3dT  , dC-2dT  , dC-1dT  , dC0dT  , dCTdT  ] */
void TrajectoryOptimisationNLP::hessianMaxMinVel (
    const double    controlPoints [4],
    const double    totalDuration, 
    double          dMaxVel[5][5], 
    double          dMinVel[5][5]) {

    //variable to store a master copy of controlPoints and totalDuration
    double  masterControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};
    //variable to store a temporary copy of controlPoints and totalDuration
    double  tempControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};

    /*  variables to temporarily store forward and backward components of each differentiation approximation
        index 0 means +diffDiff, index 1 means -diffDiff.
        first index is j(col), second index is i(row)
    */
    double  splineParams[2][2][5], maxVel[2][2], minVel[2][2];
    
    for(int j = 0; j < 5; j++) { //iterate through columns of hessian
        for(int i = j; i < 5; i++) { //start at diagonal, iterate through rows
            //calculate spline parameters at forward and backward points
            tempControlPoints[j] = masterControlPoints[j] + diffDiff;
            tempControlPoints[i] = masterControlPoints[i] + diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[0][0]);
            tempControlPoints[j] = masterControlPoints[j] - diffDiff;
            tempControlPoints[i] = masterControlPoints[i] + diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[1][0]);
            tempControlPoints[j] = masterControlPoints[j] + diffDiff;
            tempControlPoints[i] = masterControlPoints[i] - diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[0][1]);
            tempControlPoints[j] = masterControlPoints[j] - diffDiff;
            tempControlPoints[i] = masterControlPoints[i] - diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[1][1]);
            //reset tempControlPoints
            tempControlPoints[j] = masterControlPoints[j];
            tempControlPoints[i] = masterControlPoints[i];

            //find max and min vel at forward and backward points
            findMaxMinVel(splineParams[0][0], &maxVel[0][0], &minVel[0][0]);
            findMaxMinVel(splineParams[1][0], &maxVel[1][0], &minVel[1][0]);
            findMaxMinVel(splineParams[0][1], &maxVel[0][1], &minVel[0][1]);
            findMaxMinVel(splineParams[1][1], &maxVel[1][1], &minVel[1][1]);

            //calculate approximate derivative W.R.T controlPoints[j], controlPoints[i]
            dMaxVel[j][i] = (maxVel[0][0]-maxVel[1][0]-maxVel[0][1]+maxVel[1][1])/(4*diffDiff*diffDiff);
            dMinVel[j][i] = (minVel[0][0]-minVel[1][0]-minVel[0][1]+minVel[1][1])/(4*diffDiff*diffDiff);
        }
    }
}

/*  function to find the hessian of MaxMinAcc for a particular joint and interval
    only the lower left of the return matrices are filled
    order is:
    [dC-3dC-3, dC-2dC-3, dC-1dC-3, dC0dC-3, dCTdC-3]
    [dC-3dC-2, dC-2dC-2, dC-1dC-2, dC0dC-2, dCTdC-2]
    [dC-3dC-1, dC-2dC-1, dC-1dC-1, dC0dC-1, dCTdC-1]
    [dC-3dC0 , dC-2dC0 , dC-1dC0 , dC0dC0 , dCTdC0 ]
    [dC-3dT  , dC-2dT  , dC-1dT  , dC0dT  , dCTdT  ] */
void TrajectoryOptimisationNLP::hessianMaxMinAcc (
    const double    controlPoints [4],
    const double    totalDuration, 
    double          dMaxAcc[5][5], 
    double          dMinAcc[5][5]) {
    
    //variable to store a master copy of controlPoints and totalDuration
    double  masterControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};
    //variable to store a temporary copy of controlPoints and totalDuration
    double  tempControlPoints[5] = {controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], totalDuration};

    /*  variables to temporarily store forward and backward components of each differentiation approximation
        index 0 means +diffDiff, index 1 means -diffDiff.
        first index is j(col), second index is i(row)
    */
    double  splineParams[2][2][5], maxAcc[2][2], minAcc[2][2];
    
    for(int j = 0; j < 5; j++) { //iterate through columns of hessian
        for(int i = j; i < 5; i++) { //start at diagonal, iterate through rows
            //calculate spline parameters at forward and backward points
            tempControlPoints[j] = masterControlPoints[j] + diffDiff;
            tempControlPoints[i] = masterControlPoints[i] + diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[0][0]);
            tempControlPoints[j] = masterControlPoints[j] - diffDiff;
            tempControlPoints[i] = masterControlPoints[i] + diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[1][0]);
            tempControlPoints[j] = masterControlPoints[j] + diffDiff;
            tempControlPoints[i] = masterControlPoints[i] - diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[0][1]);
            tempControlPoints[j] = masterControlPoints[j] - diffDiff;
            tempControlPoints[i] = masterControlPoints[i] - diffDiff;
            calculateIntervalSplineParams(tempControlPoints, tempControlPoints[4], splineParams[1][1]);
            //reset tempControlPoints
            tempControlPoints[j] = masterControlPoints[j];
            tempControlPoints[i] = masterControlPoints[i];

            //find max and min pos at forward and backward points
            findMaxMinAcc(splineParams[0][0], &maxAcc[0][0], &minAcc[0][0]);
            findMaxMinAcc(splineParams[1][0], &maxAcc[1][0], &minAcc[1][0]);
            findMaxMinAcc(splineParams[0][1], &maxAcc[0][1], &minAcc[0][1]);
            findMaxMinAcc(splineParams[1][1], &maxAcc[1][1], &minAcc[1][1]);

            //calculate approximate derivative W.R.T controlPoints[j], controlPoints[i]
            dMaxAcc[j][i] = (maxAcc[0][0]-maxAcc[1][0]-maxAcc[0][1]+maxAcc[1][1])/(4*diffDiff*diffDiff);
            dMinAcc[j][i] = (minAcc[0][0]-minAcc[1][0]-minAcc[0][1]+minAcc[1][1])/(4*diffDiff*diffDiff);
        }
    }
}

void TrajectoryOptimisationNLP::differentiateEndPos(
    const Number *  x, 
    Index           n,
    double *        endPosDerivatives) {

    double masterControlPoints[6][4];
    double tempControlPoints[6][4];
    double tempSplineParams[6][5];
    double tempEndPos[6] = {0,0,0,0,0,0};
    std::vector<EricsTrajOpt::DenavitHartenberg> tempEndDH;
    std::vector<double> forwardEndPosNorm;
    std::vector<double> backwardEndPosNorm;
    for (int joint = 0; joint < numJoints; joint++) {
        findIntervalControlPoints(x, n, joint, m_numIntervals-1, masterControlPoints[joint]);
        findIntervalControlPoints(x, n, joint, m_numIntervals-1, tempControlPoints[joint]);
    }
    for (int joint = 0; joint < numJoints; joint++) { //iterate for each joint
        for (int controlPoint = 1; controlPoint < 4; controlPoint++) { //iterate for each of the last 3 control points
            //calculate forward component of end PosNorm for current joint and control point 
            tempControlPoints[joint][controlPoint] = masterControlPoints[joint][controlPoint] + diffDiff;
            for(int i = 0; i < numJoints; i++) {
                calculateIntervalSplineParams(tempControlPoints[i], duration, tempSplineParams[i]);
                tempEndPos[i] = calculatePos(tempSplineParams[i], 1);
            }
            tempEndDH = getDHAtJointPose(m_DHParamsAt0, tempEndPos);
            forwardEndPosNorm = DH2PosNorm(tempEndDH[5]);

            //calculate backward component of end PosNorm for current joint and control point 
            tempControlPoints[joint][controlPoint] = masterControlPoints[joint][controlPoint] - diffDiff;
            for(int i = 0; i < numJoints; i++) {
                calculateIntervalSplineParams(tempControlPoints[i], duration, tempSplineParams[i]);
                tempEndPos[i] = calculatePos(tempSplineParams[i], 1);
            }
            tempEndDH = getDHAtJointPose(m_DHParamsAt0, tempEndPos);
            backwardEndPosNorm = DH2PosNorm(tempEndDH[5]);

            //reset tempControlPoints
            tempControlPoints[joint][controlPoint] = masterControlPoints[joint][controlPoint];

            //calculate derivatives and place in correct element of endPosDerivatives
            for(int constraint = 0; constraint < 9; constraint++) {
                endPosDerivatives[(constraint*numJoints*3) + (joint*3) + controlPoint-1] = 
                    (forwardEndPosNorm[constraint]-backwardEndPosNorm[constraint])/(2*diffDiff);
            }
        }
    }
}

void TrajectoryOptimisationNLP::findIntervalCoefficients(
    const Number *  x,  
    Index           n,
    const int       jointNo, 
    const int       intervalNo, 
    double          splineParameters[5]) {

    double controlPoints[4];
    findIntervalControlPoints(x, n, jointNo, intervalNo, controlPoints);
    calculateIntervalSplineParams(controlPoints, duration, splineParameters);
    return;
}

/*  function to find the control points affecting the spline segment 
    in a particular interval for a particular joint.
    order is [C-3, C-2, C-1, C0]
*/
void TrajectoryOptimisationNLP::findIntervalControlPoints(
        const Number *  x,  
        Index           n,
        const int       jointNo, 
        const int       intervalNo, 
        double          controlPoints [4]) { //C0
        controlPoints[Cneg3] = x_VALUE(jointNo, intervalNo-3);
        controlPoints[Cneg2] = x_VALUE(jointNo, intervalNo-2);
        controlPoints[Cneg1] = x_VALUE(jointNo, intervalNo-1);
        controlPoints[C0] = x_VALUE(jointNo, intervalNo-0);
    }

/*  function to calculate the parameters of a spline segment
      (coefficients A-D and segment duration deltaT) from 
      control points C1-C4 and duration.
*/
void TrajectoryOptimisationNLP::calculateIntervalSplineParams(
    const double    controlPoints [4],
    const double    totalDuration,
    double          splineParameters[5]) {
    
    splineParameters[A] = ( controlPoints[C0] 
                            - 3 * controlPoints[Cneg1] 
                            + 3 * controlPoints[Cneg2] 
                            - controlPoints[Cneg3]
                          ) / 6;

    splineParameters[B] = ( 3 * controlPoints[Cneg1] 
                            - 6 * controlPoints[Cneg2] 
                            + 3 * controlPoints[Cneg3]
                          ) / 6;

    splineParameters[C] = ( 3 * controlPoints[Cneg1] 
                            - 3 * controlPoints[Cneg3]
                          ) / 6;

    splineParameters[D] = ( 1 * controlPoints[Cneg1]
                            + 4 * controlPoints[Cneg2]
                            + 1 * controlPoints[Cneg3]
                          ) / 6;

    splineParameters[deltaT] = totalDuration/m_numIntervals;
    }

double TrajectoryOptimisationNLP::calculatePos(const double splineParameters[5], const double tau) {
    return splineParameters[A]*std::pow(tau,3) + 
           splineParameters[B]*std::pow(tau,2) +
           splineParameters[C]*tau +
           splineParameters[D];
}

double TrajectoryOptimisationNLP::calculateVel(const double splineParameters[5], const double tau) {
    return (3*splineParameters[A]*std::pow(tau,2) + 
           2*splineParameters[B]*tau +
           splineParameters[C]) / splineParameters[deltaT];
}

double TrajectoryOptimisationNLP::calculateAcc(const double splineParameters[5], const double tau) {
    return (6*splineParameters[A]*tau + 
           2*splineParameters[B]) / (splineParameters[deltaT]*splineParameters[deltaT]);
}

    /*  method to convert end in tool space as 
        [posX, posY, posZ, rotX, rotY, rotZ] (applied in reverse order)
        aka 'pose' form
        into 
        [posX, posY, posZ, XnormX, XnormY, XnormZ, YnormX, YnormY, YnormZ] form
        aka 'posNorm' form.

        rotX, Y and Z are in degrees
    */
    std::vector<double> TrajectoryOptimisationNLP::pose2PosNorm(
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

    /*  method to convert end in tool space as 
        denavit-hartenberg form
        into 
        [posX, posY, posZ, XnormX, XnormY, XnormZ, YnormX, YnormY, YnormZ] form
        aka 'posNorm' form
    */
    std::vector<double> TrajectoryOptimisationNLP::DH2PosNorm(
        EricsTrajOpt::DenavitHartenberg pose
    ) {
        return rotTransMat2PosNorm(pose.getRotMat(), pose.getTransVec());
    }

    /*  method to convert pose in tool space as 
        rotation and translation matrices
        into 
        [posX, posY, posZ, XnormX, XnormY, XnormZ, YnormX, YnormY, YnormZ] form
        aka 'posNorm' form
    */
    std::vector<double> TrajectoryOptimisationNLP::rotTransMat2PosNorm(
        std::vector<std::vector<double>> rotMat,
        std::vector<double> transMat
    ) {
        std::vector<double> posNormVec = transMat;

        //multiplying the rotation vector by [1;0;0] and [0;1;0] is just selecting the 0th and 1st columns
        for(int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                posNormVec.push_back(rotMat[j][i]);
            }
        }
        return posNormVec;
    }

    std::vector<EricsTrajOpt::DenavitHartenberg> TrajectoryOptimisationNLP::getDHAtJointPose(
        std::vector<std::vector<double>> DHParamsAt0,
        double jointSpacePos[6]
    ) {
        //joint 3 needs to be flipped, for some reason
        jointSpacePos[2] = -jointSpacePos[2];

        std::vector<EricsTrajOpt::DenavitHartenberg> returnDH = {*new EricsTrajOpt::DenavitHartenberg(
            DHParamsAt0[0][0], 
            DHParamsAt0[0][1]+jointSpacePos[0], 
            DHParamsAt0[0][2], 
            DHParamsAt0[0][3])};
        for (int i = 1; i < 6; i++) {
            returnDH.push_back(returnDH[i-1] * (*new EricsTrajOpt::DenavitHartenberg(
            DHParamsAt0[i][0], 
            DHParamsAt0[i][1]+jointSpacePos[i], 
            DHParamsAt0[i][2], 
            DHParamsAt0[i][3])));
        }

        return returnDH;
    }

/*  calculate the values of the distance constraints.
    'distances' must be a pointer to the 1st index of the g array that holds a distance constraint.
*/
void TrajectoryOptimisationNLP::calcAllCollisions(
    const Number * x, 
    const Index n, 
    Number * distances
) {
    //some timing variables
    auto start = std::chrono::steady_clock::now();
    auto afterRotations = std::chrono::steady_clock::now();
    auto afterDistance = std::chrono::steady_clock::now();

    //variable to hold time and the size of the time step
    double time = 0;
    double dT = 1e-1;

    //variables to hold spline params for each joint for current interval
    double intervalSplineParams[6][5];

    //variables to hold arm pos in j-space and DH params for current time
    double pos[6];
    std::vector<EricsTrajOpt::DenavitHartenberg> DHParams;

    int indexOffset = 0;

    //variables to temporarily hold distance and index in distances
    double tempDistance;
    int tempIndex;

    //set all entries id distances to a high number
    for (int index = 0; index < m_numIntervals * (5+4+3+2+(6*(int)m_gripperObjects.size())); index++) {
        distances[index] = 2e19;
    }

    //iterate through all intervals
    for (int interval = 0; interval < m_numIntervals; interval++) {
        //get spline coefficients for this interval
        for (int joint = 0; joint < 6; joint++) {
            findIntervalCoefficients(x, n, joint, interval, intervalSplineParams[joint]);
        }
        //move in time steps of dT, as min distance points cannot be analytically determined
        time = 0;
        std::cout << "interval: " << interval << "\n";
        while (time < intervalSplineParams[0][5]) {
            std::cout << "time: " << time << "\n";
            start = std::chrono::steady_clock::now();
            //get arm joint positions at this time
            for (int joint = 0; joint < 6; joint++) {
                pos[joint] = calculatePos(intervalSplineParams[joint], time/intervalSplineParams[0][5]);
            }
            //rotate arm and gripper hulls to position
            DHParams = getDHAtJointPose(m_DHParamsAt0, pos);
            for (int link = 0; link < 6; link++) {
                m_armLinks[link].setHullDH(DHParams[link]);
            }
            for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                m_gripperObjects[gripperObj].setHullDH(DHParams[5]);
            }
            afterRotations = std::chrono::steady_clock::now();
            //calculate distance between static objects and arm links 2-6
            for (int link = 2; link < 7; link++) {
                tempDistance = calcStatic_2_linkCollision(link);
                tempIndex = indexOffset + 0 + link-2;
                //update distances if necessary
                if (tempDistance < distances[tempIndex]){
                    distances[tempIndex] = tempDistance;
                }
            }
            //calculate distance between arm link 1 and arm links 3-6
            for (int link = 3; link < 7; link++) {
                tempDistance = calcLink1_2_linkCollision(link);
                tempIndex = indexOffset + 5 + link-3;
                //update distances if necessary
                if (tempDistance < distances[tempIndex]){
                    distances[tempIndex] = tempDistance;
                }
            }
            //calculate distance between arm link 2 and arm links 4-6
            for (int link = 4; link < 7; link++) {
                tempDistance = calcLink2_2_linkCollision(link);
                tempIndex = indexOffset + 5+4 + link-4;
                //update distances if necessary
                if (tempDistance < distances[tempIndex]){
                    distances[tempIndex] = tempDistance;
                }
            }
            //calculate distance between arm link 3 and arm links 5-6
            for (int link = 5; link < 7; link++) {
                tempDistance = calcLink3_2_linkCollision(link);
                tempIndex = indexOffset + 5+4+3 + link-5;
                //update distances if necessary
                if (tempDistance < distances[tempIndex]){
                    distances[tempIndex] = tempDistance;
                }
            }
            //calculate distance between each gripper object and static objects plus arm links 1-5
            for (int gripperObj = 0; gripperObj < (int)m_gripperObjects.size(); gripperObj++) {
                //calculate distance between gripper object and static objects
                tempDistance = calcGripper_2_StaticCollision(gripperObj);
                tempIndex = indexOffset + 5+4+3+2 + (gripperObj*6);
                //update distances if necessary
                if (tempDistance < distances[tempIndex]){
                    distances[tempIndex] = tempDistance;
                }
                //calculate distance between gripper object and arm links 1-5
                for (int link = 1; link < 6; link++) {
                    tempDistance = calcGripper_2_linkCollision(link, gripperObj);
                    tempIndex = indexOffset + 5+4+3+2 + (gripperObj*6) + link;
                    //update distances if necessary
                    if (tempDistance < distances[tempIndex]){
                        distances[tempIndex] = tempDistance;
                    }
                }
            }
            afterDistance = std::chrono::steady_clock::now();
            auto rotateTime = afterRotations-start;
            auto distanceTime = afterDistance-afterRotations;
            std::cout << "rotation time: " << rotateTime.count() << "   distance time: " << distanceTime.count() << "\n";
            time+=dT;
        }
        indexOffset += (5+4+3+2+(6*(int)m_gripperObjects.size()));
    }
}

/*  calculate distance between specified arm link and static objects.
    arm links must first be rotated to position to be tested at.
    arm links are in the range 2-6
*/
double TrajectoryOptimisationNLP::calcStatic_2_linkCollision(
    int linkNo
) {
    double minDist = 2e19;
    for (int staticObj = 0; staticObj < (int)m_staticObstacles.size(); staticObj++) {
        obstDistSolver->OptimizeTNLP(staticObj_2_ArmLinks[linkNo-2][staticObj]);
        if (obstDistSolver->Statistics()->FinalObjective() < minDist) {
            minDist = obstDistSolver->Statistics()->FinalObjective();
        }
    }
    return minDist;
}

/*  calculate distance between specified arm link and arm link 1.
    arm links must first be rotated to position to be tested at.
    arm links are in the range 3-6
*/
double TrajectoryOptimisationNLP::calcLink1_2_linkCollision(
    int linkNo
) {
    obstDistSolver->OptimizeTNLP(link1_2_ArmLinks[linkNo-3]);
    return obstDistSolver->Statistics()->FinalObjective();
}

/*  calculate distance between specified arm link and arm link 2.
    arm links must first be rotated to position to be tested at.
    arm links are in the range 4-6
*/
double TrajectoryOptimisationNLP::calcLink2_2_linkCollision(
    int linkNo
) {
    obstDistSolver->OptimizeTNLP(link2_2_ArmLinks[linkNo-4]);
    return obstDistSolver->Statistics()->FinalObjective();
}

/*  calculate distance between specified arm link and arm link 2.
    arm links must first be rotated to position to be tested at.
    arm links are in the range 5-6
*/
double TrajectoryOptimisationNLP::calcLink3_2_linkCollision(
    int linkNo
) {
    obstDistSolver->OptimizeTNLP(link3_2_ArmLinks[linkNo-5]);
    return obstDistSolver->Statistics()->FinalObjective();
}

/*  calculate distance between specified arm link and specified gripper object.
    arm links must first be rotated to position to be tested at.
    arm links are in the range 1-5
    gripper objects are in the range 0-(size(gripperObjects)-1)
*/
double TrajectoryOptimisationNLP::calcGripper_2_linkCollision(
    int linkNo,
    int gripperObject
) {
    obstDistSolver->OptimizeTNLP(gripperObj_2_ArmLinks[linkNo-1][gripperObject]);
    return obstDistSolver->Statistics()->FinalObjective();
}

/*  calculate distance between specified arm link and static objects.
    arm links must first be rotated to position to be tested at.
    gripper objects are in the range 0-(size(gripperObjects)-1)
*/
double TrajectoryOptimisationNLP::calcGripper_2_StaticCollision(
    int gripperObject
) {
    double minDist = 2e19;
    for (int staticObj = 0; staticObj < (int)m_staticObstacles.size(); staticObj++) {
        obstDistSolver->OptimizeTNLP(gripperObj_2_staticObj[gripperObject][staticObj]);
        if (obstDistSolver->Statistics()->FinalObjective() < minDist) {
            minDist = obstDistSolver->Statistics()->FinalObjective();
        }
    }
    return minDist;
}



//function to graph the pos and vel of a solution, given the x array
void TrajectoryOptimisationNLP::graphSolution(
    const Number *  x,
    const Index     n
) {
    namespace plt = matplotlibcpp;
    double intervalDuration = duration/m_numIntervals;
    double splineParameters[5];
    std::vector<double> timeVec;
    std::vector<double> posVec[6];
    std::vector<double> velVec[6];
    std::string colors[6] = {"red", "orange", "yellow", "green", "blue", "purple"};
    for (Index interval = 0; interval < m_numIntervals; interval++) {
        for (double tau = 0; tau < 1; tau = tau + 0.001) {
            timeVec.push_back((interval+tau)*intervalDuration);
            for (Index joint = 0; joint < 6; joint ++) {
                findIntervalCoefficients(x, n, joint, interval, splineParameters);
                posVec[joint].push_back(calculatePos(splineParameters, tau));
                velVec[joint].push_back(calculateVel(splineParameters, tau));
            }
        }
    }
    plt::subplot(2,1,1);
    for (int joint = 0; joint < 6; joint ++) {
        plt::plot(timeVec, posVec[joint], colors[joint]);
    }
    plt::subplot(2,1,2);
    for (int joint = 0; joint < 6; joint ++) {
        plt::plot(timeVec, velVec[joint], colors[joint]);
    }
    plt::show();
    
}

