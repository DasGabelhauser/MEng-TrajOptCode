// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
//modified by Eric Walker

#ifndef TRAJECTORY_OPTIMISATION_NLP
#define TRAJECTORY_OPTIMISATION_NLP

#include "/usr/local/include/coin-or/IpTNLP.hpp"
#include <map>
#include <utility>
#include <vector>
#include "DenavitHartenberg.hpp"
#include "MeshOps.hpp"

#include "MeshDifferenceQP.hpp"
#include "/usr/local/include/coin-or/IpIpoptApplication.hpp"

using namespace Ipopt;

/* Construct an NLP to find the quickest route between two points, avoiding obstacles.
   */


class TrajectoryOptimisationNLP: public TNLP
{
public:
   

    /** Default constructor */
    TrajectoryOptimisationNLP(  int numIntervals, 
                                double * startPose, 
                                double * endPose, 
                                std::vector<EricsTrajOpt::Hull> staticObstacles,
                                std::vector<EricsTrajOpt::Hull> gripperObjects,
                                std::vector<EricsTrajOpt::Hull> armLinks,
                                std::vector<std::vector<double>> DHParamsAt0,
                                std::vector<std::vector<double>> &returnX,
                                double maxPos = 360, 
                                double maxVel = 180, 
                                double maxAcc = 40);

    /** Default destructor */
    virtual ~TrajectoryOptimisationNLP();

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the NLP */
    virtual bool get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
    );

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(
        Index   n,
        Number* x_l,
        Number* x_u,
        Index   m,
        Number* g_l,
        Number* g_u
    );

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(
        Index   n,
        bool    init_x,
        Number* x,
        bool    init_z,
        Number* z_L,
        Number* z_U,
        Index   m,
        bool    init_lambda,
        Number* lambda
    );

    /** Method to return the objective value */
    virtual bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
    );

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
    );

    /** Method to return the constraint residuals */
    virtual bool eval_g(
       Index         n,
       const Number* x,
       bool          new_x,
       Index         m,
       Number*       g
    );

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Index         nele_jac,
        Index*        iRow,
        Index*        jCol,
        Number*       values
    );

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    virtual bool eval_h(
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
        Number*       values
    );

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
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
    );
    //@}

    // function to calculate the pos value of a spline segment given by
    // cubic coefficients [A-D and segment duration deltaT] at time tau
    static double calculatePos(const double splineParameters[5], const double tau);

    // function to calculate the vel value of a spline segment given by
    // cubic coefficients [A-D and segment duration deltaT] at time tau
    static double calculateVel(const double splineParameters[5], const double tau);

    // function to calculate the acc value of a spline segment given by
    // cubic coefficients [A-D and segment duration deltaT] at time tau
    static double calculateAcc(const double splineParameters[5], const double tau);

private:

    //extra information about the problem we are solving
    //number of spline intervals
    const int m_numIntervals;
    //number of spline control points (will be m_numIntervals + 3)
    const int m_numControlPoints;
    //the desired start pose of the arm in joint space
    const double m_startPose[6];
    //the desired end pose of the arm in joint space //TODO make tool-space
    const double m_endPose[6];

    //vectors of hulls representing different collections of collidable objects
    std::vector<EricsTrajOpt::Hull> m_staticObstacles;
    std::vector<EricsTrajOpt::Hull> m_gripperObjects;
    std::vector<EricsTrajOpt::Hull> m_armLinks;

    //a vector containing the DH parameters from each link to the last with all joints at 0
    std::vector<std::vector<double>> m_DHParamsAt0;

    //a vector to return the design variables in
    std::vector<std::vector<double>> m_returnX;

    //limits for pos vel and acc. minimum limits are assumed to be opposite.
    const double m_maxPos;
    const double m_maxVel;
    const double m_maxAcc;
    //map to map a pair of indices representing an entry in the hessian matrix to
    //the relevant index of the triplet-formatted hessian arrays
    // key pair is <j(col), i(row)>
    std::map<std::pair<int, int>, int> hessianIndexMap;

    //vectors/arrays of meshDifferenceQPs to calculate distances between various obstacles
    std::vector<Ipopt::SmartPtr<Ipopt::TNLP>> staticObj_2_ArmLinks [5];
    Ipopt::SmartPtr<Ipopt::TNLP> link1_2_ArmLinks [4];
    Ipopt::SmartPtr<Ipopt::TNLP> link2_2_ArmLinks [3];
    Ipopt::SmartPtr<Ipopt::TNLP> link3_2_ArmLinks [2];
    std::vector<Ipopt::SmartPtr<Ipopt::TNLP>> gripperObj_2_ArmLinks [5];
    std::vector<std::vector<Ipopt::SmartPtr<Ipopt::TNLP>>> gripperObj_2_staticObj;

    //solver to solve above MeshDifferenceQPs
    SmartPtr<IpoptApplication> obstDistSolver;



    // function to find the minimum and maximum pos of a spline within an interval,
    // given spline coefficients (A-D, deltaT) 
    void findMaxMinPos (
        const double    splineParameters[5],
        double *        maxPos, 
        double *        minPos);
    
    // function to find the minimum and maximum vel of a spline within an interval,
    // given spline coefficients (A-D, deltaT) 
    void findMaxMinVel (
        const double    splineParameters[5],
        double *        maxVel, 
        double *        minVel);

    // function to find the minimum and maximum acc of a spline within an interval,
    // given spline coefficients (A-D, deltaT) 
    void findMaxMinAcc (
        const double    splineParameters[5],
        double *        maxAcc, 
        double *        minAcc);

    /*  function to differentiate the minimum and maximum pos of a spline within an interval,
        given spline control points C-3-C0
        order is:
        [dC-3, dC-2, dC-1, dC0] */
    void differentiateMaxMinPos (
        const double    controlPoints [4],
        double          dMaxPos[4], 
        double          dMinPos[4]);
    
    /*  function to differentiate the minimum and maximum vel of a spline within an interval,
        given spline control points C-3-C0 and segment duration (deltaT)
        order is:
        [dC-3, dC-2, dC-1, dC0, dT] */
    void differentiateMaxMinVel (
        const double    controlPoints [4],
        const double    totalDuration,
        double          dMaxVel[5], 
        double          dMinVel[5]);

    /*  function to differentiate the minimum and maximum acc of a spline within an interval,
        given spline control points C-3-C0 and segment duration (deltaT)
        order is:
        [dC-3, dC-2, dC-1, dC0, dT] */
    void differentiateMaxMinAcc (
        const double    controlPoints [4],
        const double    totalDuration, 
        double          dMaxAcc[5], 
        double          dMinAcc[5]);


    /*function to find the hessian of MaxMinPos for a particular joint and interval
      only the lower left of the return matrices are filled
    order is:
    [dC-3dC-3, dC-2dC-3, dC-1dC-3, dC0dC-3]
    [dC-3dC-2, dC-2dC-2, dC-1dC-2, dC0dC-2]
    [dC-3dC-1, dC-2dC-1, dC-1dC-1, dC0dC-1]
    [dC-3dC0 , dC-2dC0 , dC-1dC0 , dC0dC0 ] */

    void hessianMaxMinPos (
        const double    controlPoints [4],
        double          dMaxPos[4][4], 
        double          dMinPos[4][4]);
    
    /*function to find the hessian of MaxMinVel for a particular joint and interval
      only the lower left of the return matrices are filled
    order is:
    [dC-3dC-3, dC-2dC-3, dC-1dC-3, dC0dC-3, dCTdC-3]
    [dC-3dC-2, dC-2dC-2, dC-1dC-2, dC0dC-2, dCTdC-2]
    [dC-3dC-1, dC-2dC-1, dC-1dC-1, dC0dC-1, dCTdC-1]
    [dC-3dC0 , dC-2dC0 , dC-1dC0 , dC0dC0 , dCTdC0 ]
    [dC-3dT  , dC-2dT  , dC-1dT  , dC0dT  , dCTdT  ] */
    void hessianMaxMinVel (
        const double    controlPoints [4],
        const double    totalDuration, 
        double          dMaxVel[5][5], 
        double          dMinVel[5][5]);

    /*function to find the hessian of MaxMinAcc for a particular joint and interval
      only the lower left of the return matrices are filled
    order is:
    [dC-3dC-3, dC-2dC-3, dC-1dC-3, dC0dC-3, dCTdC-3]
    [dC-3dC-2, dC-2dC-2, dC-1dC-2, dC0dC-2, dCTdC-2]
    [dC-3dC-1, dC-2dC-1, dC-1dC-1, dC0dC-1, dCTdC-1]
    [dC-3dC0 , dC-2dC0 , dC-1dC0 , dC0dC0 , dCTdC0 ]
    [dC-3dT  , dC-2dT  , dC-1dT  , dC0dT  , dCTdT  ] */
    void hessianMaxMinAcc (
        const double    controlPoints [4],
        const double    totalDuration, 
        double          dMaxAcc[5][5], 
        double          dMinAcc[5][5]);

    void differentiateEndPos(
        const Number *  x, 
        Index           n,
        double *        endPosDerivatives);
    
    // function to find the coefficients of the spline segment 
    // in a particular interval for a particular joint.
    void findIntervalCoefficients(
        const Number *  x,  
        Index           n,
        const int       jointNo, 
        const int       intervalNo, 
        double          splineParameters[5]);

    /*  function to find the control points affecting the spline segment 
        in a particular interval for a particular joint.
        order is [C0, C-1, C-2, C-3]
    */
    void findIntervalControlPoints(
        const Number *  x,  
        Index           n,
        const int       jointNo, 
        const int       intervalNo, 
        double          controlPoints [4]);
    
    //function to calculate the parameters of a spline segment
    // (coefficients A-D and segment duration deltaT) from 
    // control points C1-C4.
    void calculateIntervalSplineParams(
        const double    controlPoints [4],
        const double    totalDuration,
        double          splineParameters[5]);


    //function to graph the pos and vel of a solution, given the x array
    void graphSolution(
        const Number *  x,
        const Index     n
    );

    /*  method to convert pose in tool space as 
        [posX, posY, posZ, rotX, rotY, rotZ] (applied in reverse order)
        aka 'pose' form
        into 
        [posX, posY, posZ, XnormX, XnormY, XnormZ, YnormX, YnormY, YnormZ] form
        aka 'posNorm' form
    */
    std::vector<double> pose2PosNorm(
        double posX, 
        double posY, 
        double posZ, 
        double rotX, 
        double rotY, 
        double rotZ
    );

    /*  method to convert pose in tool space as 
        denavit-hartenberg form
        into 
        [posX, posY, posZ, XnormX, XnormY, XnormZ, YnormX, YnormY, YnormZ] form
        aka 'posNorm' form
    */
    std::vector<double> DH2PosNorm(
        EricsTrajOpt::DenavitHartenberg pose
    );

    /*  method to convert pose in tool space as 
        rotation and translation matrices
        into 
        [posX, posY, posZ, XnormX, XnormY, XnormZ, YnormX, YnormY, YnormZ] form
        aka 'posNorm' form
    */
    std::vector<double> rotTransMat2PosNorm(
        std::vector<std::vector<double>> rotMat,
        std::vector<double> transMat
    );

    std::vector<EricsTrajOpt::DenavitHartenberg> getDHAtJointPose(
        std::vector<std::vector<double>> DHParamsAt0,
        double jointSpacePos[6]
    );

    void calcAllCollisions(
        const Number * x, 
        const Index n, 
        Number * distances
    );

    /*  calculate distance between specified arm link and static objects.
        arm links must first be rotated to position to be tested at
        arm links are in the range 2-6
    */
    double calcStatic_2_linkCollision(
        int linkNo
    );

    /*  calculate distance between specified arm link and arm link 1.
        arm links must first be rotated to position to be tested at.
        arm links are in the range 3-6
    */
    double calcLink1_2_linkCollision(
        int linkNo
    );

    /*  calculate distance between specified arm link and arm link 2.
        arm links must first be rotated to position to be tested at.
        arm links are in the range 4-6
    */
    double calcLink2_2_linkCollision(
        int linkNo
    );

    /*  calculate distance between specified arm link and arm link 2.
        arm links must first be rotated to position to be tested at.
        arm links are in the range 5-6
    */
    double calcLink3_2_linkCollision(
        int linkNo
    );

    /*  calculate distance between specified arm link and specified gripper object.
        arm links must first be rotated to position to be tested at.
        arm links are in the range 1-5
        gripper objects are in the range 0-(size(gripperObjects)-1)
    */
    double calcGripper_2_linkCollision(
        int linkNo,
        int gripperObject
    );

    /*  calculate distance between specified arm link and static objects.
        arm links must first be rotated to position to be tested at.
        gripper objects are in the range 0-(size(gripperObjects)-1)
    */
    double calcGripper_2_StaticCollision(
        int gripperObject
    );

    
    

    /**@name Methods to block default compiler methods.
     *
     *  The compiler automatically generates the following three methods.
     *  Since the default compiler implementation is generally not what
     *  you want (for all but the most simple classes), we usually
     *  put the declarations of these methods in the private section
     *  and never implement them. This prevents the compiler from
     *  implementing an incorrect "default" behavior without us
     *  knowing. (See Scott Meyers book, "Effective C++")
     */
    //@{
    TrajectoryOptimisationNLP(
        const TrajectoryOptimisationNLP&
    );

    TrajectoryOptimisationNLP& operator=(
        const TrajectoryOptimisationNLP&
    );
    //@}
};

#endif //TRAJECTORY_OPTIMISATION_NLP
