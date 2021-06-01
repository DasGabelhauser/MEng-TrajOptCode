// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
//modified by Eric Walker

#include "MeshDifferenceQP.hpp"

#include <cassert>
#include <iostream>
#include <cmath>

using namespace Ipopt;

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

//#define DEBUG_PRINT_ALL_FUNCTIONS

// constructor
Mesh_Difference_QP::Mesh_Difference_QP(EricsTrajOpt::Hull* hull1, EricsTrajOpt::Hull* hull2) {
    this->hull1 = hull1;
    this->hull2 = hull2;
}

// destructor
Mesh_Difference_QP::~Mesh_Difference_QP()
{ }

// [TNLP_get_nlp_info]
// returns the size of the problem
bool Mesh_Difference_QP::get_nlp_info(
    Index&          n,
    Index&          m,
    Index&          nnz_jac_g,
    Index&          nnz_h_lag,
    IndexStyleEnum& index_style
)
{
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP get_nlp_info start\n";
    #endif
    // The problem has 
    n = 6;
 
    // num of constraints equal to num of triangles in each hull
    m = hull1->getNumTriangles() + hull2->getNumTriangles();
 
    // nnz in jacobian is 3 * the number of triangles in each hull.
    // some might be 0 at some rotations, but this simplifies things
    nnz_jac_g = (hull1->getNumTriangles() + hull2->getNumTriangles()) * 3;
 
    // the Hessian is also dense and has 36 total nonzeros, but we
    // only need the lower left corner (since it is symmetric)
    nnz_h_lag = 21;
 
    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP get_nlp_info end\n";
    #endif
    return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool Mesh_Difference_QP::get_bounds_info(
    Index   n,
    Number* x_l,
    Number* x_u,
    Index   m,
    Number* g_l,
    Number* g_u
) {
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP get_bounds_info start\n";
    #endif
    for (int i = 0; i < n; i++) {
        //no constraints on variables
        x_l[i] = -2e19;
        x_u[i] = 2e19;
    }

    for (int i = 0; i < hull1->getNumTriangles(); i++) {
        // all constraints have a lower bound of -inf
        g_l[i] = -2e19;
        // take elements of g_u from values in meshDifferenceConstraints
        g_u[i] =    hull1->getTriangle(i).Xn * hull1->getTriangle(i).X1 +
                    hull1->getTriangle(i).Yn * hull1->getTriangle(i).Y1 +
                    hull1->getTriangle(i).Zn * hull1->getTriangle(i).Z1;
    }
    for (int i = 0; i < hull2->getNumTriangles(); i++) {
        // all constraints have a lower bound of -inf
        g_l[i + hull1->getNumTriangles()] = -2e19;
        // take elements of g_u from values in meshDifferenceConstraints
        g_u[i + hull1->getNumTriangles()] =    
                    hull2->getTriangle(i).Xn * hull2->getTriangle(i).X1 +
                    hull2->getTriangle(i).Yn * hull2->getTriangle(i).Y1 +
                    hull2->getTriangle(i).Zn * hull2->getTriangle(i).Z1;
    }
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP get_bounds_info end\n";
    #endif
    return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool Mesh_Difference_QP::get_starting_point(
    Index   n,
    bool    init_x,
    Number* x,
    bool    init_z,
    Number* z_L,
    Number* z_U,
    Index   m,
    bool    init_lambda,
    Number* lambda
) 
{
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP get_starting_point start\n";
    #endif
    // Here, we assume we only have starting values for x
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
    // initialize x to the first point of the first triangle of each hull 
    x[0] = hull1->getTriangle(0).X1;
    x[1] = hull1->getTriangle(0).Y1;
    x[2] = hull1->getTriangle(0).Z1;
    x[3] = hull2->getTriangle(0).X1;
    x[4] = hull2->getTriangle(0).Y1;
    x[5] = hull2->getTriangle(0).Z1;
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP get_starting_point end\n";
    #endif
    return true;

}
// [TNLP_get_starting_point]

// [TNLP_eval_f]
// returns the value of the objective function
bool Mesh_Difference_QP::eval_f(
    Index         n,
    const Number* x,
    bool          new_x,
    Number&       obj_value
)
{
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_f start\n";
    #endif
    obj_value = std::sqrt(  std::pow(x[3] - x[0], 2) + 
                            std::pow(x[4] - x[1], 2) + 
                            std::pow(x[5] - x[2], 2));
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_f end\n";
    #endif
    return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool Mesh_Difference_QP::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_grad_f start\n";
    #endif
    assert(n == 6);

    double denominator = std::sqrt( std::pow(x[3] - x[0], 2) + 
                                    std::pow(x[4] - x[1], 2) + 
                                    std::pow(x[5] - x[2], 2));
    
    //prevent divide by zero error
    if (denominator == 0) {
        grad_f[0] = 0;
        grad_f[1] = 0;
        grad_f[2] = 0;
    //if dividing by zero will not happen, calculate values
    } else {
        grad_f[0] = (x[0]-x[3])/std::sqrt(  std::pow(x[3] - x[0], 2) + 
                                            std::pow(x[4] - x[1], 2) + 
                                            std::pow(x[5] - x[2], 2));
        grad_f[1] = (x[1]-x[4])/std::sqrt(  std::pow(x[3] - x[0], 2) + 
                                            std::pow(x[4] - x[1], 2) + 
                                            std::pow(x[5] - x[2], 2));
        grad_f[2] = (x[2]-x[5])/std::sqrt(  std::pow(x[3] - x[0], 2) + 
                                            std::pow(x[4] - x[1], 2) + 
                                            std::pow(x[5] - x[2], 2));
    }
    grad_f[3] = -grad_f[0];
    grad_f[4] = -grad_f[1];
    grad_f[5] = -grad_f[2];
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_grad_f end\n";
    #endif
    return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool Mesh_Difference_QP::eval_g(
    Index         n,
    const Number* x,
    bool          new_x,
    Index         m,
    Number*       g
)
{
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_g start\n";
    #endif
    /*  calculate elements of g using values stored in meshDifferenceConstraints
        and current values of x. */ 
    for (int i = 0; i < hull1->getNumTriangles(); i++) {
        g[i] =  hull1->getTriangle(i).Xn * x[0] +
                hull1->getTriangle(i).Yn * x[1] +
                hull1->getTriangle(i).Zn * x[2];
    }
    int offset = hull1->getNumTriangles();
    for (int i = 0; i < hull2->getNumTriangles(); i++) {
        g[i+offset] =   hull2->getTriangle(i).Xn * x[3] +
                        hull2->getTriangle(i).Yn * x[4] +
                        hull2->getTriangle(i).Zn * x[5];
    }
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_g end\n";
    #endif
    return true;
}
// [TNLP_eval_g]

// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool Mesh_Difference_QP::eval_jac_g(
    Index         n,
    const Number* x,
    bool          new_x,
    Index         m,
    Index         nele_jac,
    Index*        iRow,
    Index*        jCol,
    Number*       values
) {
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_jac_g start\n";
    #endif
    int xOffset = 3;
    int gOffset = hull1->getNumTriangles();
    int jacOffset = hull1->getNumTriangles() * 3;

    if( values == NULL )
    {
        // copy the structure of the Jacobian
        for (int i = 0; i < hull1->getNumTriangles(); i++) {
            for (int j = 0; j < 3; j++) {
                iRow[i*3 + j] = i;
                jCol[i*3 + j] = j;
            }
        }
        for (int i = 0; i < hull2->getNumTriangles(); i++) {
            for (int j = 0; j < 3; j++) {
                iRow[i*3 + j + jacOffset] = i + gOffset;
                jCol[i*3 + j + jacOffset] = j + xOffset;
            }
        }

    }
    else
    {
        //std::cout << "evalling jac g\n";
        // copy the values of the Jacobian of the constraints
        for (int i = 0; i < hull1->getNumTriangles(); i++) {
            values[i*3 + 0] = hull1->getTriangle(i).Xn;
            values[i*3 + 1] = hull1->getTriangle(i).Yn;
            values[i*3 + 2] = hull1->getTriangle(i).Zn; 
        }
        for (int i = 0; i < hull2->getNumTriangles(); i++) {
            values[i*3 + 0 + jacOffset] = hull2->getTriangle(i).Xn;
            values[i*3 + 1 + jacOffset] = hull2->getTriangle(i).Yn;
            values[i*3 + 2 + jacOffset] = hull2->getTriangle(i).Zn; 
        }
    }
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_jac_g end\n";
    #endif
    return true;
}
// [TNLP_eval_jac_g]

// [TNLP_eval_h]
//return the structure or values of the Hessian
bool Mesh_Difference_QP::eval_h(
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
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_h start\n";
    #endif
    assert(n == 6);

    if( values == NULL )
    {
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.

        // the hessian for this problem is actually dense
        Index idx = 0;
        for( Index row = 0; row < n; row++ ) {
            for( Index col = 0; col <= row; col++ ) {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }
        assert(idx == nele_hess);
    }
    else {
        // return the values. This is a symmetric matrix, fill the lower left
        // triangle only
        double denominator = std::sqrt(std::pow(std::pow(x[3] - x[0], 2) + 
                                                std::pow(x[4] - x[1], 2) + 
                                                std::pow(x[5] - x[2], 2),
                                                3));
        //prevent divide by zero error
        if (denominator == 0) {
            for (int i = 0; i < 21; i++) {
                values[i] = 0;
            }
        //if divde by zero will not happen, calculate values
        } else {
            // fill the objective portion
            values[0] = obj_factor * 
                        (x[1]*x[1] + x[2]*x[2] + x[4]*x[4] + x[5]*x[5] - 2*x[1]*x[4] - 2*x[2]*x[5]) /
                        denominator; // 0,0

            values[1] = obj_factor * ((x[0]-x[3])*(x[4]-x[1])) / denominator;   // 1,0
            values[2] = obj_factor * 
                        (x[0]*x[0] + x[2]*x[2] + x[3]*x[3] + x[5]*x[5] - 2*x[0]*x[3] - 2*x[2]*x[5]) /
                        denominator; // 1,1

            values[3] = obj_factor * ((x[0]-x[3])*(x[5]-x[2])) / denominator;   // 2,0
            values[4] = obj_factor * ((x[1]-x[4])*(x[5]-x[2])) / denominator;   // 2,1
            values[5] = obj_factor * 
                        (x[0]*x[0] + x[1]*x[1] + x[3]*x[3] + x[4]*x[4] - 2*x[0]*x[3] - 2*x[1]*x[4]) /
                        denominator; // 2,2

            values[6] = -values[0]; // 3,0 (-0,0)
            values[7] = -values[1]; // 3,1 (-1,0, or -0,1)
            values[8] = -values[3]; // 3,2 (-2,0, or -0,2)
            values[9] = values[0];  // 3,3 (same as 0,0)

            values[10] = -values[1]; // 4,0 (-1,0)
            values[11] = -values[2]; // 4,1 (-1,1)
            values[12] = -values[4]; // 4,2 (-2,1, or -1,2)
            values[13] = -values[7]; // 4,3 (-3,1, or -1,3)
            values[14] = values[2];  // 4,4 (same as 1,1)

            values[15] = -values[3]; // 5,0 (-2,0)
            values[16] = -values[4]; // 5,1 (-2,1)
            values[17] = -values[5]; // 5,2 (-2,2)
            values[18] = -values[8]; // 5,3 (-3,2, or -2,3)
            values[19] = -values[12]; // 5,4 (-4,2, or -2,4)
            values[20] = values[5];  // 5,5 (same as 2,2)
        }
    }
    #ifdef DEBUG_PRINT_ALL_FUNCTIONS
    std::cout << "MeshDifferenceQP eval_h end\n";
    #endif
    return true;
}
// [TNLP_eval_h]

// [TNLP_finalize_solution]
void Mesh_Difference_QP::finalize_solution(
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
    //this function doesnt need to do anything for actual operation, 
    //but these lines are useful to have around for debug purposes.

    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.
 
    // For this example, we write the solution to the console
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
 
    /*std::cout << std::endl << "Final value of the constraints:" << std::endl;
    for( Index i = 0; i < m; i++ )
    {
       std::cout << "g(" << i << ") = " << g[i] << std::endl;
    }*/
} 
// [TNLP_finalize_solution]
