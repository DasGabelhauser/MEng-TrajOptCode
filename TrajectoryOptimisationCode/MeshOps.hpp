//Author: Eric Walker

#ifndef MESH_OPS
#define MESH_OPS

#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "DenavitHartenberg.hpp"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;


namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    /*  Class to store data about a triangle in 3d space */
    class Triangle {
    public:
        /* X Y and Z components of surface normal. */
        float Xn=0;
        float Yn=0;
        float Zn=0;
        /* X Y and Z components of vertices 1 2 and 3. */
        float X1=0;
        float Y1=0;
        float Z1=0;
        float X2=0;
        float Y2=0;
        float Z2=0;
        float X3=0;
        float Y3=0;
        float Z3=0;

        //Default constructor
        Triangle(){}

        /*  Construct a Triangle from 3*4 vector of floats.
            row is X, Y, Z, column is normal, vertex 1, vertex 2, vertex 3.
        */
        Triangle(std::vector<std::vector<float>> triangleData);

        /*  set Triangle data from 3*4 vector of floats.
            row is X, Y, Z, column is normal, vertex 1, vertex 2, vertex 3.
            similar to constructing a Triangle from the same, but with a preexisting Triangle.
        */
        void setTriangleData(std::vector<std::vector<float>> triangleData) {
            Xn = triangleData[0][0];
            Yn = triangleData[1][0];
            Zn = triangleData[2][0];
            X1 = triangleData[0][1];
            Y1 = triangleData[1][1];
            Z1 = triangleData[2][1];
            X2 = triangleData[0][2];
            Y2 = triangleData[1][2];
            Z2 = triangleData[2][2];
            X3 = triangleData[0][3];
            Y3 = triangleData[1][3];
            Z3 = triangleData[2][3];
        }

        std::vector<std::vector<float>> getTriangleData() {
            return *new std::vector<std::vector<float>> = {{Xn, X1, X2, X3},
                                                           {Yn, Y1, Y2, Y3},
                                                           {Zn, Z1, Z2, Z3}};
        }
    };

    /*  Class to represent a convex triangulated mesh.
        coplanar triangles are ignored. */
    class Hull {
    private:
        // a vector of Triangle objects representing the mesh after a transformation.
        std::vector<Triangle> m_TriangleList;
        // a vector of Triangle objects representing the mesh before a transformation.
        std::vector<Triangle> m_FixedTriangleList;

        //methods to add a traingle to the hull. used internally.
        bool addTriangle(char * triangleData);

        bool addTriangle(std::vector<std::vector<float>> triangleData);
    public:
        Hull(){}

        /*  Constructs a Hull object from a char array of binary STL triangle data 
            args:
            numTriangles - number of triangles in char array
            trainglesData - char array with all triangle data in binary STL format 
            (after the 80 char header and number of triangles entry). */
        Hull(int numTriangles, char * trianglesData);

        /*  Constructs a Hull object from two triangle list vectors
            args:
            triangleList - a list of triangles to initialise the hull's triangleList member with.
            fixedTriangleList - a list of triangles to initialise the hull's fixedTriangleList member with.
        */
        Hull(std::vector<Triangle> triangleList, std::vector<Triangle> fixedTriangleList) :
        m_TriangleList{triangleList}, m_FixedTriangleList{fixedTriangleList} {}

        /*  Constructs a Hull from a CGAL polyhedron object
            args:
            polyhedron - the CGAL polyhedron object to construct the Hull from
        */
        Hull(Polyhedron polyhedron);

        /*  Change the Hull's triangleList to the Hull's fixedTriangleList translated by translations  
            args:
            translations - a translation vector to translate fixedTriangleList by
        */
        void setHullTranslation(std::vector<double> translations);

        /*  Change the Hull's triangleList to the Hull's fixedTriangleList rotated by rotMat  
            args:
            rotMat - a rotation matrix to rotate fixedTriangleList by
        */
        void setHullRotation(std::vector<std::vector<double>> rotMat);

        /*  Change the Hull's triangleList to the Hull's fixedTriangleList transformed by DH  
            args:
            DH - a Denavit-Hartenberg matrix to rotate and translate fixedTriangleList by
        */
        void setHullDH(DenavitHartenberg DH);

        /*  Permanently transform the Hull to a new position and orientation, I.E., transform fixedTriangleList.
            rotations are performed in order Z, Y, X.
            Angles are measured in degrees.
            args:
            posX - the X component of the new position, relative to the previous origin.
            posY - the Y component of the new position, relative to the previous origin.
            posZ - the Z component of the new position, relative to the previous origin.
            rotX - the X component of the new orientation, relative to the previous orientation.
            rotY - the Y component of the new orientation, relative to the previous orientation.
            rotZ - the Z component of the new orientation, relative to the previous orientation.
        */
        void setHullPermPosRot(double posX, double posY, double posZ, double rotX, double rotY, double rotZ);

        Triangle getTriangle(int index) {
            return m_TriangleList[index];
        }

        Triangle getFixedTriangle(int index) {
            return m_FixedTriangleList[index];
        }

        std::vector<Triangle> getTriangleList() {
            return m_TriangleList;
        }

        int getNumTriangles() {
            return m_TriangleList.size();
        }

        /*  overloaded + operator. Concatenates the triangleLists and fixedTriangleLists of two Hulls
            into a single new hull.
        */
        Hull operator+ (Hull b) {
            // initialise two new triangleLists
            std::vector<Triangle> triangleList;
            std::vector<Triangle> fixedTriangleList;
            // add all of Hull A's triangles to the new triangleLists
            for (int i = 0; i < this->getNumTriangles(); i++) {
                triangleList.push_back(*new Triangle(this->getTriangle(i)));
                fixedTriangleList.push_back(*new Triangle(this->getFixedTriangle(i)));
            }
            // add all of Hull B's triangles to the new triangleLists
            for (int i = 0; i < b.getNumTriangles(); i++) {
                triangleList.push_back(*new Triangle(b.getTriangle(i)));
                fixedTriangleList.push_back(*new Triangle(b.getFixedTriangle(i)));
            }
            //create a new Hull adn return it
            return *new Hull(triangleList, fixedTriangleList);
        }

        /*  Exports the Hull to a binary STL file.
            args:
            filePathName - the absolute file path and file name of the STL file to write to/create.
        */
        void exportBinaryStl(std::string filePathName);
    };



    /*  Read in a binary STL and return a Hull object representing the file 
        filePathName: a string describing the absolute filepath and filename of the file to be read in.
    */
    Hull readBinaryStl(std::string filePathName);

    /*  Round a float value to the specified number of decimal places/
        value: the value to be rounded.
        decimalPlaces: the number of decimal places to round to.
    */
    float round_up(float value, int decimalPlaces);

    /*  Converts 4 chars from an array of chars into an int.
        args:
        charBuffer - the char array to read from
        position - the position in charBuffer to use as the first char.
        returns:
        the 4 chars, converted to an integer value
    */
    int readInt(char* charBuffer, int position);

    /*  Converts 4 chars from an array of chars into a float.
        args:
        charBuffer - the char array to read from
        position - the position in charBuffer to use as the first char.
        returns:
        the 4 chars, converted to a float value
    */
    float readFloat(char* charBuffer, int position);

    /*  Converts a rotation described as a rotation around Z, Y and X into a rotation matrix.
        Angles are measured in degrees.
        args:
        rotX - the X component of the rotation.
        rotY - the Y component of the rotation.
        rotZ - the Z component of the rotation.
        returns:
        a 2D vector containing the reulting rotation matrix.
    */
    std::vector<std::vector<double>> makeRotMat(double rotX, double rotY, double rotZ) ;

    
    
}

#endif //MESH_OPS