//Author: Eric Walker
#include "MeshOps.hpp"
using namespace EricsTrajOpt;
namespace EricsTrajOpt {
    /*  Construct a Triangle from 3*4 vector of floats.
        row is X, Y, Z, column is normal, vertex 1, vertex 2, vertex 3.
    */
    Triangle::Triangle(std::vector<std::vector<float>> triangleData) :
        //copy all data to correct member.
        Xn{triangleData[0][0]},
        Yn{triangleData[1][0]},
        Zn{triangleData[2][0]},
        X1{triangleData[0][1]},
        Y1{triangleData[1][1]},
        Z1{triangleData[2][1]},
        X2{triangleData[0][2]},
        Y2{triangleData[1][2]},
        Z2{triangleData[2][2]},
        X3{triangleData[0][3]},
        Y3{triangleData[1][3]},
        Z3{triangleData[2][3]}
    {}

    /*  Constructs a Hull object from a char array of binary STL triangle data 
        args:
        numTriangles - number of triangles in char array
        trainglesData - char array with all triangle data in binary STL format 
        (after the 80 char header and number of triangles entry). */
    Hull::Hull(int numTriangles, char * trianglesData) {
        // Iterate through all triangles read in from file
        for (int i = 0; i < numTriangles; i++) {
            // Add to triangleList. Each triangle has 50 bytes of data
            addTriangle(trianglesData + (50*i));
        }
    }


    /*  Constructs a Hull from a CGAL Polyhedron_3 object
        args:
        polyhedron - the CGAL polyhedron object to construct the Hull from
    */
    Hull::Hull(Polyhedron polyhedron) {
        //ensure polyhedron is triangular
        CGAL::Polygon_mesh_processing::triangulate_faces(polyhedron);
        std::vector<std::vector<float>> triangleData;
        //std::vector<std::vector<double>> vertices;
        double normalisingFactor;
        Polyhedron::Halfedge_around_facet_circulator facCirc;

        //iterate through each face
        for (Polyhedron::Facet_iterator facet = polyhedron.facets_begin(); facet != polyhedron.facets_end(); ++facet) {
            //iterate through each vertex of the face
            triangleData = *new std::vector<std::vector<float>>({   {0,0,0,0},
                                                                    {0,0,0,0},
                                                                    {0,0,0,0}   });
            //iterate and set facet coordinates
            facCirc = facet->facet_begin();
            do {
                triangleData[0][(int)std::distance(facet->facet_begin(), facCirc) + 1] = (float)CGAL::to_double(facCirc->vertex()->point().x());
                triangleData[1][(int)std::distance(facet->facet_begin(), facCirc) + 1] = (float)CGAL::to_double(facCirc->vertex()->point().y());
                triangleData[2][(int)std::distance(facet->facet_begin(), facCirc) + 1] = (float)CGAL::to_double(facCirc->vertex()->point().z());
            } while ( ++facCirc != facet->facet_begin());
            //calculate facet normal
            triangleData[0][0] =    ((triangleData[1][2]-triangleData[1][1])*(triangleData[2][3]-triangleData[2][1]))-
                                    ((triangleData[1][3]-triangleData[1][1])*(triangleData[2][2]-triangleData[2][1]));
            triangleData[1][0] =    ((triangleData[2][2]-triangleData[2][1])*(triangleData[0][3]-triangleData[0][1]))-
                                    ((triangleData[2][3]-triangleData[2][1])*(triangleData[0][2]-triangleData[0][1]));
            triangleData[2][0] =    ((triangleData[0][2]-triangleData[0][1])*(triangleData[1][3]-triangleData[1][1]))-
                                    ((triangleData[0][3]-triangleData[0][1])*(triangleData[1][2]-triangleData[1][1]));

            //normalise facet normal
            normalisingFactor = std::sqrt((triangleData[0][0]*triangleData[0][0])+
                                          (triangleData[1][0]*triangleData[1][0])+
                                          (triangleData[2][0]*triangleData[2][0]));
            triangleData[0][0] = triangleData[0][0]/normalisingFactor;
            triangleData[1][0] = triangleData[1][0]/normalisingFactor;
            triangleData[2][0] = triangleData[2][0]/normalisingFactor;

            //add triangle to hull
            addTriangle(triangleData);
        }
    }

    /*  Add a triangle to the Hull using a 50 char array of triangle data. */
    bool Hull::addTriangle(char * triangleData){
        //convert char triangleData into a 2D vector of float triangleData
        std::vector<std::vector<float>> floatTriangleData =  {
            {   readFloat(triangleData, 0 ), readFloat(triangleData, 12), readFloat(triangleData, 24), readFloat(triangleData, 36) },
            {   readFloat(triangleData, 4 ), readFloat(triangleData, 16), readFloat(triangleData, 28), readFloat(triangleData, 40) },
            {   readFloat(triangleData, 8 ), readFloat(triangleData, 20), readFloat(triangleData, 32), readFloat(triangleData, 44) }
        };
        // use float triangleData to add new triangle
        return addTriangle(floatTriangleData);
    }
    
    /*  Add a triangle to the Hull using a 3*4 float vector of triangle data. */
    bool Hull::addTriangle(std::vector<std::vector<float>> triangleData) {
        // check if triangle is coplanar with other triangle
        for (std::size_t i = 0; i < m_TriangleList.size(); i++) {
            if (m_FixedTriangleList[i].Xn == triangleData[0][0] and
                m_FixedTriangleList[i].Yn == triangleData[1][0] and
                m_FixedTriangleList[i].Zn == triangleData[2][0]) {
                return false;
            }
        }
        //if loop is escaped without returning, no triangle sharing normal is in hull already, add triangle
        m_FixedTriangleList.push_back(* new Triangle(triangleData));
        m_TriangleList.push_back(* new Triangle(triangleData));
        return true;
    }

    /*  Change the Hull's triangleList to the Hull's fixedTriangleList translated by translations  
        args:
        translations - a translation vector to translate fixedTriangleList by
    */
    void Hull::setHullTranslation(std::vector<double> translations) {
        std::vector<std::vector<float>> tempTriangleData;
        //translate each triangle
        for (std::size_t triangle = 0; triangle < m_TriangleList.size(); triangle++) {
            tempTriangleData = m_FixedTriangleList[triangle].getTriangleData();
            //translate only the vertex coords, not the normal vector
            for (int vertex = 0; vertex < 3; vertex++) {
                for (int coord = 0; coord < 3; coord++) {
                    tempTriangleData[coord][vertex+1] += translations[coord];
                }
            }
            m_TriangleList[triangle].setTriangleData(tempTriangleData);
        }
    }

    /*  Change the Hull's triangleList to the Hull's fixedTriangleList rotated by rotMat  
        args:
        rotMat - a rotation matrix to rotate fixedTriangleList by
    */
    void Hull::setHullRotation(std::vector<std::vector<double>> rotMat) {
        std::vector<std::vector<float>> oldTriangleData;
        std::vector<std::vector<float>> newTriangleData;
        //rotate each triangle
        for (std::size_t triangle = 0; triangle < m_TriangleList.size(); triangle++) {
            oldTriangleData = m_FixedTriangleList[triangle].getTriangleData();
            newTriangleData = { {0,0,0,0},
                                {0,0,0,0},
                                {0,0,0,0} };
            //rotate all vertices and normal vector
            for (int vertex = 0; vertex < 4; vertex++) { // col
                for (int coord = 0; coord < 3; coord++) { // row
                    for (int i = 0; i < 3; i++) {
                        newTriangleData[vertex][coord] += rotMat[coord][i]*oldTriangleData[i][vertex];
                    }
                }
            }
            m_TriangleList[triangle].setTriangleData(newTriangleData);
        }
    }

    /*  Change the Hull's triangleList to the Hull's fixedTriangleList transformed by DH  
        args:
        DH - a Denavit-Hartenberg matrix to rotate and translate fixedTriangleList by
    */
    void Hull::setHullDH(DenavitHartenberg DH) {
        std::vector<double> translations = DH.getTransVec();
        std::vector<std::vector<double>> rotMat = DH.getRotMat();
        std::vector<std::vector<float>> oldTriangleData;
        std::vector<std::vector<float>> newTriangleData;
        //transform each triangle
        for (std::size_t triangle = 0; triangle < m_TriangleList.size(); triangle++) {
            oldTriangleData = m_FixedTriangleList[triangle].getTriangleData();
            newTriangleData = { {0,0,0,0},
                                {0,0,0,0},
                                {0,0,0,0} };
            //rotate each vertex and normal vector
            for (int vertex = 0; vertex < 4; vertex++) { // col
                for (int coord = 0; coord < 3; coord++) { // row
                    for (int i = 0; i < 3; i++) {
                        newTriangleData[coord][vertex] += rotMat[coord][i]*oldTriangleData[i][vertex];
                    }
                }
            }
            //translate only the vertex coords, not the normal vector
            for (int vertex = 1; vertex < 4; vertex++) { // col
                for (int coord = 0; coord < 3; coord++) { // row
                    newTriangleData[coord][vertex] += translations[coord];
                }
            }
            m_TriangleList[triangle].setTriangleData(newTriangleData);
        }
    }

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
    void Hull::setHullPermPosRot(double posX, double posY, double posZ, double rotX, double rotY, double rotZ) {
        std::vector<double> translations = {posX, posY, posZ};
        std::vector<std::vector<double>> rotMat = makeRotMat(rotX, rotY, rotZ);
        std::vector<std::vector<float>> oldTriangleData;
        std::vector<std::vector<float>> newTriangleData;
        //transform each triangle
        for (std::size_t triangle = 0; triangle < m_FixedTriangleList.size(); triangle++) {
            oldTriangleData = m_FixedTriangleList[triangle].getTriangleData();
            newTriangleData = { {0,0,0,0},
                                {0,0,0,0},
                                {0,0,0,0} };
            //rotate each vertex and normal vector
            for (int vertex = 0; vertex < 4; vertex++) { // col
                for (int coord = 0; coord < 3; coord++) { // row
                    for (int i = 0; i < 3; i++) {
                        newTriangleData[coord][vertex] += rotMat[coord][i]*oldTriangleData[i][vertex];
                    }
                }
            }
            //translate only the vertex coords, not the normal vector
            for (int vertex = 1; vertex < 4; vertex++) { // col
                for (int coord = 0; coord < 3; coord++) { // row
                    newTriangleData[coord][vertex] += translations[coord];
                }
            }
            m_FixedTriangleList[triangle].setTriangleData(newTriangleData);
            m_TriangleList[triangle].setTriangleData(newTriangleData);
        }
    }

    /*  Exports the Hull to a binary STL file.
        args:
        filePathName - the absolute file path and file name of the STL file to write to/create.
    */
    void Hull::exportBinaryStl(std::string filePathName) {
        int numTriangles = this->getNumTriangles();
        std::vector<std::vector<float>> triangleData;
        std::ofstream stl (filePathName, std::ios::out | std::ios::binary);
        //write the 80 char header of the STL file
        for (int i = 0; i < 80; i++) {
            stl.put((char) 0);
        }
        //write the number of triangles integer
        stl.write( reinterpret_cast<const char*>( & numTriangles), sizeof( int ));

        //write each triangle's data to the file
        for (int i = 0; i < numTriangles; i++) {
            triangleData = getTriangle(i).getTriangleData();
            //write the triangle normal and vertex data
            for (int vertex = 0; vertex < 4; vertex++) {
                for (int coord = 0; coord < 3; coord++) {
                    stl.write( reinterpret_cast<const char*>( &triangleData[coord][vertex]), sizeof( float ));
                }
            }
            //write the 'Attribute byte count' chars
            stl.put((char) 0);
            stl.put((char) 0);
        }
        stl.close();
    }

    /*  Read in a binary STL and return a Hull object representing the file 
        filePathName: a string describing the absolute filepath and filename of the file to be read in.
    */
    Hull readBinaryStl(std::string filePathName) {
        std::ifstream stl (filePathName, std::ios::in | std::ios::binary);
        //ignore the 80 byte STL file header
        stl.ignore(80);
        //read the number of triangles integer
        char charNumTriangles [4];
        stl.read(charNumTriangles, 4);
        int numTriangles = readInt(charNumTriangles,0);
        //read the triangle data. each triangle is stored in 50 bytes of data
        char * trianglesData = new char [numTriangles * 50];
        stl.read(trianglesData, numTriangles * 50);
        stl.close();
        //construct a Hull from triangleData and return
        return Hull(numTriangles, trianglesData);
    }

    /*  Round a float value to the specified number of decimal places/
        value: the value to be rounded.
        decimalPlaces: the number of decimal places to round to.
    */
    float round_up(float value, int decimalPlaces) {
        const float multiplier = std::pow(10.0, decimalPlaces);
        return std::roundf(value * multiplier) / multiplier;
    }

    /*  Converts 4 chars from an array of chars into an int.
        args:
        charBuffer - the char array to read from
        position - the position in charBuffer to use as the first char.
        returns:
        the 4 chars, converted to an integer value
    */
    int readInt(char* charBuffer, int position) {
        char tempBuffer [4] = {charBuffer[position+0], charBuffer[position+1], charBuffer[position+2], charBuffer[position+3]};
        return * (int*) tempBuffer;
    }

    /*  Converts 4 chars from an array of chars into a float.
        args:
        charBuffer - the char array to read from
        position - the position in charBuffer to use as the first char.
        returns:
        the 4 chars, converted to a float value
    */
    float readFloat(char* charBuffer, int position) {
        char tempBuffer [4] = {charBuffer[position+0], charBuffer[position+1], charBuffer[position+2], charBuffer[position+3]};
        return * (float*) tempBuffer;
    }

    /*  Converts a rotation described as a rotation around Z, Y and X into a rotation matrix.
        Angles are measured in degrees.
        args:
        rotX - the X component of the rotation.
        rotY - the Y component of the rotation.
        rotZ - the Z component of the rotation.
        returns:
        a 2D vector containing the reulting rotation matrix.
    */
    std::vector<std::vector<double>> makeRotMat(double rotX, double rotY, double rotZ) {
        //convert from degrees to radians
        double gamma = (rotX * M_PI) / 180;
        double beta = (rotY * M_PI) / 180;
        double alpha = (rotZ * M_PI) / 180;
        //construct rotation matrix
        std::vector<std::vector<double>> rotMat = {
            {   std::cos(beta)*std::cos(alpha),
                -std::cos(beta)*std::sin(alpha),
                std::sin(beta)
            },{
                std::cos(gamma)*std::sin(alpha) + std::sin(gamma)*std::sin(beta)*std::cos(alpha),
                std::cos(gamma)*std::cos(alpha) - std::sin(gamma)*std::sin(beta)*std::sin(alpha),
                -std::sin(gamma)*std::cos(beta)
            },{
                std::sin(gamma)*std::sin(alpha) - std::cos(gamma)*std::sin(beta)*std::cos(alpha),
                std::sin(gamma)*std::cos(alpha) + std::cos(gamma)*std::sin(beta)*std::sin(alpha),
                std::cos(gamma)*std::cos(beta)
            }
        };
        return rotMat;
    }
}
