#include "ChassisImportHandler.hpp"

#include <fstream>
#include <vector>
#include <string>

#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/Polyhedron_builder_from_STL.h>
#include <CGAL/Point_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include "coppeliaGeometricRoutines-master/geom.h"

namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    typedef CGAL::Simple_cartesian<double>     Kernel;
    typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
    typedef Polyhedron::HalfedgeDS             HalfedgeDS;

    /*  Read an STL file into a CGAL Polyhedron object.
        args:
        filePathName - the absolute filepath and filename of the STL file.
        P            - a pointer to a Polyhedron object to return the mesh data in.
        returns:
        true if operation was successful.
        false if operation was unsuccessful.
    */
    bool readSTL2Polyhedron(std::string filePathName, Polyhedron * P) {
        std::ifstream stl(filePathName, std::ios::in | std::ios::binary);
        if(!stl) {
            return false;
        }
        CGAL::Polyhedron_builder_from_STL<HalfedgeDS> stlDelegate (stl);

        P->delegate(stlDelegate);

        stl.close();
        //CGAL::Polygon_mesh_processing::triangulate_faces(P);
        return true;
    }


    CObbStruct* readSTL2CoppMesh(std::string filePathName) {
        Polyhedron P;
        if(readSTL2Polyhedron(filePathName, &P)){
            std::vector<simReal> vertices;
            std::vector<int> indices;
            Polyhedron::Halfedge_around_facet_circulator facCirc;

            //iterate through each vertex
            for (Polyhedron::Vertex_iterator vertex = P.vertices_begin(); vertex != P.vertices_end(); ++vertex) {
                vertices.push_back(CGAL::to_double(vertex->point().x())/1000);
                vertices.push_back(CGAL::to_double(vertex->point().y())/1000);
                vertices.push_back(CGAL::to_double(vertex->point().z())/1000);
            }
            //iterate through each face
            for (Polyhedron::Facet_iterator facet = P.facets_begin(); facet != P.facets_end(); ++facet) {
                facCirc = facet->facet_begin();
                do {
                    indices.push_back((int)std::distance(P.vertices_begin(), facCirc->vertex()));
                } while ( ++facCirc != facet->facet_begin());
            }
            std::cout << "Vertices: " << vertices.size() << ", indices: " << indices.size() << "\n";
            return geom_createMesh(&vertices[0],vertices.size(),&indices[0],indices.size(), nullptr, 10000.0, 30);
        }
        return nullptr;
    }

    CObbStruct* readMultipleSTLs2CoppMesh(std::vector<std::string> filepathnames) {
        std::vector<simReal> vertices;
        std::vector<int> indices;
        Polyhedron::Halfedge_around_facet_circulator facCirc;
        
        for (int i = 0; i < (int)filepathnames.size(); i++) {
            Polyhedron P;
            readSTL2Polyhedron(filepathnames[i], &P);

            //iterate through each vertex
            for (Polyhedron::Vertex_iterator vertex = P.vertices_begin(); vertex != P.vertices_end(); ++vertex) {
                vertices.push_back(CGAL::to_double(vertex->point().x()));
                vertices.push_back(CGAL::to_double(vertex->point().y()));
                vertices.push_back(CGAL::to_double(vertex->point().z()));
            }
            //iterate through each face
            for (Polyhedron::Facet_iterator facet = P.facets_begin(); facet != P.facets_end(); ++facet) {
                facCirc = facet->facet_begin();
                do {
                    indices.push_back((int)std::distance(P.vertices_begin(), facCirc->vertex()));
                } while ( ++facCirc != facet->facet_begin());
            }
            std::cout << "Vertices: " << vertices.size() << ", indices: " << indices.size() << "\n";
        }
        std::cout << "Vertices: " << vertices.size() << ", indices: " << indices.size() << "\n";
        return geom_createMesh(&vertices[0],vertices.size(),&indices[0],indices.size());
    }

    std::vector<blueprintLine> readChassisBlueprint(std::string filePathName) {
        std::vector<blueprintLine> blueprint = {};
        std::ifstream blueprintFile(filePathName, std::ios::in);
        std::string line;
        std::string::size_type linePos;
        blueprintLine tempBPLine;

        std::getline(blueprintFile, line);
        while(!blueprintFile.eof()) {
            std::cout << line << "\n";
            tempBPLine.skeleton = std::stoi(line, &linePos);

            line = line.substr(linePos+1);
            tempBPLine.organType = std::stoi(line, &linePos);
            
            line = line.substr(linePos+1);
            tempBPLine.posX = (std::stod(line, &linePos)*1000);
            line = line.substr(linePos+1);
            tempBPLine.posY = (std::stod(line, &linePos)*1000);
            line = line.substr(linePos+1);
            tempBPLine.posZ = (std::stod(line, &linePos)*1000);
            line = line.substr(linePos+1);
            tempBPLine.rotX = (std::stod(line, &linePos)*180)/M_PI;
            line = line.substr(linePos+1);
            tempBPLine.rotY = (std::stod(line, &linePos)*180)/M_PI;
            line = line.substr(linePos+1);
            tempBPLine.rotZ = (std::stod(line, &linePos)*180)/M_PI;

            blueprint.push_back(tempBPLine);

            std::getline(blueprintFile, line);
        }
        return blueprint;
    }

}