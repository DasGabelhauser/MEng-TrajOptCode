#include "ChassisImportHandler.hpp"

#include <fstream>
#include <vector>
#include <string>

#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/Polyhedron_builder_from_STL.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace EricsTrajOpt {
    using namespace EricsTrajOpt;

    typedef Polyhedron::HalfedgeDS             HalfedgeDS;
    typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
    typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;

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
        return true;
    }

    /*  Read an STL file into a Hull object.
        args:
        filePathName - the absolute filepath and filename of the STL file.
        hull - a Hull object to return the mesh data in.
        returns:
        true if the operation succeeded.
        false if the operation did not succeed
    */
    bool readSTL2Hull(std::string filePathName, Hull * hull) {
        //read STL file into CGAL polyhedron
        Polyhedron P;
        if(readSTL2Polyhedron(filePathName, &P)){
            //create a new HUll object from the CGAl polyhedron
            *hull = *new Hull(P);
            return true;
        }
        return false;
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


    std::vector<Hull> importConvexDecomposedChassis(std::string filePathName, std::vector<blueprintLine> blueprint) {
        //read in polyhedron from file
        Polyhedron P;
        readSTL2Polyhedron(filePathName, &P);
        //convert to nef_polyhedron, as needed for convex decomposition
        Nef_polyhedron_3 nefPoly = Nef_polyhedron_3(P);
        //convex decompose nef_polyhedron
        CGAL::convex_decomposition_3(nefPoly);

        //return vector for covex sub-hulls
        std::vector<Hull> convexParts;
        Polyhedron subPoly;

        // the first volume is the outer volume, which is
        // ignored in the decomposition
        Volume_const_iterator ci = ++nefPoly.volumes_begin();

        //iterate through all convex sub parts
        for( ; ci != nefPoly.volumes_end(); ++ci) {
            if(ci->mark()) {
                //convert to polyhedron
                nefPoly.convert_inner_shell_to_polyhedron(ci->shells_begin(), subPoly);
                //convert to Hull and push to return vector
                Hull hull = *new Hull(subPoly);
                convexParts.push_back(hull);

            }
        }

        Hull headSlot;
        readSTL2Hull("/home/eric/Documents/MENG-Robot-Fab-Project/Optimisation Code/Models/Chassis/convexClips/head_slot_CH.stl", &headSlot);
        Hull organSlot;
        readSTL2Hull("/home/eric/Documents/MENG-Robot-Fab-Project/Optimisation Code/Models/Chassis/convexClips/slide_clip_CH.stl", &organSlot);
        headSlot.setHullPermPosRot( blueprint[0].posX, 
                                    blueprint[0].posY, 
                                    blueprint[0].posZ, 
                                    blueprint[0].rotX, 
                                    blueprint[0].rotY, 
                                    blueprint[0].rotZ);
        convexParts.push_back(headSlot);
        Hull tempOrganSlot;
        for (std::vector<blueprintLine>::size_type clip = 1; clip < blueprint.size(); clip++) {
            tempOrganSlot = *new Hull(organSlot);
            tempOrganSlot.setHullPermPosRot( blueprint[clip].posX, 
                                        blueprint[clip].posY, 
                                        blueprint[clip].posZ, 
                                        blueprint[clip].rotX, 
                                        blueprint[clip].rotY, 
                                        blueprint[clip].rotZ);
            convexParts.push_back(tempOrganSlot);
        }



        return convexParts;
    }

    Hull importConvexHullChassis(std::string filePathName) {
        //read in polyhedron from file
        Polyhedron P;
        readSTL2Polyhedron(filePathName, &P);

        Hull convexHull;
        Polyhedron convexPoly;

        //create convex hull
        CGAL::convex_hull_3(P.points_begin(), P.points_end(), convexPoly);

        //convert to hull
        convexHull = *new EricsTrajOpt::Hull(convexPoly);

        return convexHull;
    }


}