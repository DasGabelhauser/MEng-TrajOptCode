#include "RobotBuilder.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "coppeliaGeometricRoutines-master/geom.h"
#include <chrono>
using namespace Ipopt;





int main(int, char**) {
    //initialise a RobotBuilder
    //RobotBuilder is only used as this code was modified from main trajectory optimisation
    //code, in which RobotBuilder was used to load all necessary mesh files.
    EricsTrajOpt::RobotBuilder("3134");
    
    return(0);
}