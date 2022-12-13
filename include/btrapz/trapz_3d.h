#pragma once

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <eigen3/Eigen/Dense>

class OccupancyGrid3D{
    public:
        enum CellStatus{
            FREE = 0,
            OCCUPIED = 1,
            UNKNOWN = 0
        };

        
}