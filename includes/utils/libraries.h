//libraries.h
#ifndef LIBRARIES_H
#define LIBRARIES_H

#include "Custom_Constrained_Delaunay_triangulation_2.h"

//CGAL headers
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Line_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/intersection_2.h>

//Standard C++ libraries
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>
#include <cmath>
#include <optional>
#include <random>
#include <fstream>
#include <boost/json.hpp>
#include <algorithm>
#include <chrono>

//Other libraries
#include <stdio.h>


enum SteinerMethod {
    CIRCUMCENTER,
    MIDPOINT_LONGEST_EDGE,
    PROJECTION,
    ADJACENT_FACES,
    NUM_METHODS,
    CENTROID
};

enum InputCategory {
    CONVEX_NO_CONSTRAINTS,
    CONVEX_OPEN_CONSTRAINTS,
    CONVEX_CLOSED_CONSTRAINTS,
    NOT_CONVEX_PARALLEL_N0_CONSTRAINTS,
    UNSPECIFIED_BOUNDARY
};
#endif


