#include "includes/utils/functions_task1.h"
#include "includes/utils/Custom_Constrained_Delaunay_triangulation_2.h"

using namespace boost::json;
using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
typedef CGAL::Polygon_2<K> Polygon;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Point_2 = K::Point_2;
using Line_2 = K::Line_2;
using Face_handle = Custom_CDT::Face_handle;
using Vertex_handle = Custom_CDT::Vertex_handle;
using Segment_2 = K::Segment_2;

//Calculate squared distance between two points
double squared_distance_1(const Point_2 &a, const Point_2 &b)
{
    return CGAL::to_double(CGAL::squared_distance(a, b));
}

bool is_obtuse(const Point_2 &a, const Point_2 &b, const Point_2 &c)
{
    return (CGAL::angle(a, b, c) == CGAL::OBTUSE ||
            CGAL::angle(b, a, c) == CGAL::OBTUSE ||
            CGAL::angle(a, c, b) == CGAL::OBTUSE);
}

//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles_1(CDT &cdt, const Polygon &polygon)
{
    int obtuse_count = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point_2 p1 = fit->vertex(0)->point();
        Point_2 p2 = fit->vertex(1)->point();
        Point_2 p3 = fit->vertex(2)->point();

        if (is_obtuse(p1, p2, p3) &&
            polygon.bounded_side(CGAL::midpoint(p1, p2)) != CGAL::ON_UNBOUNDED_SIDE &&
            polygon.bounded_side(CGAL::midpoint(p1, p3)) != CGAL::ON_UNBOUNDED_SIDE &&
            polygon.bounded_side(CGAL::midpoint(p2, p3)) != CGAL::ON_UNBOUNDED_SIDE)
        {
            obtuse_count++;
        }
    }
    return obtuse_count;
}

// Return true if 2 faces (two triangles) form a convex polygon
bool is_convex_1(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3, const Point_2 &p4)
{
    vector<Point_2> points = {p1, p2, p3, p4};

    // Check if the polygon is convex
    return CGAL::is_convex_2(points.begin(), points.end(), K());
}

// Return true if approves the flip
bool can_flip(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3, const Point_2 &p4)
{
    // Check if the quadrilateral (p2, p4, p3, p1) is convex
    if (is_convex_1(p1, p2, p3, p4))
    {
        // Count obtuse angles before the flip
        int obtuse_before_cnt = 0;
        if (is_obtuse(p1, p2, p3))
            obtuse_before_cnt++;
        if (is_obtuse(p1, p3, p4))
            obtuse_before_cnt++;

        // Count obtuse angles after the potential flip
        int obtuse_after_cnt = 0;
        if (is_obtuse(p1, p2, p4))
            obtuse_after_cnt++;
        if (is_obtuse(p2, p3, p4))
            obtuse_after_cnt++;

        // Check collinearity of all possible triplets among p1, p2, p3, p4
        if (CGAL::orientation(p1, p2, p3) == CGAL::COLLINEAR ||
            CGAL::orientation(p1, p2, p4) == CGAL::COLLINEAR ||
            CGAL::orientation(p1, p3, p4) == CGAL::COLLINEAR ||
            CGAL::orientation(p2, p3, p4) == CGAL::COLLINEAR)
        {
            // Skip the flip if any three points are collinear
            return false;
        }

        // If its worth flipping, do it!
        if (obtuse_after_cnt < obtuse_before_cnt)
            return true;
        else
            return false;
    }
    // Case quadrilateral is not convex
    return false;
}

// Midpoint Insertion:
// Finds the longest edge of the obtuse triangle and calculates its midpoint.
// Simulates inserting the midpoint to see if it reduces obtuse angles. If successful, inserts it into the actual CDT.
void insert_midpoint_1(Custom_CDT &custom_cdt, Polygon &polygon)
{
    bool progress = true;
    while (progress)
    {
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face)
        {
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            if (is_obtuse(p1, p2, p3))
            {
                Custom_CDT simulation = custom_cdt;
                CGAL::Segment_2 longest_edge = find_longest_edge(p1, p2, p3);
                Point_2 midpoint = CGAL::midpoint(longest_edge.source(), longest_edge.target());

                if (is_point_inside_region(midpoint, polygon) && is_face_inside_region_1(face, polygon))
                {
                    int obtuses_before = count_obtuse_triangles_1(simulation, polygon);
                    simulation.insert_no_flip(midpoint);
                    start_the_flips_1(simulation, polygon);
                    if (obtuses_before > count_obtuse_triangles_1(simulation, polygon))
                    {
                        custom_cdt.insert_no_flip(midpoint);
                        start_the_flips_1(custom_cdt, polygon);
                        progress = true;
                        if((polygon.bounded_side(midpoint) == CGAL::ON_BOUNDARY))
                                update_polygon_1(polygon, midpoint, longest_edge.source(), longest_edge.target());
                        break;
                    }
                }
            }
            if (count_obtuse_triangles_1(custom_cdt, polygon) == 0)
                progress = 0;
        }
    }
}

void insert_orthocenter_1(Custom_CDT &custom_cdt, const Polygon &polygon)
{
    bool progress = true;
    while (progress)
    {
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face)
        {
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            //Check if the triangle is obtuse
            if (is_obtuse(p1, p2, p3))
            {
                Custom_CDT simulation = custom_cdt;
                Point_2 obtuse_vertex = find_obtuse_vertex_1(p1, p2, p3);

                //Calculate the orthocenter of the obtuse triangle
                Point_2 orthocenter = find_orthocenter(p1, p2, p3);

                //Verify if the orthocenter point can be inserted
                if (is_point_inside_region(orthocenter, polygon) && is_face_inside_region_1(face, polygon))
                {
                    int obtuses_before = count_obtuse_triangles_1(simulation, polygon);
                    simulation.insert_no_flip(orthocenter);
                    start_the_flips_1(simulation, polygon);

                    //Check if insertion of orthocenter reduces obtuse triangles
                    if (obtuses_before > count_obtuse_triangles_1(simulation, polygon))
                    {
                        custom_cdt.insert_no_flip(orthocenter);
                        start_the_flips_1(custom_cdt, polygon);
                        progress = true;
                        break;
                    }
                }
            }
            if (count_obtuse_triangles_1(custom_cdt, polygon) == 0)
                progress = 0;
        }
    }
}

void insert_projection_1(Custom_CDT &custom_cdt, Polygon& polygon)
{
    //Projection case
    //Progress flag to terminate the loop if after from 1 for loop without progress (no reduce the obtuses) to reminate the insert_projection
    bool progress = true;
    bool fill_boundary = true;
    while (progress)
    {
        //If we end with the loops for the boundary faces, set the progress of main case false
        if (!fill_boundary)
            progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face)
        {
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();
            Point_2 opposite1, opposite2;

            if (is_obtuse(p1, p2, p3))
            {
                Point_2 obtuse_angle_vertex = find_obtuse_vertex_1(p1, p2, p3);
                //Face_handle triangleA = face;
                //Find the obtuse point, and the 2 opposites
                if (obtuse_angle_vertex == p1)
                {
                    opposite1 = p2;
                    opposite2 = p3;
                }
                else if (obtuse_angle_vertex == p2)
                {
                    opposite1 = p1;
                    opposite2 = p3;
                }
                else
                {
                    opposite1 = p1;
                    opposite2 = p2;
                }
                Line_2 line(opposite1, opposite2);
                Point_2 projected_point = line.projection(obtuse_angle_vertex);
                //Copy the cdt for simulation
                Custom_CDT simulation = custom_cdt;
                bool insert_projection = is_point_inside_region(projected_point, polygon);

                int obtuses_before = count_obtuse_triangles_1(simulation, polygon);
                //Itarate all the faces of the boundary. If the point in side (not on) of the boundary is obtuse, insert projection
                if (fill_boundary && insert_projection && is_face_inside_region_1(face, polygon))
                {
                    if (is_face_on_boundary(custom_cdt, face) && !(polygon.bounded_side(obtuse_angle_vertex) == CGAL::ON_BOUNDARY))
                    {
                        simulation.insert_no_flip(projected_point);
                        start_the_flips_1(simulation, polygon);
                        if (obtuses_before > count_obtuse_triangles_1(simulation, polygon))
                        {
                            /*Original insertion of Projection*/
                            custom_cdt.insert_no_flip(projected_point);
                            start_the_flips_1(custom_cdt, polygon);
                            if((polygon.bounded_side(projected_point) == CGAL::ON_BOUNDARY))
                                update_polygon_1(polygon, projected_point, opposite2, opposite1);
                            
                        }
                    }
                }
                //main case, not only the boundary faces
                else if (insert_projection)
                {
                    simulation.insert(projected_point);
                    start_the_flips_1(simulation, polygon);

                    if (obtuses_before > count_obtuse_triangles_1(simulation, polygon))
                    {
                        /*Original insertion of Projection*/
                        custom_cdt.insert_no_flip(projected_point);
                        start_the_flips_1(custom_cdt, polygon);
                        progress = true;
                        if((polygon.bounded_side(projected_point) == CGAL::ON_BOUNDARY))
                                update_polygon_1(polygon, projected_point, opposite2, opposite1);
                        break;
                    }
                    // Case the projection point reduce -1 the obtuse angle, but make new obtuse angle in the neighbor
                    // If this neighbor is boundary face (at least 2 points lies in the boundary polygon) and the projection point
                    // is part of this face but is in (not on) the boundary, add the projection
                    else if (obtuses_before == count_obtuse_triangles_1(simulation, polygon))
                    {

                        for (auto sim_face = simulation.finite_faces_begin(); sim_face != simulation.finite_faces_end(); ++sim_face)
                        {
                            if (simulation.is_infinite(sim_face))
                                continue;

                            // Check if the projected point is part of a boundary face
                            if (sim_face->vertex(0)->point() == projected_point ||
                                sim_face->vertex(1)->point() == projected_point ||
                                sim_face->vertex(2)->point() == projected_point)
                            {
                                if (is_face_on_boundary(simulation, sim_face))
                                {
                                    if (is_obtuse(sim_face->vertex(0)->point(), sim_face->vertex(1)->point(), sim_face->vertex(2)->point()))
                                    {
                                        // Insert into the original triangulation
                                        custom_cdt.insert_no_flip(projected_point);
                                        start_the_flips_1(custom_cdt, polygon);
                                        progress = true;
                                        //"call" the fill boundary case
                                        fill_boundary = true;
                                        if((polygon.bounded_side(projected_point) == CGAL::ON_BOUNDARY))
                                            update_polygon_1(polygon, projected_point, opposite2, opposite1);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            //when the case of itaration of boundary faces ends
            if ((++face == custom_cdt.finite_faces_end()) && (fill_boundary))
            {
                fill_boundary = false;
            }
            --face;

            if (count_obtuse_triangles_1(custom_cdt, polygon) == 0)
            {
                progress = false;
                break;
            }
            continue;
        }
    }
}

//Function to check if an edge is part of the boundary of the polygon
bool is_edge_on_boundary_1(const Point_2 &p1, const Point_2 &p2, const Polygon &polygon)
{
    for (auto edge_it = polygon.edges_begin(); edge_it != polygon.edges_end(); ++edge_it)
    {
        if ((edge_it->source() == p1 && edge_it->target() == p2) ||
            (edge_it->source() == p2 && edge_it->target() == p1))
        {
            return true;
        }
    }
    return false;
}

//Ιnsert Steiner points at circumcenters
void insert_circumcenter_centroid_1(Custom_CDT &custom_cdt, const Polygon &polygon)
{
    bool progress = true;
    while (progress)
    {
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face)
        {
            //Get the vertices of the current triangle
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            if (is_obtuse(p1, p2, p3))
            {
                //Compute the circumcenter of the triangle
                Face_handle triangleA = face;
                Point circumcenter = CGAL::circumcenter(p1, p2, p3);
                if (is_point_inside_region(circumcenter, polygon))
                {
                    Custom_CDT simulation = custom_cdt;
                    Face_handle locate_face = simulation.locate(circumcenter);
                    if (is_convex_1(p1, p2, p3, circumcenter) && is_face_inside_region_1(locate_face, polygon))
                    {
                        int initial_obtuse_count = count_obtuse_triangles_1(custom_cdt, polygon);
                        /*Simulate circumcenter insertion*/
                        simulation.insert_no_flip(circumcenter);
                        start_the_flips_1(simulation, polygon);
                        int final_obtuse_count = count_obtuse_triangles_1(simulation, polygon);
                        //Check if the flip resolved obtuse angles in the two faces
                        if (final_obtuse_count < initial_obtuse_count)
                        {
                            /*Original circumcenter insertion*/
                            custom_cdt.insert_no_flip(circumcenter);
                            start_the_flips_1(custom_cdt, polygon);
                            progress = true;
                            final_obtuse_count = count_obtuse_triangles_1(custom_cdt, polygon);
                            break;
                        }
                    }
                }
                else
                {
                    //Cetroid case (if circumcenter went out of bounds)
                    Point_2 centroid = CGAL::centroid(p1, p2, p3);
                    bool insert_centroid = can_insert_centroid(custom_cdt, triangleA, centroid, polygon);
                    //Just check if the insertion of centroid has a benefit
                    if (insert_centroid && is_face_inside_region_1(face, polygon))
                    {
                        custom_cdt.insert(centroid);
                        start_the_flips_1(custom_cdt, polygon);
                        progress = true;
                        break;
                    }
                }
            }
        }
    }
}

//Just simulate if insert centroid (and auto flips) we reduce the obtuses
bool can_insert_centroid(Custom_CDT &custom_cdt, Face_handle &triangleA, const Point_2 &centroid, const Polygon &polygon)
{
    Custom_CDT simulation = custom_cdt;
    //Get the vertices of triangle A
    Point_2 p1 = triangleA->vertex(0)->point();
    Point_2 p2 = triangleA->vertex(1)->point();
    Point_2 p3 = triangleA->vertex(2)->point();

    int initial_obtuse_count = count_obtuse_triangles_1(custom_cdt, polygon);
    /*Simulate centroid insertion*/
    simulation.insert(centroid);
    start_the_flips_1(simulation, polygon);
    int final_obtuse_count = count_obtuse_triangles_1(simulation, polygon);

    //Check if the number of obtuse triangles decreased or stayed the same
    if (final_obtuse_count < initial_obtuse_count)
        return true;
    else
        return false;
}

void start_the_flips_1(Custom_CDT &cdt, const Polygon &polygon)
{
    bool progress = true;
    while (progress)
    {
        progress = false;
        for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge)
        {
            CDT::Face_handle f1 = edge->first;
            int i = edge->second;
            CDT::Face_handle f2 = f1->neighbor(i);

            if (cdt.is_infinite(f1) || cdt.is_infinite(f2))
                continue;

            Point_2 p1 = f1->vertex(cdt.ccw(i))->point(); // First vertex on the shared edge (Counter-Clock Wise)
            Point_2 p3 = f1->vertex(cdt.cw(i))->point();  // Second vertex on the shared edge (Clock Wise)
            Point_2 p2 = f1->vertex(i)->point();          // Opposite vertex in the first triangle
            //Mirror index gets the opposite vertex of the second triangle (f2)
            int mirror_index = cdt.mirror_index(f1, i);
            Point_2 p4 = f2->vertex(mirror_index)->point();
            //if the eddge is constraints or is on
            if (cdt.is_constrained(*edge) || is_edge_on_boundary_1(p1, p3, polygon))
                continue;
            
            if(!is_face_inside_region_1(f1, polygon) || !is_face_inside_region_1(f2, polygon)) continue;
            if (can_flip(p1, p2, p3, p4))
            {
                cdt.flip(f1, i);
                progress = true;
                break;
            }
            else
                continue;
        }
    }
}

//If 1 point is on the boundary
bool is_point_inside_region(const Point_2 &point, const Polygon polygon)
{

    //Check if the point is inside the polygon
    return (polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE) || (polygon.bounded_side(point) == CGAL::ON_BOUNDARY);
}

Point_2 find_obtuse_vertex_1(const Point_2 &v1, const Point_2 &v2, const Point_2 &v3)
{
    //Calculate squared distances
    double ab2 = CGAL::to_double(squared_distance_1(v1, v2));
    double ac2 = CGAL::to_double(squared_distance_1(v1, v3));
    double bc2 = CGAL::to_double(squared_distance_1(v2, v3));

    if (ab2 + ac2 < bc2)
        return v1; //obtuse angle at v1
    if (ab2 + bc2 < ac2)
        return v2; //obtuse angle at v2
    if (ac2 + bc2 < ab2)
        return v3; //obtuse angle at v3

    throw std::logic_error("No obtuse angle found in the triangle.");
}

CGAL::Segment_2<K> find_longest_edge(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    //Create segments for each edge of the triangle
    CGAL::Segment_2<K> edge1(p1, p2);
    CGAL::Segment_2<K> edge2(p2, p3);
    CGAL::Segment_2<K> edge3(p3, p1);

    //Compare the squared lengths of the segments to find the longest one
    auto length1 = CGAL::squared_distance(p1, p2);
    auto length2 = CGAL::squared_distance(p2, p3);
    auto length3 = CGAL::squared_distance(p3, p1);

    //Determine which edge is the longest
    if (length1 >= length2 && length1 >= length3)
    {
        return edge1;
    }
    else if (length2 >= length1 && length2 >= length3)
    {
        return edge2;
    }
    else
    {
        return edge3;
    }
}

Point_2 find_medial_of_longest_side(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    //Calculate the squared lengths of each side
    double length1_sq = CGAL::to_double(CGAL::squared_distance(p1, p2)); // Side p1 - p2
    double length2_sq = CGAL::to_double(CGAL::squared_distance(p2, p3)); // Side p2 - p3
    double length3_sq = CGAL::to_double(CGAL::squared_distance(p1, p3)); // Side p1 - p3

    //Find the longest side
    if (length1_sq >= length2_sq && length1_sq >= length3_sq)
    {
        //Longest side is p1 - p2
        return CGAL::midpoint(p1, p2);
    }
    else if (length2_sq >= length1_sq && length2_sq >= length3_sq)
    {
        //Longest side is p2 - p3
        return CGAL::midpoint(p2, p3);
    }
    else
    {
        //Longest side is p1 - p3
        return CGAL::midpoint(p1, p3);
    }
}

Point_2 find_orthocenter(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    //Calculate slopes and midpoints of sides
    Line_2 line1(p1, p2);
    Line_2 line2(p2, p3);
    Line_2 line3(p1, p3);

    //Find the orthogonal bisectors
    Line_2 perp1 = line1.perpendicular(p3);
    Line_2 perp2 = line2.perpendicular(p1);

    //Intersection of two perpendiculars gives the orthocenter
    CGAL::Object result = CGAL::intersection(perp1, perp2);
    const Point_2 *orthocenter = CGAL::object_cast<Point_2>(&result);

    //Ensure result validity
    if (orthocenter)
        return *orthocenter;
    else
        throw std::runtime_error("Orthocenter calculation failed.");
}

//Found if this face is on the boundary!!
bool is_face_on_boundary(const Custom_CDT &cdt, Face_handle face)
{
    //Loop through each edge of the face
    for (int i = 0; i < 3; ++i)
    {
        Custom_CDT::Edge edge(face, i);
        if (cdt.is_infinite(face->neighbor(i)))
        {
            return true; // Edge is on the boundary
        }
    }
    return false;
}

void run_task1(Custom_CDT& custom_cdt, Polygon& polygon){
    int init_obtuses = count_obtuse_triangles_1(custom_cdt, polygon);
    cout<<"Initial number of obtuses: "<<init_obtuses<<endl;
    int end;
    //If we have progress (reduce obtuses) run again
    bool progress = true;
    while(progress){
        progress = false;
        int start = count_obtuse_triangles_1(custom_cdt, polygon);
        int num_obtuses_before = start;
        //Flips
        start_the_flips_1(custom_cdt, polygon);
        int num_obtuses_after = count_obtuse_triangles_1(custom_cdt, polygon);
        if(num_obtuses_after == 0) break;

        //Circumcenter - Centroid
        num_obtuses_before = count_obtuse_triangles_1(custom_cdt, polygon);
        insert_circumcenter_centroid_1(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles_1(custom_cdt, polygon);
        if(num_obtuses_after == 0) break;

        //Midpoint
        num_obtuses_before = count_obtuse_triangles_1(custom_cdt, polygon);
        insert_midpoint_1(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles_1(custom_cdt, polygon);
        if(num_obtuses_after == 0) break;

        //Projection
        num_obtuses_before = count_obtuse_triangles_1(custom_cdt, polygon);
        insert_projection_1(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles_1(custom_cdt, polygon);
        if(num_obtuses_after == 0) break;
        
        //Orthocenter
        num_obtuses_before = count_obtuse_triangles_1(custom_cdt, polygon);
        insert_orthocenter_1(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles_1(custom_cdt, polygon);
        if(num_obtuses_after == 0) break;
        
        end = num_obtuses_after;
        if(end < start ) progress = true;
        else progress = false;
        if(end == 0) break;
    }
}

//If face is inside of region boundary
bool is_face_inside_region_1(const Face_handle& face, const Polygon& polygon) {
    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point(); 

    //Check if vertices are inside or on the boundary of the polygon
    bool vertices_inside = 
        (polygon.bounded_side(p1) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(p1) == CGAL::ON_BOUNDARY) &&
        (polygon.bounded_side(p2) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(p2) == CGAL::ON_BOUNDARY) &&
        (polygon.bounded_side(p3) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(p3) == CGAL::ON_BOUNDARY);

    if (!vertices_inside) return false;

    //Check if edges are fully inside or on the boundary
    bool edges_inside = 
        (polygon.bounded_side(CGAL::midpoint(p1, p2)) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(CGAL::midpoint(p1, p2)) == CGAL::ON_BOUNDARY) &&
        (polygon.bounded_side(CGAL::midpoint(p1, p3)) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(CGAL::midpoint(p1, p3)) == CGAL::ON_BOUNDARY) &&
        (polygon.bounded_side(CGAL::midpoint(p2, p3)) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(CGAL::midpoint(p2, p3)) == CGAL::ON_BOUNDARY);
    
    if (!edges_inside) return false;

    //The lines of the face
    Segment_2 line1(p1, p2), line2(p1, p3), line3(p2, p3);
    //bool edge_intersects_boundary = true;
    //For every edge of the polygon (boundary)
    for (auto edge = polygon.edges_begin(); edge != polygon.edges_end(); ++edge) {
        //if(!edge_intersects_boundary) break;
        Segment_2 poly_edge(*edge);
        //Check intersections for each face edge
        vector<Segment_2> face_edges = {line1, line2, line3};
        //Check every line of the face if it has intersection with the curent polygon edge
        for (const auto& face_edge : face_edges) {
            auto result = CGAL::intersection(face_edge, poly_edge);
            //If we have intersection between edge face and edge polygon (curent)
            if (result) {
                //Check if the intersection return point (not segment, line etc)
                if (const Point_2* point = boost::get<Point_2>(&*result)) {
                    //Check also for clean intersection, not intersection like the polygon edge intersect vertex of face
                    if ((*point != face_edge.source()) && (*point != face_edge.target()) && face_edge.has_on(*point)){
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void update_polygon_1(Polygon& polygon, const Point_2& steiner_point, const Point_2& p1, const Point_2& p2) {
    //Ensure the Steiner point is on the boundary of the polygon
    if (polygon.bounded_side(steiner_point) != CGAL::ON_BOUNDARY) {
        return; // The Steiner point is not on the polygon boundary, no update needed
    }

    //Iterate through the polygon's vertices to find the edge (p1, p2)
    auto it = polygon.vertices_begin();
    for (; it != polygon.vertices_end(); ++it) {
        auto next_it = next(it);
        if (next_it == polygon.vertices_end()) next_it = polygon.vertices_begin(); //Wrap around for closed polygon

        if ((*it == p1 && *next_it == p2) || (*it == p2 && *next_it == p1)) {
            //Insert the Steiner point between the vertices
            polygon.insert(next_it, steiner_point); //Insert before next_it
            return; //Exit after insertion
        }
    }
}