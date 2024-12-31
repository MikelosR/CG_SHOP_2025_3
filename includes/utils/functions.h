#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "libraries.h"
#include "ant.h"

using namespace boost::json;
using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Polygon = CGAL::Polygon_2<K>;
using Point_2 = K::Point_2;
using Line_2 = K::Line_2;
using Segment_2 = K::Segment_2;
using Face_handle = Custom_CDT::Face_handle;
using Vertex_handle = Custom_CDT::Vertex_handle;
using std_string = std::string;
typedef K::FT FT;

//Steiner methods
//void insert_circumcenter_centroid(Custom_CDT& custom_cdt, const Face_handle& face, const Polygon& polygon, Point_2& circum_or_centroid);
void insert_projection(Custom_CDT& custom_cdt, const Face_handle& face, Polygon& polygon, Point_2& in_projection, Segment_2& opposide_edge);
void insert_midpoint(Custom_CDT& custom_cdt, const Face_handle& face, Polygon& polygon, Point_2& in_midpoint, Segment_2& longest_edge);
bool insert_adjacent_steiner(Custom_CDT& custom_cdt, const Face_handle& face, const Polygon& polygon, Point_2& adjacent_steiner);
void insert_adjacent_steiner_local_search(Custom_CDT& custom_cdt, const Face_handle& face1, const Polygon& polygon, Point_2& adjacent_steiner);
bool insert_circumcenter(Custom_CDT& circumcenter_cdt, const Face_handle& face, const Polygon& polygon, Point_2& circumcenter_steiner);
void insert_centroid(Custom_CDT& centroid_cdt, const Face_handle& face, const Polygon& polygon, Point_2& centroid_steiner);

//Helper function for circumcenter, checking if the circumcenter was placed in neighbor face
bool is_circumcenter_in_neighbor(const Custom_CDT& cdt, const Face_handle& face, const Point_2& circumcenter);
//Helper function for Midpoint (find the longest edge of a face)
Segment_2 find_longest_edge(const Face_handle& face);
//Helper function for Adjacent Steiner
bool has_obtuse_neighbors(const Custom_CDT& custom_cdt, const Face_handle& face, const Polygon& polygon);

//Algorithms
void local_search(Custom_CDT& custom_cdt, Polygon& polygon, int& L, const std_string& name_of_instance);
void simulated_annealing(Custom_CDT& custom_cdt, Polygon& polygon, int max_iterations, const double& alpha, const double& beta, const int& batch_size, const std_string& name_of_instance);
void ant_colony(Custom_CDT& custom_cdt, Polygon& polygon, const double& alpha, const double& beta, const double& chi, const double& psi, const double& lamda, const int& L, const int& kappa, const std_string& name_of_instance);

//Helper functions for Simulated Annealing
bool should_accept_bad_steiner(const double deltaE, const double T);
double calculate_energy(const int obtuse_faces, const int steiner_points, const double alpha, const double beta);

//Helper functions for Ant Colony
double calculate_radius_to_height(const Face_handle& face, const Custom_CDT& cdt);
double hta_vertex_projection(double rho);
double hta_circumcenter(double rho);
double hta_midpoint(double rho);
double hta_mean_adjacent(bool has_obtuse_neighbors);
Face_handle give_random_obtuse(Custom_CDT& custom_cdt, Polygon& polygon);
SteinerMethod selectSteinerMethod(const double& ro, const vector<double>& taf, vector<double>& hta, double chi, double psi, bool obtuse_neighbors);
void updatePheromones(vector<double>& taf, vector<double>& delta_taf, vector<Ant> selected_ants, double lamda);
bool are_faces_equal(const Face_handle& face1, const Face_handle& face2);
vector<Ant> save_the_best(vector<Ant>& ants);
//Check for conflict between 2 ants
bool have_conflict(Ant& ant1, Ant& ant2);
void printAntDetails(vector<Ant>& ants);
void affected_faces(Custom_CDT& best_cdt, Ant& ant);

/*General purpose functions*/
void start_the_flips(Custom_CDT& cdt, const Polygon& polygon);
//Return true if approves the flip
bool is_it_worth_flip(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);
void update_polygon(Polygon& polygon, const Point_2& steiner_point, const Point_2& point1, const Point_2& point2);
bool is_point_inside_region(const Point_2& point, const Polygon& polygon);
bool is_face_inside_region(const Face_handle& face, const Polygon& polygon);
bool is_edge_inside_region(const Point_2& point1, const Point_2& point2, const Polygon& polygon);
bool is_edge_on_boundary(const Point_2& p1, const Point_2& p2, const Polygon& polygon);
Point_2 compute_centroid(const vector<Point_2>& points);
int count_vertices(const Custom_CDT& cdt);
void print_polygon_edges(const Polygon& polygon);
bool is_polygon_convex(const set<Point_2>& unique_points);
//Check if a face is obtuse
bool is_obtuse(const Face_handle& face);
//Check if a face (3 points) is obtuse
bool is_obtuse2(const Point_2& p1, const Point_2& p2, const Point_2& p3);
//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles(CDT& cdt, const Polygon& polygon);
//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);
Point_2 find_obtuse_vertex(const Point_2& v1, const Point_2& v2, const Point_2& v3);
////////////////////////////////////////////////////////////////////////////////////////////////////////

//JSON INPUT - OUTPUT METHODS
void read_json(const std_string& filename, value& jv);
void output(value jv, Custom_CDT custom_cdt, vector<Point_2> points, int obtuse_count, std_string output_path);
bool is_steiner_point(Vertex_handle vertex, const vector<Point_2>& original_points);
std_string convert_to_string(const FT& coord);
std_string format_double(double value);

//3rd task
bool boundary_straight_lines(const Polygon& polygon);
bool are_constraints_closed(const vector<pair<int, int>>& additional_constraints, int num_points, const vector<Point_2>& points, 
                            const Polygon& polygon);
bool is_closed_from_boundary(const vector<pair<int, int>>& additional_constraints, const vector<Point_2>& points, 
                            const vector<int> degree, const Polygon& polygon);
bool is_closed_from_circle(const vector<pair<int, int>>& additional_constraints, const vector<Point_2>& points, 
                            const vector<int> degree, const Polygon& polygon);
void method_output(const vector<int> count_steiners, std_string method_name, const std_string& name_of_instance, 
                    const int num_steiners, const int num_obtuses);
bool is_same_edge(const Segment_2& seg1, const Segment_2& seg2);

vector<Face_handle> find_faces_intersecting_polygon_edges(const Custom_CDT& cdt, const Polygon& polygon);
vector<Face_handle> find_faces_inside_boundary(const Custom_CDT& cdt, const Polygon& polygon);
bool are_constraints_open(const vector<pair<int, int>>& additional_constraints, int num_points);
vector<Segment_2> find_shared_edges(CDT& cdt, const Polygon& polygon);
vector<Segment_2> find_non_touching_boundary_edges(const Polygon& polygon, const CDT& cdt, vector<Segment_2>& shared_edges);
vector<Segment_2> edges_new_boundary(const Polygon& polygon, vector<Segment_2>& shared_edges);


#endif