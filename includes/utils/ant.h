#ifndef ANT_H
#define ANT_H

#include "libraries.h"

using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Polygon = CGAL::Polygon_2<K>;
using Point_2 = K::Point_2;
using Face_handle = Custom_CDT::Face_handle;
using Segment_2 = K::Segment_2;



class Ant {
public:
    //Constructor without arguments
    Ant();  
    //With explicit accept and: Ant ant1 = cdt;
    explicit Ant(const Custom_CDT& initial_cdt);
    Ant(const Ant& other_ant);
    void set_steiner(const Point_2& in_ant_steiner_point);
    void set_face_in_ant_affect_faces(const Face_handle& face);
    void set_steiner_method(SteinerMethod in_method);
    void set_energy(double in_energy);
    void set_Custom_CDT(Custom_CDT& new_cdt);
    void set_DeltaE(double in_DeltaE);
    void set_conflict(bool in_conflict);
    //Static because we want to call it without an instance of Ant
    static void initialize_Ants(vector<Ant>& ants, Custom_CDT& best_cdt);
    void set_reduce_obtuses(bool in_ant_reduce_obtuses);
    void set_num_of_obtuses(const int in_num_of_obtuses);
    void set_longest_edge_midpoint(Segment_2 in_longest_edge);
    void set_opposite_edge_projection(Segment_2 in_opposite_edge);


    void clear_ant_affect_faces();
    std::set<Face_handle>& get_affected_faces();
    SteinerMethod get_steiner_method() const;
    Custom_CDT& get_Custom_CDT();
    const Point_2& get_steiner_point() const;
    bool get_reduce_obtuses();
    bool get_conflict() const;
    double get_energy() const;
    double get_DeltaE() const;
    int get_num_of_obtuses() const;
    Segment_2 get_longest_edge_midpoint() const;
    Segment_2 get_opposite_edge_projection() const;

private:
    std::set<Face_handle> ant_affect_faces;
    SteinerMethod ant_steiner_method;
    Custom_CDT ant_cdt;
    Point_2 ant_steiner_point;
    //Midpoint edge: We need this edge to check if the steiner was entered on the boundary
    Segment_2 longest_edge;
    //Projection edge: We need this edge to check if the steiner was entered on the boundary
    Segment_2 opposite_edge;
    double ant_energy;
    double DeltaE;
    bool ant_conflict;
    bool ant_reduce_obtuses;
    int num_of_obtuses;
};

#endif
