#include "includes/utils/ant.h"

//Constructor
Ant::Ant(const Custom_CDT& initial_cdt): ant_cdt(initial_cdt), ant_energy(0.0), DeltaE(0.0), ant_conflict(false),
    ant_reduce_obtuses(false), num_of_obtuses(-1), ant_steiner_point(0,0), ant_steiner_method(NUM_METHODS),
    longest_edge(Point_2(0, 0), Point_2(0, 0)), opposite_edge(Point_2(0, 0), Point_2(0, 0)){}

Ant::Ant() :
    ant_cdt(Custom_CDT()),
    ant_steiner_method(NUM_METHODS),  //Member initializer list for initialization
    ant_energy(0.0),
    DeltaE(0.0),
    num_of_obtuses(-1),
    ant_conflict(false),
    ant_reduce_obtuses(false),
    longest_edge(Point_2(0, 0), Point_2(0, 0)),
    opposite_edge(Point_2(0, 0), Point_2(0, 0)) {
}

//Copy Constructor
Ant::Ant(const Ant& other): ant_cdt(other.ant_cdt), ant_energy(other.ant_energy), DeltaE(other.DeltaE), 
    ant_conflict(other.ant_conflict), ant_reduce_obtuses(other.ant_reduce_obtuses), 
    num_of_obtuses(other.num_of_obtuses), ant_steiner_point(other.ant_steiner_point), 
    ant_steiner_method(other.ant_steiner_method), ant_affect_faces(other.ant_affect_faces),
    longest_edge(other.longest_edge), opposite_edge(other.opposite_edge) {   
}

void Ant::set_Custom_CDT(Custom_CDT& new_cdt) {
    ant_cdt = new_cdt;
}

void Ant::set_face_in_ant_affect_faces(const Face_handle& face) {
    ant_affect_faces.insert(face);
}

void Ant::initialize_Ants(vector<Ant>& ants, Custom_CDT& best_cdt){
    int count_ants = ants.size();
    Point_2 temp_steiner_point(0,0);
    Segment_2 default_edge(Point_2(0, 0), Point_2(0, 0));
    for (int i = 0; i < count_ants; ++i) {
        ants[i] = Ant(best_cdt);
        ants[i].set_energy(0.0);
        ants[i].set_DeltaE(0.0);
        ants[i].set_reduce_obtuses(false);
        ants[i].set_steiner(temp_steiner_point);
        ants[i].set_num_of_obtuses(-1);
        ants[i].set_steiner_method(NUM_METHODS);
        ants[i].clear_ant_affect_faces();
        ants[i].set_conflict(false);
        ants[i].set_longest_edge_midpoint(default_edge);
        ants[i].set_opposite_edge_projection(default_edge);
    }
}

void Ant::clear_ant_affect_faces() {
    ant_affect_faces.clear();
}

//Set Steiner method
void Ant::set_steiner_method(SteinerMethod in_method) {
    ant_steiner_method = in_method;
}

void Ant::set_conflict(bool in_conflict){
    ant_conflict = in_conflict;
}

//Set Steiner point
void Ant::set_steiner(const Point_2& in_ant_steiner_point){
    ant_steiner_point = in_ant_steiner_point;
}

//Set DeltaE
void Ant::set_DeltaE(double in_DeltaE) {
    DeltaE = in_DeltaE;
}
//Set energy
void Ant::set_energy(double in_energy) {
    ant_energy = in_energy;
}
void Ant::set_reduce_obtuses(bool in_ant_reduce_obtuses){
    ant_reduce_obtuses = in_ant_reduce_obtuses;
}

void Ant::set_num_of_obtuses(const int in_num_of_obtuses){
    num_of_obtuses = in_num_of_obtuses;
}

void Ant::set_longest_edge_midpoint(Segment_2 in_longest_edge){
    longest_edge = in_longest_edge;
}

void Ant::set_opposite_edge_projection(Segment_2 in_opposite_edge){
    opposite_edge = in_opposite_edge;
}


///////////////Getters
std::set<Face_handle>& Ant::get_affected_faces(){
    return ant_affect_faces;
}

SteinerMethod Ant::get_steiner_method() const {
    return ant_steiner_method;
}

Custom_CDT& Ant::get_Custom_CDT(){
    return ant_cdt;
}

Segment_2 Ant::get_longest_edge_midpoint() const{
    return longest_edge;
}

Segment_2 Ant::get_opposite_edge_projection() const{
    return opposite_edge;
}

const Point_2& Ant::get_steiner_point() const {
    return ant_steiner_point;
}

double Ant::get_energy() const {
    return ant_energy;
}

double Ant::get_DeltaE() const {
    return DeltaE;
}

bool Ant::get_reduce_obtuses(){
    return ant_reduce_obtuses;
}

bool Ant::get_conflict() const{
    return ant_conflict;
}

int Ant::get_num_of_obtuses() const{
    return num_of_obtuses;
}