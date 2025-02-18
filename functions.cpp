#include "includes/utils/functions.h"
#include "includes/utils/extra_graphics.h"

using namespace boost::json;
using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
typedef CGAL::Polygon_2<K> Polygon;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Point_2 = K::Point_2;
using Line_2 = K::Line_2;
using Segment_2 = K::Segment_2;
using Face_handle = Custom_CDT::Face_handle;
using Vertex_handle = Custom_CDT::Vertex_handle;
using std_string = std::string;

//Is obtuse face
bool is_obtuse(const Face_handle& face) {
    Point_2 a = face->vertex(0)->point();
    Point_2 b = face->vertex(1)->point();
    Point_2 c = face->vertex(2)->point();
    return (CGAL::angle(a, b, c) == CGAL::OBTUSE || 
            CGAL::angle(b, a, c) == CGAL::OBTUSE || 
            CGAL::angle(a, c, b) == CGAL::OBTUSE);
}

//Is obtuse 3 points (1 face)
bool is_obtuse2(const Point_2& a, const Point_2& b, const Point_2& c) {
    return (CGAL::angle(a, b, c) == CGAL::OBTUSE || 
            CGAL::angle(b, a, c) == CGAL::OBTUSE || 
            CGAL::angle(a, c, b) == CGAL::OBTUSE);
}

//Read JSON file
void read_json(const std_string& filename, value& jv) {
    ifstream file(filename);
    if (!file) {
        cerr<<"Error opening file: "<<filename<<endl;
        return;
    }
    std_string json_str((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
    jv = parse(json_str);
}

//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles(CDT& cdt, const Polygon& polygon) {
    int obtuse_count = 0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        Point_2 p1 = face->vertex(0)->point();
        Point_2 p2 = face->vertex(1)->point();
        Point_2 p3 = face->vertex(2)->point();
        if (is_obtuse(face) && is_face_inside_region(face, polygon)) obtuse_count++;
    }
    return obtuse_count;
}


//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4) {
    vector<Point_2> points = {p1, p2, p3, p4};

    //Check if the polygon is convex
    return CGAL::is_convex_2(points.begin(), points.end(), K()); 
}

//Return true if approves the flip
bool is_it_worth_flip(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4){
    //Check if the quadrilateral (p2, p4, p3, p1) is convex
    if (is_convex(p1, p2, p3, p4)) {
            //Count obtuse angles before the flip
            int obtuse_before_cnt = 0;
            if(is_obtuse2(p1, p2, p3)) obtuse_before_cnt++;
            if(is_obtuse2(p1, p3, p4)) obtuse_before_cnt++;

            //Count obtuse angles after the potential flip
            int obtuse_after_cnt = 0;
            if(is_obtuse2(p1, p2, p4)) obtuse_after_cnt++;
            if(is_obtuse2(p2, p3, p4)) obtuse_after_cnt++;

            //Check collinearity of all possible triplets among p1, p2, p3, p4
            if (CGAL::orientation(p1, p2, p3) == CGAL::COLLINEAR ||
                CGAL::orientation(p1, p2, p4) == CGAL::COLLINEAR ||
                CGAL::orientation(p1, p3, p4) == CGAL::COLLINEAR ||
                CGAL::orientation(p2, p3, p4) == CGAL::COLLINEAR) {
                //Skip the flip if any three points are collinear
                return false;
            }

            //If its worth flipping, do it!
            if(obtuse_after_cnt < obtuse_before_cnt) return true;
            else return false;
    }
    //Case quadrilateral is not convex
    return false;      
}


//Helper function to compute the centroid of a set of points (polygon)
Point_2 compute_centroid(const vector<Point_2>& points) {
    //Initialize bounding box values with the first point, converted to double
    double x_min = CGAL::to_double(points[0].x());
    double x_max = x_min;
    double y_min = CGAL::to_double(points[0].y());
    double y_max = y_min;

    //Loop through each point, updating bounding box coordinates
    for (const auto& p : points) {
        x_min = min(x_min, CGAL::to_double(p.x()));
        x_max = max(x_max, CGAL::to_double(p.x()));
        y_min = min(y_min, CGAL::to_double(p.y()));
        y_max = max(y_max, CGAL::to_double(p.y()));
    }

    //Calculate the centroid as the midpoint of the bounding box
    return Point_2((x_min + x_max) / 2, (y_min + y_max) / 2);
}

//Projection case
void insert_projection(Custom_CDT& custom_cdt, const Face_handle& face, Polygon& polygon, Point_2& in_projection, Segment_2& opposide_edge){
    
    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point();
    Point_2 opposite1, opposite2;

    Point_2 obtuse_angle_vertex = find_obtuse_vertex(p1, p2, p3);
    //Find the obtuse point, and the 2 opposites
    if (obtuse_angle_vertex == p1) {
        opposite1 = p2;
        opposite2 = p3;
    } else if (obtuse_angle_vertex == p2) {
        opposite1 = p1;
        opposite2 = p3;
    } else if (obtuse_angle_vertex == p3){
        opposite1 = p1;
        opposite2 = p2;
    }
    else{
        cout<<"There is no obtuse angle"<<endl;
    }
    
    //Take the projection Steiner point
    Line_2 line(opposite1, opposite2);
    Point_2 projected_point = line.projection(obtuse_angle_vertex);
    //"Return" and the projected_point
    in_projection = projected_point;
    //And the opposide edge of the projected_point
    opposide_edge = Segment_2(opposite1, opposite2);
    bool insert_projection = is_point_inside_region(projected_point, polygon);
    int obtuses_before = count_obtuse_triangles(custom_cdt, polygon);

    //standard check
    if(insert_projection){
        custom_cdt.insert_no_flip(projected_point);
        start_the_flips(custom_cdt, polygon);
    }
    else {
        cout<<"PROJECTION DIDN'T INSERTED"<<endl;
        cout<<"projected_point.x: "<<projected_point.x()<<" projected_point.y: "<<projected_point.y()<<endl;
    }
}

//Midpoint Insertion:
//Finds the longest edge of the obtuse triangle and calculates its midpoint.
void insert_midpoint(Custom_CDT& custom_cdt, const Face_handle& face, Polygon& polygon, Point_2& in_midpoint, Segment_2& longest_edge) {
    longest_edge = find_longest_edge(face);
    Point_2 midpoint = CGAL::midpoint(longest_edge.source(), longest_edge.target());
    in_midpoint = midpoint;

    if (is_point_inside_region(midpoint, polygon)) {
        custom_cdt.insert_no_flip(midpoint);
        start_the_flips(custom_cdt, polygon);
    }  
}

bool insert_adjacent_steiner(Custom_CDT& custom_cdt, const Face_handle& face1, const Polygon& polygon, Point_2& adjacent_steiner) {
    //Before calling insert_adjacent_steiner, we know that the face1 is obtuse face
    if (!has_obtuse_neighbors(custom_cdt, face1, polygon)) return false;
    set<Point_2> unique_points;        //Collect all unique vertices of obtuse neighbors
    queue<Face_handle> face_queue;     //For BFS traversal
    set<Face_handle> visited_faces;    //To avoid revisiting faces

    //Start BFS with the given face
    face_queue.push(face1);
    visited_faces.insert(face1);

    //Add the vertices of the initial face to the unique points
    for (int i = 0; i < 3; ++i) {
        unique_points.insert(face1->vertex(i)->point());
    }

    //While we don't have anymore face handles in the queue
    while (!face_queue.empty()) {
        Face_handle curent_face = face_queue.front();
        face_queue.pop();

        //Check if the current face is inside the polygon
        if (!is_face_inside_region(curent_face, polygon)) continue;
        
        //For every edge of the curent face, find obtuse neighbors and insert them into face_queue
        for (int i = 0; i < 3; ++i) {
            Face_handle neighbor = curent_face->neighbor(i);

            //Get one edge of the curent_face
            Point_2 p1 = curent_face->vertex((i + 1) % 3)->point();
            Point_2 p2 = curent_face->vertex((i + 2) % 3)->point();

            //Skip already visited faces or unsuitable neighbors
            if (visited_faces.count(neighbor) > 0 || custom_cdt.is_infinite(neighbor) || 
                !is_obtuse(neighbor) || custom_cdt.is_constrained(make_pair(curent_face, i)) || is_edge_on_boundary(p1, p2, polygon)) {
                continue;
            }

            //Check if the neighbor is obtuse
            if (is_obtuse(neighbor)) {
                //Add all vertices of this obtuse neighbor to the unique points
                for (int j = 0; j < 3; ++j) {
                    unique_points.insert(neighbor->vertex(j)->point());
                }
                //Mark this neighbor as visited and add it to the queue
                visited_faces.insert(neighbor);
                face_queue.push(neighbor);
            }
            else continue;
        }
    }

    //Compute the midpoint of all collected vertices
    if (!unique_points.empty()) {
        vector<Point_2> points(unique_points.begin(), unique_points.end());
        adjacent_steiner = compute_centroid(points);
    } else {
        //Default to the centroid of the original face if no neighbors are found
        Point_2 v0 = face1->vertex(0)->point();
        Point_2 v1 = face1->vertex(1)->point();
        Point_2 v2 = face1->vertex(2)->point();
        adjacent_steiner = CGAL::centroid(v0, v1, v2);
    }
    //Check if the polygon is convex
    if(is_polygon_convex(unique_points)){
        custom_cdt.insert_no_flip(adjacent_steiner);
        start_the_flips(custom_cdt, polygon);
        return true;
    }
    return false;
    
}

//Adhjacent steiner method only for local search
void insert_adjacent_steiner_local_search(Custom_CDT& custom_cdt, const Face_handle& face1, const Polygon& polygon, Point_2& adjacent_steiner) {
    set<Point_2> unique_points;
    unsigned int initial_obtuse_count = count_obtuse_triangles(custom_cdt, polygon);
    unsigned int best_obtuse_count = initial_obtuse_count;
    //Initialize best_steiner_point as the centroid of face1
    Point_2 v0 = face1->vertex(0)->point();
    Point_2 v1 = face1->vertex(1)->point();
    Point_2 v2 = face1->vertex(2)->point();
    Point_2 best_steiner_point = CGAL::centroid(v0, v1, v2);
    adjacent_steiner = best_steiner_point;
    //Use a queue to perform BFS-like traversal
    queue<Face_handle> face_queue;
    set<Face_handle> visited_faces;
    //Start with the initial face
    face_queue.push(face1);
    visited_faces.insert(face1);
    //Add the vertices of the curent face to the unique points
    for (int i = 0; i < 3; ++i) {
        unique_points.insert(face1->vertex(i)->point());
    }

    while (!face_queue.empty()) {
        Face_handle curent_face = face_queue.front();
        face_queue.pop();
        
        //Check if the curent face is inside the region boundary
        if (!is_face_inside_region(curent_face, polygon))  continue;

        //Traverse the neighbors of the curent face
        for (int i = 0; i < 3; ++i) {
            Face_handle neighbor = curent_face->neighbor(i);

            //Skip already visited faces or unsuitable neighbors
            if (visited_faces.count(neighbor) > 0 || custom_cdt.is_infinite(neighbor) || 
                !is_obtuse(neighbor) || custom_cdt.is_constrained(make_pair(curent_face, i))) {
                continue;
            }
            //We have just an iterator and this iterator may go out of bounds
            //because he doesn't have the supervision of for loop
            if(!is_face_inside_region(neighbor, polygon)) continue;
            
            //Add the vertices of the curent face to the unique points
            for (int i = 0; i < 3; ++i) {
                unique_points.insert(neighbor->vertex(i)->point());
            }
            //Compute the centroid of all collected points as a Steiner candidate
            vector<Point_2> polygon_points(unique_points.begin(), unique_points.end());
            Point_2 curent_steiner_point = compute_centroid(polygon_points);

            //Simulate inserting this Steiner point in a temporary CDT
            Custom_CDT simulate_cdt = custom_cdt;
            simulate_cdt.insert_no_flip(curent_steiner_point);
            start_the_flips(simulate_cdt, polygon);
            unsigned int simulated_obtuse_count = count_obtuse_triangles(simulate_cdt, polygon);
        
            //Change the best_steiner_point and update the best_obtuse_count if worth it
            if((simulated_obtuse_count < best_obtuse_count) && is_polygon_convex(unique_points)){
                //Update best_obtuse_count, best_steiner_point, adjacent_steiner
                best_obtuse_count = simulated_obtuse_count;
                best_steiner_point = curent_steiner_point;
                adjacent_steiner = best_steiner_point;
                custom_cdt.insert_no_flip(best_steiner_point);
                start_the_flips(custom_cdt, polygon);
            }
            //Mark the neighbor as visited and add it to the queue
            visited_faces.insert(neighbor);
            face_queue.push(neighbor);
        }
    }
}

//The local Search method
void local_search(Custom_CDT& custom_cdt, Polygon& polygon, int& L, const std_string& name_of_instance, 
                bool& in_randomization, const double& alpha, const double& beta, vector<int> subset, std_string category,
                const bool& run_auto_method){
    vector<int> count_steiners(6, 0);
    vector<Point_2> random_steiners;
    Point_2 temp_random_steiner;
    unsigned int num_of_obtuses = 0, num_of_steiners = 0, obtuse_custom = 0, num_of_obtuses_before = 0;
    bool progress = true, try_randomization = false;
    int obtuse_best_cdt = 0;
    //3rd task
    int init_vertices = custom_cdt.number_of_vertices();
    int init_num_obtuses = count_obtuse_triangles(custom_cdt, polygon);
    double p_sum = 0.0;
    time_t start_time, end_time; 
    time(&start_time);
    Custom_CDT best_cdt = custom_cdt;

    while(L > 0){
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
            if (!is_obtuse(face)) continue;
            if (!is_face_inside_region(face, polygon)) continue;
            
            //Create temporary CDT copies for each method
            Custom_CDT cdt_circum = custom_cdt;
            Custom_CDT cdt_midpoint = custom_cdt;
            Custom_CDT cdt_projection = custom_cdt;
            Custom_CDT cdt_adjacent = custom_cdt;
            Custom_CDT cdt_centroid = custom_cdt;
            vector<Custom_CDT> cdt_variants = {cdt_circum, cdt_midpoint, cdt_projection, cdt_adjacent, 
                                                cdt_centroid};
            vector<Point_2> steiner_points(6);
            //Midpoint edge: We need this edge to check if the steiner was entered on the boundary
            Segment_2 longest_edge;
            //Projection edge: We need this edge to check if the steiner was entered on the boundary
            Segment_2 opposide_edge;
            //Apply Steiner point insertion methods
            insert_circumcenter(cdt_variants[0], face, polygon, steiner_points[0]);// dont_use_circumcenter = true;
            insert_midpoint(cdt_variants[1], face, polygon, steiner_points[1], longest_edge);
            insert_projection(cdt_variants[2], face, polygon, steiner_points[2], opposide_edge);
            insert_adjacent_steiner_local_search(cdt_variants[3], face, polygon, steiner_points[3]);
            insert_centroid(cdt_variants[4], face, polygon, steiner_points[4]);

            //Vector to store obtuse counts
            vector<unsigned int> obtuses_after(5);

            for (int i = 0; i < cdt_variants.size(); ++i) {
                obtuses_after[i] = count_obtuse_triangles(cdt_variants[i], polygon);
            }
            //Find the method with the minimum obtuse triangles
            auto min_iter = std::min_element(obtuses_after.begin(), obtuses_after.end());
            unsigned int min_index = std::distance(obtuses_after.begin(), min_iter);
            obtuse_best_cdt = count_obtuse_triangles(best_cdt, polygon);
            //Apply the best method
            if (obtuses_after[min_index] < obtuse_best_cdt) {
                num_of_obtuses_before = count_obtuse_triangles(custom_cdt, polygon);
                custom_cdt.insert_no_flip(steiner_points[min_index]);
                start_the_flips(custom_cdt, polygon);
                custom_cdt = cdt_variants[min_index];
                best_cdt = custom_cdt;
                obtuse_best_cdt = count_obtuse_triangles(best_cdt, polygon);
                if(run_auto_method){                
                    //3rd task
                    if(try_randomization){
                        in_randomization = true;
                        cout<<"Random steiner inserted: "<<temp_random_steiner<<endl;
                        random_steiners.emplace_back(temp_random_steiner);
                        count_steiners[5]++;
                        try_randomization = false;
                    }
                    else progress = true;
                }

                count_steiners[min_index]++;
                //3rd task
                if(run_auto_method){
                    num_of_steiners = custom_cdt.number_of_vertices() - init_vertices;
                    p_sum += p_sum_function(num_of_steiners - 1, num_of_obtuses_before, obtuses_after[min_index]);
                }

                //For projection or midpoint check if the steiner inserted in the boundary of polygon and update the polygon
                if(min_index == 1) update_polygon(polygon, steiner_points[min_index], longest_edge.source(), longest_edge.target());
                if(min_index == 2) update_polygon(polygon, steiner_points[min_index], opposide_edge.source(), opposide_edge.target());
                break; //Restart iteration
            }
        }
        
        if(!progress){
            L--;
            custom_cdt = best_cdt;
            if(run_auto_method){
                try_steiner_around_centroid(custom_cdt, polygon, temp_random_steiner);
                obtuse_custom = count_obtuse_triangles(custom_cdt, polygon);
                obtuse_best_cdt = count_obtuse_triangles(best_cdt, polygon);
                
                try_randomization = true;
                if(obtuse_custom < obtuse_best_cdt) {
                    in_randomization = true;
                    best_cdt = custom_cdt;
                    count_steiners[5]++;
                    random_steiners.emplace_back(temp_random_steiner);
                    try_randomization = false;
                    cout<<"Random steiner inserted: "<<temp_random_steiner<<endl;
                }
                //Try this new cdt
            }
        }
    }

    custom_cdt = best_cdt;
    if(run_auto_method){
        double front;
        num_of_steiners = best_cdt.number_of_vertices() - init_vertices;
        if (num_of_steiners > 1)
            front = abs(1.0/(num_of_steiners - 1.0));
        else front = 0.0;

        double rate_of_convergence = front * p_sum;
        obtuse_best_cdt = count_obtuse_triangles(best_cdt, polygon);
        double Energy = calculate_energy(obtuse_best_cdt, num_of_steiners, alpha, beta);
        std_string method_name = "Local Search";
        int num_of_steiner = 0;

        //Sum the total number of steiners
        for(int i = 0; i < count_steiners.size(); ++i) {
            num_of_steiner += count_steiners[i];
        }

        num_of_obtuses = count_obtuse_triangles(custom_cdt, polygon);
        method_output(count_steiners, method_name, name_of_instance, num_of_steiner, init_num_obtuses, num_of_obtuses, 
                    in_randomization, random_steiners, rate_of_convergence, Energy, subset, category);
    }  
    time(&end_time);
    double time_taken = double(end_time - start_time); 
    cout<<"Time taken by program: "<<name_of_instance<<" is : "<<time_taken<<" sec "<<endl;
}


double calculate_energy(const int obtuse_faces, const int steiner_points, const double alpha, const double beta) {
    return alpha * obtuse_faces + beta * steiner_points;
}

bool should_accept_bad_steiner(const double deltaE, const double T) {
    //Compute e^(-∆E / T)
    double probability = exp(-deltaE / T);
    
    //Generate R uniformly in [0, 1]
    double R = static_cast<double>(rand()) / RAND_MAX;

    //Accept transition if e^(-∆E / T) ≥ R
    return probability >= R;
}

//Return p(n) with given that n_steiner, n_steiner + 1, previous_obtuses, current_obtuses
double p_sum_function(int n_steiner, int previous_obtuses, int obtuse_faces){
    //Ensure we do not encounter log(0) or log of negative values
    if (previous_obtuses <= 0 || obtuse_faces <= 0 || n_steiner <= 0) return 0;
    //ln(1) == 0
    if (obtuse_faces == previous_obtuses)  return 0;

    //Calculate p(n)
    double p_n = abs((log(static_cast<double>(obtuse_faces) / previous_obtuses))) / abs((log(static_cast<double>(n_steiner + 1) / n_steiner)));
    return p_n;
}


//Simualated annealing method
void simulated_annealing(Custom_CDT& custom_cdt, Polygon& polygon, int max_iterations, const double& alpha, 
                        const double& beta, const int& batch_size, const std_string& name_of_instance, 
                        bool& randomization, vector<int> subset, std_string category, const bool& run_auto_method){
    int obtuse_faces = count_obtuse_triangles(custom_cdt, polygon);
    
    int init_vertices = custom_cdt.number_of_vertices();
    double T = 1.0, delta_E = 0.0, E_new = 0.0, cooling_rate = 0.99, min_temp = 1e-6;
    double best_E = calculate_energy(obtuse_faces, 0, alpha, beta);
    int num_of_transition = 0, random_steiner = 0;
    vector<int> count_steiners(6, 0), temp_counter_steiner(6,0);
    Custom_CDT simulate_cdt = custom_cdt, best_cdt = custom_cdt, curent_cdt = custom_cdt;
   
    //These edges are Midpoint and opposite Projection edges, we need these edges to check, if these steiners was entered on the boundary
    Segment_2 longest_edge, opposite_edge;
    int best_num_steiner = 0, best_obtuse_faces = obtuse_faces, counter_steiner = 0;
    bool obtuse_neighbors = false, is_polygon_convex = false, try_randomization = false;

    //From the random number generator choose method from values vector
    std::mt19937 rng(std::random_device{}()); // Initialize RNG
    vector<int> values = subset; // Define possible values
    std::uniform_int_distribution<int> dist(0, values.size() - 1);

    //3rd task
    int init_num_obtuses = obtuse_faces;
    double p_sum = 0.0, p_sum_best = 0.0, temp_p_sum = 0;
    int previous_obtuses = obtuse_faces, previous_num_steiner = 0;

    fill(temp_counter_steiner.begin(), temp_counter_steiner.end(), 0);
    time_t start_time, end_time; 
    time(&start_time);
    vector<Point_2> vector_random_steiners;
    Point_2 temp_random_steiner, steiner_point;

    for (int i = 0; i < max_iterations && T > min_temp; ++i) {
        if (obtuse_faces == 0) break;
        //After from brake, or finite_faces_end(), keep the simulate_cdt as curent_cdt
        curent_cdt = simulate_cdt;
        for (auto face = curent_cdt.finite_faces_begin(); face != curent_cdt.finite_faces_end(); ++face){
            if (!is_obtuse(face)) continue;
            if (!is_face_inside_region(face, polygon)) continue;
            //Choose a random steiner from vector
            random_steiner = values[dist(rng)];
            switch(random_steiner){
                //If circumcenter steiner is outside of the boundary, continue
                case 0: 
                    if(!insert_circumcenter(simulate_cdt, face, polygon, steiner_point)) continue;
                    
                    break;
                case 1: insert_midpoint(simulate_cdt, face, polygon, steiner_point, longest_edge); break;
                case 2: insert_projection(simulate_cdt, face, polygon, steiner_point, opposite_edge); break;
                case 3:
                    //If the polygon of the adjacent steiner is not convex or if the face has no obtuse neighbors, skip the face
                    is_polygon_convex = insert_adjacent_steiner(simulate_cdt, face, polygon, steiner_point);
                    if((!is_polygon_convex)){
                        insert_projection(simulate_cdt, face, polygon, steiner_point, opposite_edge);
                        random_steiner = PROJECTION;
                    }
                    break;
                case 4: insert_centroid(simulate_cdt, face, polygon, steiner_point); break;
                default: break;
            }
            obtuse_faces = count_obtuse_triangles(simulate_cdt, polygon);
            counter_steiner = simulate_cdt.number_of_vertices() - init_vertices;
            E_new = calculate_energy(obtuse_faces, counter_steiner, alpha, beta);
            delta_E = E_new - best_E;

            //For any undetectable program error
            if (delta_E == 0) {
                simulate_cdt = curent_cdt;
                //cout<<"EROOR delta_E == 0"<<endl;
                continue;
            }
            //Trick to insert into should_accept_bad_steiner(delta_E,T) to reintroduce triangulation as best_cdt because we have increase the obtuses by 3
            if (delta_E >= (3*alpha)) delta_E = 0.000001;
            
            if(delta_E < 0){
                //Update the iterator
                curent_cdt = simulate_cdt;
                //Update the best value
                best_cdt = simulate_cdt;
                best_E = E_new;

                //Optional for prints
                best_obtuse_faces = obtuse_faces;
                //3rd task
                if(run_auto_method){                
                    p_sum += p_sum_function(counter_steiner - 1, previous_obtuses, obtuse_faces) + temp_p_sum;
                    previous_obtuses = obtuse_faces;
                    temp_p_sum = 0;
                    p_sum_best = p_sum;
                }
                
                //Restart the counter of bad steiner insertions
                num_of_transition = 0;
                //Update the counters for steiners
                temp_counter_steiner[random_steiner]++;
                //3rd task for the output stats file
                for(int i = 0; i < temp_counter_steiner.size(); ++i) {
                    count_steiners[i] += temp_counter_steiner[i];
                }

                if(try_randomization && run_auto_method){
                    cout<<"Random steiner inserted: "<<temp_random_steiner<<endl;
                    vector_random_steiners.emplace_back(temp_random_steiner);
                    randomization = true;
                    try_randomization = false;
                }
               
                //For projection or midpoint check if the steiner inserted in the boundary of polygon and update the polygon
                if(random_steiner == 1) update_polygon(polygon, steiner_point, longest_edge.source(), longest_edge.target());
                if(random_steiner == 2) update_polygon(polygon, steiner_point, opposite_edge.source(), opposite_edge.target());
                fill(temp_counter_steiner.begin(), temp_counter_steiner.end(), 0);
                break;
            }
            else if(should_accept_bad_steiner(delta_E,T)){
                num_of_transition++;
                //3rd task
                if(run_auto_method){
                    temp_p_sum += p_sum_function(counter_steiner - 1, previous_obtuses, obtuse_faces);
                    previous_obtuses = obtuse_faces;
                }
                
                //Run the for loop with the simulated_cdt
                curent_cdt = simulate_cdt;
                temp_counter_steiner[random_steiner]++;
                //If we havn't improve after from 5 steiner insertion or if we have increase the obtuses by 3, reset the simulated_cdt
                if(num_of_transition >= batch_size || delta_E >= (3*alpha) || delta_E == 0.000001){
                    simulate_cdt = best_cdt; //Reset to the best triangulation
                    curent_cdt = best_cdt;
                    num_of_transition = 0;
                    temp_p_sum = 0;
                    fill(temp_counter_steiner.begin(), temp_counter_steiner.end(), 0);
                    //Try to insert insert_steiner_around_centroid (3rd task)
                    if(i > max_iterations/1.5 && best_obtuse_faces > 1 && run_auto_method) {
                        try_steiner_around_centroid(simulate_cdt, polygon, temp_random_steiner);
                        temp_counter_steiner[5]++;
                        try_randomization = true;
                        num_of_transition++;
                        curent_cdt = simulate_cdt;
                        obtuse_faces = count_obtuse_triangles(simulate_cdt, polygon);
                        if(best_obtuse_faces > obtuse_faces) {
                            cout<<"Random steiner inserted: "<<temp_random_steiner<<endl;
                            vector_random_steiners.emplace_back(temp_random_steiner);
                            randomization = true;
                            count_steiners[5]++;
                            fill(temp_counter_steiner.begin(), temp_counter_steiner.end(), 0);
                            best_cdt = curent_cdt;
                            simulate_cdt = curent_cdt;
                            best_obtuse_faces = obtuse_faces;
                            num_of_transition = 0;
                            try_randomization = false;  
                        }
                    }
                }
                break;
            }
            //Go up to the "valley", Update temperature (increase)
            if (T < 1.0 && ((i > 180 && i < 190) || (i > 320 && i < 330))) T = T*1.4;
            if (T < 1.0 && ((i > 440 && i < 450) || (i > 560 && i < 570))) T = T*1.4;
            if (T < 1.0 && ((i > 680 && i < 690) || (i > 830 && i < 840))) T = T*1.4;
            if (T < 1.0 && ((i > 940 && i < 950))) T = T*1.4;
            
            //Case that we didn't insert this steiner into simulate_cdt. So, take back the previous simulate_cdt (custom_cdt)
            simulate_cdt = curent_cdt;
        }
        //Update temperature (decrease)
        T = T*(cooling_rate);
        //cout<<"Iteration: " <<i<< ", T: "<<T<<", best_obtuse_faces: "<<best_obtuse_faces<<" best_E: "<<best_E<<endl; 
    }

    //"Return" the best cdt
    custom_cdt = best_cdt;
    if(run_auto_method){
        double front;
        best_num_steiner = best_cdt.number_of_vertices() - init_vertices;
        if (best_num_steiner > 1)
            front = abs(1.0/(best_num_steiner - 1.0));
        else front = 0.0;    
        double rate_of_convergence = front * p_sum_best;
        
        std_string method_name = "SA";
        
        obtuse_faces = count_obtuse_triangles(best_cdt, polygon);
        best_E  = calculate_energy(obtuse_faces, best_num_steiner, alpha, beta);
        
        method_output(count_steiners, method_name, name_of_instance, best_num_steiner, init_num_obtuses, 
                        best_obtuse_faces, randomization, vector_random_steiners, rate_of_convergence, best_E,
                        subset, category);
    }
    time(&end_time);
    double time_taken = double(end_time - start_time); 
    cout<<"Time taken by  : "<<name_of_instance<<" is : "<<" sec "<<time_taken<<endl;    
}


//Ant colony method
void ant_colony(Custom_CDT& custom_cdt, Polygon& polygon, const double& alpha, const double& beta, 
                const double& chi, const double& psi, const double& lamda, const int& L, const int& kappa, 
                const std_string& name_of_instance, bool& randomization, vector<int> subset, std_string category,
                const bool& run_auto_method){
    int init_vertices = custom_cdt.number_of_vertices();
    int obtuse_faces = count_obtuse_triangles(custom_cdt, polygon);
    int init_num_obtuses = obtuse_faces;
    int new_obtuse_faces = obtuse_faces; 
    int best_obtuses = new_obtuse_faces;
    //int best_vertices = init_vertices;
    double best_E = calculate_energy(new_obtuse_faces, 0, alpha, beta);
    int counter_steiner = 0;
    int count_ants = kappa;
    bool obtuse_neighbors = false, can_we_use_adjacent = false, progress = false, try_randomization = false;
    //If we dont have progress for 4 loops, use random steiner
    int progress_counter = 0, progress_obtuses = obtuse_faces;
    double ro = 0.0;
    //3rd task
    Point_2 random_steiner;
    vector<Point_2> vector_random_steiners;
    vector<Point_2> inserted_steiners;
    double p_sum = 0.0;
    time_t start_time, end_time;
    //Let the program running for 15 loops without progress
    int non_progress_counter = 15;
    time(&start_time);
    int num_of_steiners = 0;
    vector<int> count_steiners(6, 0);
    //Midpoint edge: We need this edge to check if the steiner was entered on the boundary
    Segment_2 longest_edge;
    //Projection edge: We need this edge to check if the steiner was entered on the boundary
    Segment_2 opposite_edge;
    //Initialize the vector of ants
    vector<Ant> ants(count_ants);
    //The last ants after checking for conflicts
    vector<Ant> ant_last_winners_vector;
    //Vector to keep the ants have reduced obtuses
    vector<Ant> ant_reduce_obtuses_vector;

    //Ant best_cycle_ant;
    Ant::initialize_Ants(ants, custom_cdt);

    //Initialize pheromones    
    vector<double> taf(NUM_METHODS);
    vector<double> hta(NUM_METHODS);
    vector<double> delta_taf(NUM_METHODS);
    
    for (int i = 0; i < NUM_METHODS; ++i) {
        taf[i] = 0.5;
        hta[i] = 0.5;
        delta_taf[i] = 0.0;
    }

    Point_2 curent_steiner_point;
    Custom_CDT curent_cdt = custom_cdt;
    Custom_CDT best_cdt = custom_cdt;
    Custom_CDT random_cdt = custom_cdt;
    SteinerMethod curent_method;

    //3rd task, choose steiner from values
    std::mt19937 rng(std::random_device{}()); //Initialize RNG
    vector<int> values = subset; //Define possible values
    std::uniform_int_distribution<int> dist(0, values.size() - 1); //Generate index
    bool choose_auto_method = false;  
    //Start the L cycles
    for (int cycle = 0; cycle < L; ++cycle) {
        if (new_obtuse_faces == 0) break;
        //Clean the vectors
        ant_reduce_obtuses_vector.clear();
        ant_last_winners_vector.clear();      
        //Ants
        for (int ant_index = 0; ant_index < count_ants; ++ant_index) {
            ants[ant_index].set_Custom_CDT(curent_cdt);
            //Chose obtuse face and check it, give_random_obtuse has check is_obtuse(face), is_face_inside_region(face, polygon)
            Face_handle face = give_random_obtuse(ants[ant_index].get_Custom_CDT(), polygon);
            /*Improve triangulation*/
            ro = calculate_radius_to_height(face, curent_cdt);
            obtuse_neighbors = has_obtuse_neighbors(curent_cdt, face, polygon);

            //If value vector is [0,1,2,3,4] use the selectSteinerMethod, else, choose method from values
            if(choose_auto_method) curent_method = (SteinerMethod)values[dist(rng)];
            else curent_method = selectSteinerMethod(ro, taf, hta, chi, psi, obtuse_neighbors);
            
            switch(curent_method){
                //If circumcenter steiner is outside of the boundary or the opposite edge of obtuse vertex is constraint, use the centroid
                case 0: 
                    if(!insert_circumcenter(curent_cdt, face, polygon, curent_steiner_point)){
                        insert_centroid(curent_cdt, face, polygon, curent_steiner_point); 
                        curent_method = CENTROID;
                        break;
                    }
                    else break;
                case 1: insert_midpoint(curent_cdt, face, polygon, curent_steiner_point, longest_edge); break;
                case 2: insert_projection(curent_cdt, face, polygon, curent_steiner_point, opposite_edge); break;
                case 3:
                    //If the face has obtuse neighbor(s) and the polygon of adjacent points is convex, then insert the adjacent steiner
                    can_we_use_adjacent = insert_adjacent_steiner(curent_cdt, face, polygon, curent_steiner_point);
                    if((!can_we_use_adjacent)){
                        insert_projection(curent_cdt, face, polygon, curent_steiner_point, opposite_edge);
                        curent_method = PROJECTION;
                    }
                    break;
                default: break;
            }
            //Save the No of method into Ant
            ants[ant_index].set_steiner_method(curent_method);
            //Save the curent cdt into Ant
            ants[ant_index].set_Custom_CDT(curent_cdt);
            //Save the steiner into Ant   
            ants[ant_index].set_steiner(curent_steiner_point);
            
            ants[ant_index].set_num_of_obtuses(count_obtuse_triangles(ants[ant_index].get_Custom_CDT(), polygon));
            new_obtuse_faces = ants[ant_index].get_num_of_obtuses();
            counter_steiner = ants[ant_index].get_Custom_CDT().number_of_vertices() - init_vertices;
            //Save the energy into Ant
            ants[ant_index].set_energy( calculate_energy(new_obtuse_faces, counter_steiner, alpha, beta) );
            //Save the DeltaE into Ant
            ants[ant_index].set_DeltaE( (ants[ant_index].get_energy() - best_E) );
            /*Evaluate the resulting triangulation*/
            //Update if this steiner improve the triangulation
            if (ants[ant_index].get_DeltaE() < 0) {
                ants[ant_index].set_reduce_obtuses(true);
                //Save and the longest_edge, opposite_edge because we need it to update the polygon, if the steiner is on boundary
                if((ants[ant_index].get_steiner_method() == 1) && polygon.bounded_side(ants[ant_index].get_steiner_point()) == CGAL::ON_BOUNDARY) 
                    ants[ant_index].set_longest_edge_midpoint(longest_edge);
                if((ants[ant_index].get_steiner_method() == 2) && polygon.bounded_side(ants[ant_index].get_steiner_point()) == CGAL::ON_BOUNDARY) 
                    ants[ant_index].set_opposite_edge_projection(opposite_edge);
            }
            else ants[ant_index].set_reduce_obtuses(false);
            
            if(run_auto_method){
                //If the try_randomization is activated, use other cdt
                if(!try_randomization) curent_cdt = best_cdt;
                else curent_cdt = random_cdt;
            }
            else curent_cdt = best_cdt;
        }
        
        //Save the bests ants (not the last winners)
        for (int ant_index = 0; ant_index < count_ants; ++ant_index){
            //If this ant didnt reduce the obtuses faces of cdt, ignore it
            if(!ants[ant_index].get_reduce_obtuses()) continue;
            //Compare the best_cdt and ant_cdt to find the differnces in the faces.
            /*If there are faces in best_cdt that ant_cdt does not have, this means that those faces were affected by ant*/
            affected_faces(best_cdt, ants[ant_index]);
            //Add this ant into ant_reduce_obtuses_vector
            ant_reduce_obtuses_vector.emplace_back(ants[ant_index]);
        }
        
        //If we had only 1 ant that improve the triangulation
        if(ant_reduce_obtuses_vector.size() == 1) {
            ant_last_winners_vector.emplace_back(ant_reduce_obtuses_vector[0]);
        }
        //Or >=2
        if(ant_reduce_obtuses_vector.size() >= 2) {
            ant_last_winners_vector = save_the_best(ant_reduce_obtuses_vector);
        }

        //3rd task, 
        int num_of_obtuses_before = 0, num_of_obtuses_after = 0;
        /*Save the best triangulation*/
        for(int i = 0; i < ant_last_winners_vector.size(); i++){
           
            //We add in the place 4 of vector the centroid
            if(ant_last_winners_vector[i].get_steiner_method() == CENTROID){
                ant_last_winners_vector[i].set_steiner_method(NUM_METHODS);
            }
            num_of_obtuses_before = count_obtuse_triangles(best_cdt, polygon);
            //3rd task
            if(try_randomization && run_auto_method) {
                best_cdt.insert_no_flip(random_steiner);
                start_the_flips(best_cdt, polygon);
                cout<<fixed<<"Random steiner inserted: "<<random_steiner<<endl;
                inserted_steiners.emplace_back(random_steiner);
                
            }
            inserted_steiners.emplace_back(ant_last_winners_vector[i].get_steiner_point());
            
            count_steiners[ant_last_winners_vector[i].get_steiner_method()]++;
            best_cdt.insert_no_flip(ant_last_winners_vector[i].get_steiner_point());
            start_the_flips(best_cdt, polygon);
            if(run_auto_method){
                //3rd task, p_sum
                num_of_steiners = best_cdt.number_of_vertices() - init_vertices;
                num_of_obtuses_after = count_obtuse_triangles(best_cdt, polygon);
                p_sum += p_sum_function(num_of_steiners - 1, num_of_obtuses_before, num_of_obtuses_after);
            }
           
            curent_steiner_point = ant_last_winners_vector[i].get_steiner_point();
            longest_edge = ant_last_winners_vector[i].get_longest_edge_midpoint();
            opposite_edge = ant_last_winners_vector[i].get_opposite_edge_projection();
            curent_method = ant_last_winners_vector[i].get_steiner_method();

            if((curent_method == 1) && (polygon.bounded_side(curent_steiner_point) == CGAL::ON_BOUNDARY))
                update_polygon(polygon, curent_steiner_point, longest_edge.source(), longest_edge.target());
            if((curent_method == 2) && (polygon.bounded_side(curent_steiner_point) == CGAL::ON_BOUNDARY))
                update_polygon(polygon, curent_steiner_point, opposite_edge.source(), opposite_edge.target());
        }
        
        //Update the best triangulation and the best_E
        best_obtuses = count_obtuse_triangles(best_cdt, polygon);
        if(run_auto_method){
            //3rd task. Check for 15 loops if we dont have progress, to activate the random steiner
            if(best_obtuses < progress_obtuses){
                progress = true;
                
                progress_obtuses = best_obtuses;
                progress_counter = 0;
                if(try_randomization) {
                    randomization = true;
                    vector_random_steiners.emplace_back(random_steiner);
                }
                try_randomization = false;
                
            }
            else  progress_counter++;

            //Reset the curent cdt and random_cdt
            curent_cdt = best_cdt;
            random_cdt = best_cdt;
            //Try to insert insert_steiner_around_centroid (3rd task)
            if (progress_counter >= non_progress_counter) try_randomization = true;
            if (try_randomization) {
                non_progress_counter = 10;
                try_randomization = true;
                progress_counter = 0;
                try_steiner_around_centroid(random_cdt, polygon, random_steiner);
                int obtuses = count_obtuse_triangles(random_cdt, polygon);
                curent_cdt = random_cdt;
                if(obtuses < best_obtuses) {
                    best_cdt = random_cdt;
                    curent_cdt = best_cdt;
                    best_obtuses = obtuses;
                    progress_obtuses = obtuses;
                    counter_steiner = best_cdt.number_of_vertices() - init_vertices;
                    randomization = true;
                    vector_random_steiners.emplace_back(random_steiner);
                    cout<<fixed<<"Random steiner inserted: "<<random_steiner<<endl;
                }
            }
        }
        
        counter_steiner = best_cdt.number_of_vertices() - init_vertices;
        best_E = calculate_energy(best_obtuses, counter_steiner, alpha, beta);
        
        /*Update pheromones*/
        if(ant_reduce_obtuses_vector.size() > 0) updatePheromones(taf, delta_taf, ant_reduce_obtuses_vector, lamda);
        ///Restart the ants
        Ant::initialize_Ants(ants, best_cdt);
    }    
    cout<<endl;
    //Return the best cdt
    custom_cdt = best_cdt;
    if(run_auto_method){
        //3rd task
        double front;
        if (counter_steiner > 1)
            front = abs(1.0/(counter_steiner - 1.0));
        else front = 0.0;

        double rate_of_convergence = front * p_sum;
        std_string method_name = "Ant";
        new_obtuse_faces = count_obtuse_triangles(best_cdt, polygon);

        //Update the index 5 "how many times random steiner choosed"
        count_steiners[5] = vector_random_steiners.size();
        best_E = calculate_energy(new_obtuse_faces, counter_steiner, alpha, beta);
        method_output(count_steiners, method_name, name_of_instance, counter_steiner, init_num_obtuses, new_obtuse_faces, 
                    randomization, vector_random_steiners, rate_of_convergence, best_E, subset, category);
    }
    time(&end_time);
    double time_taken = double(end_time - start_time); 
    cout<<"Time taken by program: "<<name_of_instance<<" is : "<<time_taken<<" sec "<<endl;
}

//Check for conflict between 2 ants
bool have_conflict(Ant& ant1, Ant& ant2){
  
    //Take the affected faces of the
    set<Face_handle> face1 = ant1.get_affected_faces();
    set<Face_handle> face2 = ant2.get_affected_faces();
    //For every face of every of 2 ants check if they have same face(s)
    for (const auto& temp_face1 : face1){
        for (const auto& temp_face2 : face2){
            if(are_faces_equal(temp_face1, temp_face2)){
                return true;
            }
        }
    }
    return false;
}

//Keep the final ants that we take for the best triangulation
vector<Ant> save_the_best(vector<Ant>& ants){
    vector<Ant> winners;
    //Check each ant with the next ones
    for (int i = 0; i < (ants.size()-1); ++i){
        for (int j = i+1; j < ants.size(); ++j){
            
            if(ants[i].get_reduce_obtuses() && ants[j].get_reduce_obtuses()){
                //Check these ants if they have conflict
                if(have_conflict(ants[i], ants[j])){
                    //If an ant is loser, skip this conflict
                    if(ants[i].get_conflict_loser() || ants[j].get_conflict_loser()) continue;
                    //Update that these ants had a conflict
                    ants[i].set_conflict(true); 
                    ants[j].set_conflict(true);
                    //Choose the winner (less energy)
                    if(ants[i].get_energy() <= ants[j].get_energy()) {
                        ants[j].set_conflict_loser(true); //j ant is the loser, i the winner
                        continue; //continue in case another j with less energy is found
                    }
                    else{
                        ants[i].set_conflict_loser(true); //i ant is the loser, j the winner
                        break; //there is no point in looking any further because the i lost
                    }
                }
            }
        }
    }
    
    //The other ants without conflict and with negative Delta_e, we have to insert them into winners
    for (int i = 0; i < ants.size(); ++i){
        //Insert the ants without conflicts
        if(ants[i].get_reduce_obtuses() && !ants[i].get_conflict()){
            winners.emplace_back(ants[i]);
            continue;
        }
        //Insert the winners after conflict (not losers = winners)
        if(!ants[i].get_conflict_loser()){
            winners.emplace_back(ants[i]);
        } 
    }
    return winners;
}

//Check if two faces are equals
bool are_faces_equal(const Face_handle& face1, const Face_handle& face2) {
    set<Point_2> vertices_face1 = {
        face1->vertex(0)->point(),
        face1->vertex(1)->point(),
        face1->vertex(2)->point()
    };

    set<Point_2> vertices_face2 = {
        face2->vertex(0)->point(),
        face2->vertex(1)->point(),
        face2->vertex(2)->point()
    };
    //Copy the sets into vectors for sorting and comparison
    vector<Point_2> points_vector_face1(vertices_face1.begin(), vertices_face1.end());
    vector<Point_2> points_vector_face2(vertices_face2.begin(), vertices_face2.end());

    //Lexicographical Order of the face points. First compare x-coordinates.
    sort(points_vector_face1.begin(), points_vector_face1.end());
    sort(points_vector_face2.begin(), points_vector_face2.end());

    return points_vector_face1 == points_vector_face2;
}

//Give a random obtuse face
Face_handle give_random_obtuse(Custom_CDT& custom_cdt, Polygon& polygon) {
    //Container to store faces with obtuse angles
    vector<Face_handle> obtuse_faces;
    obtuse_faces.clear();
    //Initialize the random number generator with a random device and engine
    static std::mt19937 generator(std::random_device{}()); //Only initialize once

    for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
        if (!is_face_inside_region(face, polygon)) continue;
        if (is_obtuse(face)) {
            obtuse_faces.push_back(face);
        }
    }

    //Check if there are no obtuse faces found
    if (obtuse_faces.empty()) {
        cerr<<"No obtuse faces found!"<<endl;
        return Face_handle(); //Return an invalid handle if no obtuse face is found
    }
    //Generate a random index in [0, obtuse_faces.size() - 1]
    uniform_int_distribution<> distribution(0, obtuse_faces.size() - 1);
    //Random index in [0, obtuse_faces.size() - 1]
    int index = distribution(generator);
    
    return obtuse_faces[index];
}

//Update affected faces
void affected_faces(Custom_CDT& best_cdt, Ant& ant) {
    ant.clear_ant_affect_faces();
    Custom_CDT temp_ant_cdt = ant.get_Custom_CDT();
    for (auto best_face = best_cdt.finite_faces_begin(); best_face != best_cdt.finite_faces_end(); ++best_face) {
        bool found_face = false;
        for (auto ant_face = temp_ant_cdt.finite_faces_begin(); ant_face != temp_ant_cdt.finite_faces_end(); ++ant_face) {
            //Check the faces if are equals
            if (are_faces_equal(best_face, ant_face)) {
                found_face = true;
                break;
            }
        }
        //If best_face != ant_face means that we have conflict
        if (!found_face) ant.set_face_in_ant_affect_faces(best_face);
    }
}

void printAntDetails(vector<Ant>& ants) {
    cout<<"Number of ants: "<<ants.size() <<endl;
    for (size_t i = 0; i < ants.size(); ++i) {
        cout<<"Ant " << i << " details:"<<endl;
        cout<<"Steiner method: "<<static_cast<int>(ants[i].get_steiner_method())<<endl;
        cout<<"Energy: "<<ants[i].get_energy()<<endl;
        cout<<"DeltE: "<<ants[i].get_DeltaE()<<endl;
        cout<<"Steiner point: ("<< ants[i].get_steiner_point().x()<<", "<< ants[i].get_steiner_point().y() << ")"<<endl;
        cout<<"Num of obtuses: "<<ants[i].get_num_of_obtuses()<<endl;
        cout<<"Ant reduces the obtuses? : "<<ants[i].get_reduce_obtuses()<<endl;
        cout<<"Conflict? : "<<ants[i].get_conflict()<<endl;
        
        if(ants[i].get_reduce_obtuses()){
            set<Face_handle> temp_face = ants[i].get_affected_faces();
            cout<<"size of vector: "<<ants[i].get_affected_faces().size()<<endl;
            for (const auto& face : temp_face) {
                for (int j = 0; j < 3; ++j) {  //A triangle has 3 vertices
                    cout<<endl;
                    auto vertex = face->vertex(j);
                    if (vertex != nullptr) cout<<"Point "<<j<<": ("<<vertex->point().x()<<", "<<vertex->point().y()<< ") =>";
                    else cout<<"Point "<<j<< ": (nullptr)"<<endl;            
                }
                cout<<endl;
            }
        }
        cout<<"*******************************"<<endl;
        
    }
}

SteinerMethod selectSteinerMethod(const double& ro, const vector<double>& taf, vector<double>& hta, double chi, double psi, bool obtuse_neighbors) {
    //Ensure inputs are valid
    if (taf.size() != SteinerMethod::NUM_METHODS || hta.size() != SteinerMethod::NUM_METHODS) {
        cerr<<"Error: taf or hta size does not match NUM_METHODS."<<endl;
        return static_cast<SteinerMethod>(0); 
    }

    //CIRCUMCENTER, MIDPOINT_LONGEST_EDGE, PROJECTION, ADJACENT_FACES
    hta[0] = hta_circumcenter(ro); 
    hta[1] = hta_midpoint(ro);    
    hta[2] = hta_vertex_projection(ro);
    hta[3] = hta_mean_adjacent(obtuse_neighbors);
   
    //Calculate probabilities for each Steiner point method
    vector<double> probabilities(SteinerMethod::NUM_METHODS, 0.0);
    double denominator = 0.0;
    //Calculate the (τ_sp)^χ and (η_sp)^ψ and Σ_i...
    for (int i = 0; i < SteinerMethod::NUM_METHODS; i++) {
        probabilities[i] = pow(taf[i], chi) * pow(hta[i], psi);
        denominator += probabilities[i];
    }

    //Check for zero denominator
    if (denominator == 0.0) {
        cerr<<"Error: Denominator for probability normalization is zero."<<endl;
        return static_cast<SteinerMethod>(0); //Or an alternative fallback
    }
    
    //Calculate the Σ...
    for (int i = 0; i < SteinerMethod::NUM_METHODS; i++) {
        probabilities[i] /= denominator;
    }
    //Sum up the elements in the probabilities vector, takes a range of elements (beginning and end iterators) and an initial value, 
    //and it returns the sum of the elements in that range plus the initial value
    double sum_probabilities = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    double correction = 1.0 - sum_probabilities;
    probabilities[SteinerMethod::NUM_METHODS - 1] += correction;

    //Select a method based on the computed probabilities
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double random_value = dis(gen);

    double cumulative_probability = 0.0;
    //Roulette wheel selection
    for (int i = 0; i < SteinerMethod::NUM_METHODS; ++i) {
        cumulative_probability += probabilities[i];
        if (random_value <= cumulative_probability) {
            return static_cast<SteinerMethod>(i);
        }
    }

    //Should not reach here if probabilities are normalized correctly
    cerr<<"Warning: Random value out of cumulative probability range."<<endl;
    return static_cast<SteinerMethod>(SteinerMethod::NUM_METHODS - 1);
}

void updatePheromones(vector<double>& taf, vector<double>& delta_taf, vector<Ant> selected_ants, double lamda) {
    SteinerMethod sp;
    //If a method has been selected at least once, set a value of 1 in the same index of steinerMethod
    vector<int> num_of_methods(taf.size(), 0);
    //Reinforce pheromones (delta_taf) based on ant improvements
    for (int i = 0; i < selected_ants.size(); ++i) {
        sp = selected_ants[i].get_steiner_method();
        if (sp < 0) {
            cerr<<"Error: Invalid SteinerMethod index "<<sp<<endl;
            continue;
        }
        //We dont have pheromones for centroid
        if(sp == CENTROID) continue;
        
        if(num_of_methods[sp] == 1) {
            delta_taf[sp] = delta_taf[sp] + (1.0 / (1.0 + selected_ants[i].get_energy()));
            continue;
        }
        //Only update if improvement is true (i.e., reduces obtuse triangles)
        if (selected_ants[i].get_reduce_obtuses()) delta_taf[sp] = 1.0 / (1.0 + selected_ants[i].get_energy());
        num_of_methods[sp] = 1;
    }

    //For each method not selected by ants, set the ΔΕ zero
    for (int i = 0; i < num_of_methods.size(); ++i) if(num_of_methods[i] == 0) delta_taf[i] = 0;

    //Update the τ
    for (int i = 0; i < delta_taf.size(); ++i) taf[i] = (1.0 - lamda) * taf[i] + delta_taf[i];
}

//Function to compute ρ (radius-to-height ratio)
double calculate_radius_to_height(const Face_handle& face, const Custom_CDT& cdt) {
    //Vertices of the triangle
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();

    //Compute circumradius (R)
    Point circumcenter = CGAL::circumcenter(p1, p2, p3);
    double R = sqrt(CGAL::to_double(CGAL::squared_distance(p1, circumcenter)));

    //Compute lengths of edges
    double d1 = sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));
    double d2 = sqrt(CGAL::to_double(CGAL::squared_distance(p2, p3)));
    double d3 = sqrt(CGAL::to_double(CGAL::squared_distance(p3, p1)));

    //Find the longest edge (base)
    double longest_edge = max({d1, d2, d3});

    //Compute the height from the opposite vertex
    double s = (d1 + d2 + d3) / 2.0; //Semi-perimeter
    double area = sqrt(s * (s - d1) * (s - d2) * (s - d3)); //Triangle area (Heron's Formula)
    double height = (2 * area) / longest_edge;

    //Calculate ρ
    return R / height;
}

//True if the face has at least 1 obtuse neighbor
bool has_obtuse_neighbors(const Custom_CDT& custom_cdt, const Face_handle& face, const Polygon& polygon){
    for (int i = 0; i < 3; ++i) {
        Face_handle neighbor = face->neighbor(i);
        //Take the opposite edge from the vertex i 
        if (custom_cdt.is_constrained(make_pair(face, i))) continue;
        //Because we have just an iterator, this iterator may go out of bounds
        //because he doesn't have the supervision of for loop
        if(custom_cdt.is_infinite(neighbor) || !is_face_inside_region(neighbor, polygon)) continue;
        if(is_obtuse(neighbor)) return true;
    }
    return false;
}

//Function to compute η_vertex_projection - Vertex Projection: strategy most effective when ρ > 2.0
double hta_vertex_projection(double rho) {
    //Avoid the next calculate with this if case
    if (rho <= 1.0) return 0.0;
    return max(0.0, ((rho - 1.0) / rho));
}

//Function to compute η_circumcenter- Circumcenter: Probable when ρ ∈ [1.0, 2.0]
double hta_circumcenter(double rho) {
    return rho / (2.0 + rho);
}

//Function to compute η_midpoint
double hta_midpoint(double rho) {
    //Avoid the next calculate with this if case
    if (rho >= 1.5) return 0.0;
    return max(0.0, ((3.0 - (2.0 * rho))) / 3.0);
}

// Function to compute η_mean_adjacent (returns 1 if has_obtuse_neighbors is true)
double hta_mean_adjacent(bool has_obtuse_neighbors) {
    return (has_obtuse_neighbors) ? 1.0 : 0.0;
}

//Ιnsert Steiner points at circumcenter
bool insert_circumcenter(Custom_CDT& circumcenter_cdt, const Face_handle& face, const Polygon& polygon, Point_2& circumcenter_steiner) {    
    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point();
    //Find the obtuse vertex
    Point_2 obtuse_vertex = find_obtuse_vertex(p1, p2, p3);

    //Identify the opposite edge of the obtuse vertex
    int opposite_edge_index;
    if (obtuse_vertex == p1) opposite_edge_index = 0;       //Edge (p2, p3)
    else if (obtuse_vertex == p2) opposite_edge_index = 1;  //Edge (p1, p3)
    else if (obtuse_vertex == p3) opposite_edge_index = 2;  //Edge (p1, p2)

    auto edge = make_pair(face, opposite_edge_index);
    //Check if the opposite edge is constrained
    if (circumcenter_cdt.is_constrained(edge)) {
        return false; //Do not insert if the edge is constrained
    }

    //Compute the circumcenter of the triangle
    Point_2 circumcenter = CGAL::circumcenter(p1, p2, p3);
    //Circumcenter must be inside of the boundary
    if (is_point_inside_region(circumcenter, polygon) && is_circumcenter_in_neighbor(circumcenter_cdt, face, circumcenter)){
        if(is_convex(p1, p2, p3, circumcenter)){
            //Check the face that the circumcenter will enter must be inside of the boundary
            Face_handle locate_face = circumcenter_cdt.locate(circumcenter);
            if(is_face_inside_region(locate_face, polygon)){
                circumcenter_steiner = circumcenter;
                circumcenter_cdt.insert_no_flip(circumcenter);
                start_the_flips(circumcenter_cdt, polygon);
                return true;
            }
        }
    }
    return false;
}

//Check if the circumcenter is on side of a neighbor or not
bool is_circumcenter_in_neighbor(const Custom_CDT& cdt, const Face_handle& face, const Point_2& circumcenter) {
    for (int i = 0; i < 3; ++i) {
        Face_handle neighbor_face = face->neighbor(i);

        //Skip invalid or infinite neighbors
        if (cdt.is_infinite(neighbor_face)) continue;

        //Get vertices of the neighbor face
        Point_2 v1 = neighbor_face->vertex(0)->point();
        Point_2 v2 = neighbor_face->vertex(1)->point();
        Point_2 v3 = neighbor_face->vertex(2)->point();

        //Create a Triangle_2 object
        CGAL::Triangle_2<Custom_CDT::Geom_traits> triangle(v1, v2, v3);

        //Check if the circumcenter is inside or on the boundary of the triangle
        CGAL::Bounded_side side = triangle.bounded_side(circumcenter);
        if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY) {
            return true; //Circumcenter is inside or on the boundary of the neighbor face
        }
    }
    return false; //Circumcenter is not in or on any neighbor face
}

//Ιnsert Steiner points at centroid
void insert_centroid(Custom_CDT& centroid_cdt, const Face_handle& face, const Polygon& polygon, Point_2& centroid_steiner) {       
    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point();
    //Compute the centroid of the triangle
    Point_2 centroid = CGAL::centroid(p1, p2, p3);
    centroid_steiner = centroid;
    centroid_cdt.insert_no_flip(centroid);
    start_the_flips(centroid_cdt, polygon);
}

//Ιf we added steiner on boundary of the polygon, update the new edges of the polygon
void update_polygon(Polygon& polygon, const Point_2& steiner_point, const Point_2& p1, const Point_2& p2) {
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

//Return the number of vertices in a cdt
int count_vertices(const Custom_CDT& cdt) {
    int count = 0;
    for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
        ++count;
    }
    return count;
}

//Flips method
void start_the_flips(Custom_CDT& cdt, const Polygon& polygon){
    bool progress = true;
    while(progress){
        progress = false;
        for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
            //The face that contains this edge
            Face_handle f1 = edge->first;
            //The index of the edge within the face edge->first
            int i = edge->second;
            Face_handle f2 = f1->neighbor(i);
            
            if(cdt.is_infinite(f1) || cdt.is_infinite(f2)) continue;
            if(cdt.is_constrained(*edge)) continue;
            if(!(is_face_inside_region(f1, polygon)) || !(is_face_inside_region(f2, polygon))) continue;

            Point_2 p1 = f1->vertex(cdt.ccw(i))->point(); //First vertex on the shared edge (Counter-Clock Wise)
            Point_2 p3 = f1->vertex(cdt.cw(i))->point();  //Second vertex on the shared edge (Clock Wise)
            Point_2 p2 = f1->vertex(i)->point();          //Opposite vertex in the first triangle

            //Check if edge is on boundary or inside of the boundary
            if(!is_edge_inside_region(p1, p3, polygon)) continue;
            if(is_edge_on_boundary(p1, p3, polygon)) continue;

            //Mirror index gets the index of the vertex in f2 that is opposite to this shared edge
            int mirror_index = cdt.mirror_index(f1, i);
            Point_2 p4 = f2->vertex(mirror_index)->point();         
            if(is_it_worth_flip(p1, p2, p3, p4)){
                cdt.flip(f1, i);
                progress = true;
                break;
            }
            else continue;
        }
    }
}

//If 1 point is on the boundary
bool is_point_inside_region(const Point_2& point, const Polygon& polygon) {
    //Check if the point is inside the polygon
    return (polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE) || (polygon.bounded_side(point) == CGAL::ON_BOUNDARY);
}

//If face is inside of region boundary
bool is_face_inside_region(const Face_handle& face, const Polygon& polygon) {
    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point();
    //Calculate the centroid of the triangle
    Point_2 centroid = CGAL::centroid(p1, p2, p3);

    //Check if the centroid is inside or on the boundary of the polygon
    bool centroid_inside = (polygon.bounded_side(centroid) == CGAL::ON_BOUNDED_SIDE ||
                            polygon.bounded_side(centroid) == CGAL::ON_BOUNDARY);

    if (!centroid_inside) return false;
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
    return true;
}

//If edge is inside of region boundary
bool is_edge_inside_region(const Point_2& p1, const Point_2& p2, const Polygon& polygon){
    bool mids_inside_region = 
        (polygon.bounded_side(CGAL::midpoint(p1, p2)) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(CGAL::midpoint(p1, p2)) == CGAL::ON_BOUNDARY);
    bool points_inside_region = 
        (polygon.bounded_side(p1) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(p1) == CGAL::ON_BOUNDARY) && 
        (polygon.bounded_side(p2) == CGAL::ON_BOUNDED_SIDE || polygon.bounded_side(p2) == CGAL::ON_BOUNDARY);
    
    //Check if this edge indersect with an edge of polygon (boundary)
    /*Segment_2 edge(p1, p2);
    for (auto edge_it = polygon.edges_begin(); edge_it != polygon.edges_end(); ++edge_it) {
        Segment_2 poly_edge(*edge_it);
        auto result = CGAL::intersection(edge, poly_edge);
        if (result) {
            if (const Point_2* point = boost::get<Point_2>(&*result)) {
                //If intersection is a point, check if it's not an endpoint
                if (*point != p1 && *point != p2 && edge.has_on(*point)) {
                    return false; //Edge intersects at a point other than endpoints
                }
            }
        }
    }*/
    return mids_inside_region && points_inside_region;
}

//Function to check if an edge is part of the boundary of the polygon
bool is_edge_on_boundary(const Point_2& p1, const Point_2& p2, const Polygon& polygon) {
    for (auto edge_it = polygon.edges_begin(); edge_it != polygon.edges_end(); ++edge_it) {
        if ((edge_it->source() == p1 && edge_it->target() == p2) ||
            (edge_it->source() == p2 && edge_it->target() == p1)) {
            return true;
        }
    }
    return false;
}

Point_2 find_obtuse_vertex(const Point_2& v1, const Point_2& v2, const Point_2& v3) {
    //Calculate squared distances
    double ab2 = CGAL::to_double(squared_distance(v1, v2));
    double ac2 = CGAL::to_double(squared_distance(v1, v3));
    double bc2 = CGAL::to_double(squared_distance(v2, v3));

    if (ab2 + ac2 < bc2) return v1; //obtuse angle at v1
    if (ab2 + bc2 < ac2) return v2; //obtuse angle at v2
    if (ac2 + bc2 < ab2) return v3; //obtuse angle at v3

    throw logic_error("No obtuse angle found in the triangle.");
}

//Find the longest edge of a face
Segment_2 find_longest_edge(const Face_handle& face) {

    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point();
    //Create segments for each edge of the triangle
    Segment_2 edge1(p1, p2);
    Segment_2 edge2(p2, p3);
    Segment_2 edge3(p3, p1);

    //Compare the squared lengths of the segments to find the longest one
    auto length1 = CGAL::squared_distance(p1, p2);
    auto length2 = CGAL::squared_distance(p2, p3);
    auto length3 = CGAL::squared_distance(p3, p1);

    //Determine which edge is the longest
    if (length1 >= length2 && length1 >= length3) {
        return edge1;
    } else if (length2 >= length1 && length2 >= length3) {
        return edge2;
    } else {
        return edge3;
    }
}

void print_polygon_edges(const Polygon& polygon){
    cout << "Polygon edges:\n";
    for (auto it = polygon.edges_begin(); it != polygon.edges_end(); ++it) {
        Point_2 p1 = it->source();
        Point_2 p2 = it->target();
        cout << "(" << p1.x() << ", " << p1.y() << ") -> "<< "(" << p2.x() << ", " << p2.y() << ")\n";
    }
}

bool is_polygon_convex(const set<Point_2>& unique_points) {
    //Convert the set to a vector
    vector<Point_2> points(unique_points.begin(), unique_points.end());
    //Ensure we have at least 3 points for a polygon
    if (points.size() < 3) {
        return false; //A polygon cannot be formed
    }

    //Compute the convex hull of the points
    vector<Point_2> convex_hull_points;
    CGAL::convex_hull_2(points.begin(), points.end(), back_inserter(convex_hull_points));

    //A polygon is convex if all its points are part of its convex hull
    if (convex_hull_points.size() != points.size()) {
        return false;
    }    
    return true;
}

//Function to check if a vertex in custom_cdt is a Steiner point
bool is_steiner_point(Vertex_handle vertex, const vector<Point_2>& original_points) {
    const Point_2& vertex_point = vertex->point();
    //Check if vertex_point is in the original_points vector
    return find(original_points.begin(), original_points.end(), vertex_point) == original_points.end();
}


void output(value jv, Custom_CDT custom_cdt, vector<Point_2> original_points, int obtuse_count, std_string output_path, bool randomization){
    //Get instance_uid from input JSON
    std_string instance_uid = jv.as_object().at("instance_uid").as_string().c_str();

    //Extract method & parameters
    std_string method = jv.as_object().at("method").as_string().c_str();
    boost::json::object parameters = jv.as_object().at("parameters").as_object();

    //Prepare steiner points lists
    vector<std_string> steiner_points_x;
    vector<std_string> steiner_points_y;

    for (auto vertex = custom_cdt.finite_vertices_begin(); vertex != custom_cdt.finite_vertices_end(); ++vertex){
        // if a vertex is not found in the initial points then it must be a steiner point
        if (is_steiner_point(vertex, original_points)){
            const Point_2& p = vertex->point();
            steiner_points_x.push_back(convert_to_string(p.x())); //Store x-coordinate as string
            steiner_points_y.push_back(convert_to_string(p.y())); //Store y-coordinate as string
        }
    }
    
    //Convert steiner_points_x and steiner_points_y to JSON-compatible format
    boost::json::array steiner_points_x_json;
    boost::json::array steiner_points_y_json;

    for (const auto &x : steiner_points_x){
        steiner_points_x_json.push_back(boost::json::value(x));
    }

    for (const auto &y : steiner_points_y){
        steiner_points_y_json.push_back(boost::json::value(y));
    }

    //Prepare edges list
    vector<pair<int, int>> edges;
    //vector<Vertex_handle, int> vertex_index_map;
    std::map<Vertex_handle, int> vertex_index_map;
    int index = 0;

    //Map vertices to unique indices
    for (auto vertex = custom_cdt.finite_vertices_begin(); vertex != custom_cdt.finite_vertices_end(); ++vertex){
        vertex_index_map[vertex] = index++;
    }

    //Collect edges
    for (auto edge = custom_cdt.finite_edges_begin(); edge != custom_cdt.finite_edges_end(); ++edge){
        auto v1 = edge->first->vertex((edge->second + 1) % 3);
        auto v2 = edge->first->vertex((edge->second + 2) % 3);
        edges.emplace_back(vertex_index_map[v1], vertex_index_map[v2]);
    }

    //Convert to JSON arrays
    boost::json::array edges_json, parameters_json;

    for (const auto &edge : edges)
        edges_json.push_back(boost::json::array{edge.first, edge.second});

    //Convert parameters to JSON
    for (const auto &[key, value] : parameters){
        if (value.is_double()) {
            //Convert double values to a formatted string
            parameters_json.push_back(boost::json::object{{key, format_double(value.as_double())}});
        } 
        else {
            //Handle non-double values as usual
            parameters_json.push_back(boost::json::object{{key, value}});
        }
    }
  
    //make the json file look pretty (not in a single line)
    ostringstream json_output;
    json_output<<"{\n";
    json_output<<"  \"content_type\": \"CG_SHOP_2025_Solution\",\n";
    json_output<<"  \"instance_uid\": \""<< instance_uid<<"\",\n";
    json_output<<"  \"steiner_points_x\": "<<boost::json::serialize(steiner_points_x_json)<<",\n";
    json_output<<"  \"steiner_points_y\": "<<boost::json::serialize(steiner_points_y_json)<<",\n";
    json_output<<"  \"edges\": "<<boost::json::serialize(edges_json)<<",\n";
    json_output<<"  \"obtuse_count\": \"" <<obtuse_count<<"\",\n";
    json_output<<"  \"method\": \"" <<method<<"\",\n";
    json_output<<"  \"parameters\": "<<parameters_json<<",\n";
    json_output<<"  \"randomization\":  "<<(randomization ? "true" : "false")<<"\n";
    json_output<<"}\n";
    
    //Write to file
    ofstream output_file(output_path);
    output_file<<json_output.str();
    output_file.close();

    cout<<"Solution JSON file written as solution_output.json"<<endl;
}

//Function to format double values as strings
std_string format_double(double value) {
    ostringstream oss;
    oss.precision(1); //Adjust precision if necessary
    oss<<fixed<<value; //Use fixed-point notation
    return oss.str();
}

//Convert coord into string like franction. If it is int, return string like int
std_string convert_to_string(const K::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    ostringstream oss;
    //Check if the coordinate is an integer
    if (exact_coord.get_den() == 1.0) {
        oss<<fixed<<setprecision(0)<<coord;
        return oss.str();   
    }
    oss<<fixed<<setprecision(0)<<exact_coord.get_num()<<"/"<<exact_coord.get_den();
    
    //Return the formatted string
    return oss.str();
}

/*3rd Task*/
bool are_constraints_closed(const vector<pair<int, int>>& additional_constraints, int num_points, 
                            const vector<Point_2>& points, const Polygon& polygon) {
    if(additional_constraints.empty()) return false;
    
    vector<int> degree(num_points, 0);
    //Counting how many times each vertex appears in the constraints
    for (const auto& constraint : additional_constraints) {
        degree[constraint.first]++;
        degree[constraint.second]++;
    }
    //We have to cases for closed constraints. Circle with adge of boundary and classic circle
    bool closed_from_boundary = is_closed_from_boundary(additional_constraints, points, degree, polygon);
    if(closed_from_boundary) {
        return true;
    }
    bool closed_from_circle = is_closed_from_circle(additional_constraints, points, degree, polygon);
    if(closed_from_circle) {
        return true;
    }
    return false;
}

//If exists closed constrains with circle
bool is_closed_from_circle(const vector<pair<int, int>>& additional_constraints, const vector<Point_2>& points, 
                            const vector<int> degree, const Polygon& polygon){
    Point first_point;
    bool found_first_pair = false;
    for (const auto& constraint : additional_constraints) {
        int a = constraint.first;
        int b = constraint.second;
                
        //Check if both vertices have even degrees or the vertex is on boundary
        if ((degree[a] == 2 && degree[b] == 2)) {
            if (!found_first_pair) {
                //Save the first pair
                first_point = points[a];
                found_first_pair = true;
            }
            //Found a closed shape
            if(points[b] == first_point) return true;
            continue;
        }
        if(degree[a] == 1 || degree[b] == 1) {
            found_first_pair = false;
        }
    }
    return false; //No closed shapes found

}

//If exists closed constrains with boundary edge
bool is_closed_from_boundary(const vector<pair<int, int>>& additional_constraints, const vector<Point_2>& points, 
                            const vector<int> degree, const Polygon& polygon){
    bool touch_boundary = false, non_boundary = false;
    for (const auto& constraint : additional_constraints) {
        int a = constraint.first;
        int b = constraint.second;
        
        //Check if both vertices have even degrees or the vertex is on boundary
        if ((!touch_boundary && polygon.has_on_boundary(points[a])) && degree[b] == 2){
            touch_boundary = true;
            //Avoid parallel edges on the boundary to be considered closed_from_boundary
            if(!polygon.has_on_boundary(points[b])) non_boundary = true;
            continue;
        }
        if(touch_boundary && (!polygon.has_on_boundary(points[a]) || !polygon.has_on_boundary(points[b]))) non_boundary = true;
        //If the edge touch again the boundary and touch_boundary == true, we have closed constraints
        if(polygon.has_on_boundary(points[b]) && touch_boundary && non_boundary) return true;
        //Restart
        if(degree[b] == 1) {
            touch_boundary = false;
            non_boundary = false;
        }
    }
    return false;                        

}

bool boundary_straight_lines(const Polygon& polygon){
    if (polygon.size() < 2) return true; //A polygon with fewer than 2 points has no edges to check
    
    //Iterate through the edges of the polygon
    auto prev = polygon.vertices_end() - 1; //Last vertex
    for (auto curr = polygon.vertices_begin(); curr != polygon.vertices_end(); ++curr) {
        //Check if the edge (prev, curr) is axis-aligned
        if (prev->x() != curr->x() && prev->y() != curr->y()) {
            return false; //Edge is neither vertical nor horizontal
        }
        prev = curr; //Move to the next edge
    }

    return true; //All edges are axis-aligned
}

void method_output(const vector<int> count_steiners, std_string method_name, const std_string& name_of_instance, 
                    const int num_steiners, const int init_num_obtuses, const int num_obtuses, bool randomization, 
                    vector<Point_2>& random_steiners, const double rate_of_convergence, double Energy, vector<int> subset,
                    std_string category){
    ofstream outFile("output_simple-polygon-exterior.md", std::ios::app); //Open file for writing

    if (!outFile) {
        cerr<<"Error: Could not open output_best_instances.md for writing!"<<endl;
        return;
    }

    outFile<<"## "<<method_name<<" method, with number of steiners: "<<num_steiners<<" for the instance "<<
            name_of_instance<<" category: "<<category<<endl;
    outFile<<"Number of obtuses: "<<num_obtuses<<endl;
    outFile<<"Used Random steiner: "<<(randomization ? "true" : "false")<<endl;
    if(randomization){
        outFile<<"The Random steiner inserted at: ";
        for (const auto& p : random_steiners) {
            outFile<<fixed<<setprecision(2)<<p<<" => ";
        }
        outFile<<endl;
    }
    outFile<<endl;
    for(int a = 0; a < count_steiners.size(); a++){
        if(a == 0) outFile<<"Circumcenter was selected: "<<count_steiners[a]<<" times"<<endl;
        if(a == 1) outFile<<"Midpoint was selected: "<<count_steiners[a]<<" times"<<endl;
        if(a == 2) outFile<<"Projection was selected: "<<count_steiners[a]<<" times"<<endl;
        if(a == 3) outFile<<"Polygon method was selected: "<<count_steiners[a]<<" times"<<endl;
        if(a == 4) outFile<<"Centroid was selected: "<<count_steiners[a]<<" times"<<endl;
        if(a == 5) outFile<<"Random steiner was selected: "<<count_steiners[a]<<" times"<<endl;
    }
    outFile<<endl;
    outFile<<"Rate of convergence: "<<rate_of_convergence<<endl;
    outFile<<"Energy: "<<Energy<<endl;

    outFile<<"Choosing methods: [";
    for (int i = 0; i < subset.size(); ++i) {
        outFile<<subset[i];
        if (i != subset.size() - 1) outFile<<", ";
    }
    outFile<<"]\n";
    
    double success;
    if(init_num_obtuses > 0) success = ((double)num_obtuses/(double)init_num_obtuses)*100;
    outFile<<100-success<<"%"<<" obtuse triangles reduction success"<<endl;
    outFile<<"## "<<endl<<endl;

    outFile.close();
}

//Write in file the category for every input (oprn -closed constraints etc)
void stats_output(const std_string& name_of_instance, const std_string& category){
    ofstream outFile("stats_output.md", std::ios::app); //Open file for writing

    if (!outFile) {
        cerr<<"Error: Could not open stats_output.md for writing!"<<endl;
        return;
    }

    outFile<<"instance : "<<name_of_instance<<" belongs to the category : "<<category<<endl;   
    outFile.close();
}

//Take the minimum distance of centroid, that the circle touch an edge of the face
double compute_bounding_circle_radius(Face_handle face, const Point_2& centroid) {
    double min_distance = numeric_limits<double>::max();

    for (int i = 0; i < 3; ++i) {
        //Get the edge endpoints
        Point_2 p1 = face->vertex(i)->point();
        Point_2 p2 = face->vertex((i + 1) % 3)->point();

        //Project the centroid onto the edge and compute the distance
        double distance = CGAL::to_double(CGAL::squared_distance(centroid, Segment_2(p1, p2)));
        min_distance = min(min_distance, distance);
    }

    return sqrt(min_distance); //Return the actual radius
}

//Insert steiner point around the centroid with gaussian distribution
void insert_steiner_around_centroid(Custom_CDT& custom_cdt, Face_handle& face, Polygon& polygon, Point_2& steiner_around_centroid) {
    double stddev_ratio = 0.33;
    //Compute the centroid of the triangle
    Point_2 centroid = CGAL::centroid(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());

    //Compute the bounding circle's radius
    double max_distance = compute_bounding_circle_radius(face, centroid);

    //Use a fraction of the bounding radius as the standard deviation, determines the spread of the distribution
    double stddev = stddev_ratio * max_distance;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist_x(CGAL::to_double(centroid.x()), stddev);
    std::normal_distribution<> dist_y(CGAL::to_double(centroid.y()), stddev);

    Point_2 random_point;
    Face_handle f;

    do {
        //Generate a random point around the centroid
        double random_x = dist_x(gen);
        double random_y = dist_y(gen);
        random_point = Point_2(random_x, random_y);
        steiner_around_centroid = random_point;
        f = custom_cdt.locate(random_point);
        
        if(!is_point_inside_region(random_point, polygon)) continue;
    } while (!are_faces_equal(f, face)); //Ensure the point is inside or on the face

    if (is_point_inside_region(random_point, polygon)) {
        custom_cdt.insert_no_flip(random_point);
        start_the_flips(custom_cdt, polygon);
    }  
}

//Try to insert insert_steiner_around_centroid (3rd task)
void try_steiner_around_centroid(Custom_CDT& best_cdt, Polygon& polygon, Point_2& random_steiner){

    Face_handle random_oobtuse_face = give_random_obtuse(best_cdt, polygon);
    /*Point p1 = random_oobtuse_face->vertex(0)->point();
    Point p2 = random_oobtuse_face->vertex(1)->point();
    Point p3 = random_oobtuse_face->vertex(2)->point();*/
    Point_2 steiner_temp;
    insert_steiner_around_centroid(best_cdt, random_oobtuse_face, polygon, steiner_temp);
    random_steiner = steiner_temp;
}

vector<vector<int>> generateSubsetsWith2(int start, int end) {
    vector<vector<int>> subsets;
    int n = end - start + 1; //Number of elements in the range
    int subsetCount = pow(2, n); //Total subsets = 2^n

    //Generate subsets using bitmasking
    for (int mask = 0; mask < subsetCount; ++mask) {
        vector<int> subset;
        bool has2 = false;

        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) { //Check if the i-th bit is set
                int num = start + i;
                subset.push_back(num);
                if (num == 2) { //Check if the number 2 is included
                    has2 = true;
                }
            }
        }
        //Add subset only if it includes the number 2
        if (has2) subsets.push_back(subset);
    }
    return subsets;
}