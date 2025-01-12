#include "includes/utils/functions.h"
#include "includes/utils/extra_graphics2.h"
#include "includes/utils/functions_task1.h"

using namespace boost::json;
using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Point_2 = K::Point_2;
using Vertex_handle = CDT::Vertex_handle;
namespace bj = boost::json;
using boost_string = bj::string;
using std_string = std::string;

int main(int argc, char** argv) {

    bool run_Simulated_Annealing = false, run_Local_Search = false, run_Ant_Colony = false;
    bool has_constraints= false, is_polygon_convex = false, has_closed_constraints = false, has_open_constraints = false;
    bool has_boundary_straight_lines = false, unspecified = false;
    double alpha = 2.2, beta = 0.1, chi = 3.0, psi = 1.0, lamda = 0.5, kappa = 5;
    int L = 1230, batch_size = 5;
    int sum_convex_no_constraints = 0, sum_convex_open = 0, sum_convex_closed = 0, sum_convex_parallel = 0, sum_unspecified_boundary = 0; 
    value jv;

    std_string input_path, output_path;
    //Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (std_string(argv[i]) == "-i" && i + 1 < argc) {
            input_path = argv[++i];
        } else if (std_string(argv[i]) == "-o" && i + 1 < argc) {
            output_path = argv[++i];
        }
    }

    if (input_path.empty() || output_path.empty()) {
        cerr<<"Empty input path or output path."<<endl;
        cout<<"Check this pattern of terminal order: ./opt_triangulation -i /path/to/input.json -o /path/to/output.json"<<endl;
        return 1;
    }
    //Check the names of the test cases in folder tests
    //f.e. ./opt_triangulation -i tests/test_SA.json -o solution_output.json
    //./opt_triangulation -i tests/test_Ants.json -o solution_output.json
    //./opt_triangulation -i tests/test_Local.json -o solution_output.json
    
    read_json(input_path, jv);
    Custom_CDT custom_cdt;
    vector<Point_2> points;
    vector<pair<int, int>> additional_constraints;
    vector<int> region_boundary;
    Polygon polygon;
    Polygon simulated_polygon;
    std_string method, instance_uid;
    bool delaunay = true;

//////////// PHASE 1: INITIALIZATION //////////////////////////////
    
    if (jv.is_object()) {
        const auto& obj = jv.as_object();
        const auto& x_array = obj.at("points_x").as_array();
        const auto& y_array = obj.at("points_y").as_array();
        const auto& boundary_array = obj.at("region_boundary").as_array();
        const auto& constraints_array = obj.at("additional_constraints").as_array();
        method = std_string(obj.at("method").as_string());
        instance_uid = std_string(obj.at("instance_uid").as_string());
        const auto& parameters_obj = obj.at("parameters").as_object();
        delaunay = obj.at("delaunay").as_bool();

        //Output method and delaunay
        //cout<<"method: "<<method<<endl;
        L = parameters_obj.at("L").as_int64();
        
        //Chosen method
        if(method == "local") {
            run_Local_Search = true;
        }
        else if(method == "sa") {
            run_Simulated_Annealing = true;
            alpha = parameters_obj.at("alpha").as_double();
            beta = parameters_obj.at("beta").as_double();
            //How many "bad" steiners we accept to insert, until we will try again to add steiners in the best_cdt
            batch_size = parameters_obj.at("batch_size").as_int64();
        }
        else if(method == "ant"){
            run_Ant_Colony = true;
            
            alpha = parameters_obj.at("alpha").as_double();
            beta = parameters_obj.at("beta").as_double();
            lamda = parameters_obj.at("lambda").as_double();
            chi = parameters_obj.at("xi").as_double();
            psi = parameters_obj.at("psi").as_double();
            kappa = parameters_obj.at("kappa").as_int64();
        }
        else {
            cerr<<"Error: wrong method"<<endl;
            return 0;
        }

        for (int i = 0; i < x_array.size(); ++i) {
            double x = x_array[i].is_double() ? x_array[i].as_double() : static_cast<double>(x_array[i].as_int64());
            double y = y_array[i].is_double() ? y_array[i].as_double() : static_cast<double>(y_array[i].as_int64());
            points.emplace_back(x, y);
        }

        for (const auto& idx : boundary_array) {
            region_boundary.push_back(idx.as_int64());
        }
        
        //Add the additional constraints in vector
        for (const auto& constraint : constraints_array) {
            int idx1 = constraint.as_array()[0].as_int64();
            int idx2 = constraint.as_array()[1].as_int64();
            if (idx1 < points.size() && idx2 < points.size()) {
                additional_constraints.emplace_back(idx1, idx2);
            }
        }

        //Create a polygon from region boundary
        for (int index : region_boundary) {
            polygon.push_back(points[index]);
        }

         // Insert region boundary as constraints
        for (size_t i = 0; i < region_boundary.size(); ++i) {
            int idx1 = region_boundary[i];
            int idx2 = region_boundary[(i + 1) % region_boundary.size()]; // Wrap around to form a loop
            custom_cdt.insert_constraint(points[idx1], points[idx2]);
        }

        //Make the cdt
        for (const auto& point : points) {
            custom_cdt.insert(point);
        }

        /*for (auto fit = custom_cdt.finite_faces_begin(); fit != custom_cdt.finite_faces_end(); ++fit) {
            for (int i = 0; i < 3; ++i) { // Each triangle has 3 edges
                Point_2 p1 = fit->vertex((i + 1) % 3)->point();
                Point_2 p2 = fit->vertex((i + 2) % 3)->point();

                // Check if the edge is constrained and on the boundary
                if (custom_cdt.is_constrained(CDT::Edge(fit, i)) && is_edge_on_boundary(p1, p2, polygon)) {
                    // Remove the constraint from the edge
                    custom_cdt.remove_constraint(fit, i);
                }
            }
        }*/

        //Insert additional constraints
        for (const auto& constraint : additional_constraints) {
            custom_cdt.insert_constraint(points[constraint.first], points[constraint.second]);
        }

        /*********************************/
        //Check if the polygon is convex (3rd Task)
        if(polygon.is_convex()) is_polygon_convex = true;
        //Check if the polygon (boundary) has straight lines (3rd Task)
        if(boundary_straight_lines(polygon)) has_boundary_straight_lines = true;
        if(!additional_constraints.empty()) has_constraints = true;

        //cout<<"the instance has constraints: "<<has_constraints<<" and is convex: "<<is_polygon_convex<<endl;
        //cout<<"has_boundary_straight_lines: "<<has_boundary_straight_lines<<endl;

        //3rd Task
        if(!has_constraints && !is_polygon_convex && !has_boundary_straight_lines) {
            unspecified = true;
            //cout<<"AkANoNiStO"<<endl;
        }
        /*********************************/        
        //Check if the instance has opened or closed constraints
        if(has_constraints){
            if(are_constraints_closed(additional_constraints, points.size(), points, polygon)) {
                has_closed_constraints = true;
                //cout<<"has_closed_constraints: "<<has_closed_constraints<<endl;
            }
            //If has not closed constraints and we have constraint, so we have open constraints
            else{
                has_open_constraints = true;
                //cout<<"has_open_constraints: "<<has_open_constraints<<endl;
            }
        }
    }
    else {
        cerr<<"Jv is not object: safe exit"<<endl;
        return 0;
    }

    std_string category = "Null";
    if(is_polygon_convex && !has_constraints){
        category = "CONVEX_NO_CONSTRAINTS";
        
        //cout<<"instance : "<<instance_uid<<" belongs to the category : CONVEX_NO_CONSTRAINTS"<<endl;
        sum_convex_no_constraints++;
    }
    if(is_polygon_convex && has_open_constraints){
        category = "CONVEX_OPEN_CONSTRAINTS";
        
        //cout<<"instance : "<<instance_uid<<" belongs to the category : CONVEX_OPEN_CONSTRAINTS"<<endl;
        sum_convex_open++;
    }
    if(is_polygon_convex && has_closed_constraints){
        category = "CONVEX_CLOSED_CONSTRAINTS";
        
        //cout<<"instance : "<<instance_uid<<" belongs to the category : CONVEX_CLOSED_CONSTRAINTS"<<endl;
        sum_convex_closed++;
    }
    if(!is_polygon_convex && has_boundary_straight_lines){
        category = "NOT_CONVEX_PARALLEL_N0_CONSTRAINTS";
        
        //cout<<"instance : "<<instance_uid<<" belongs to the category : NOT_CONVEX_PARALLEL_N0_CONSTRAINTS"<<endl;
        sum_convex_parallel++;
    }
    if(!is_polygon_convex && unspecified){
        category = "UNSPECIFIED_BOUNDARY";
        
        //cout<<"instance : "<<instance_uid<<" belongs to the category : UNSPECIFIED_BOUNDARY"<<endl;
        sum_unspecified_boundary++;
    }
    //stats_output(instance_uid, category);
    //////////// PHASE 2: FLIPS & STEINER POINTS //////////////////////////////
    int obtuses_faces = count_obtuse_triangles(custom_cdt, polygon);
    int init_obtuse_faces = obtuses_faces;
    int initial_vertexes = custom_cdt.number_of_vertices();
    cout<<"Initial number of obtuses: "<<init_obtuse_faces<<endl;
    cout<<"Initial number of vertexes: "<<custom_cdt.number_of_vertices()<<endl;
    //CGAL::draw(custom_cdt);
    double success;
    bool randomization = false;
    Custom_CDT simulated_cdt = custom_cdt;
    //Run task1 if delaunay parameter is false
    if(!delaunay) {
        cout<<"**Run task1**"<<endl;
        run_task1(simulated_cdt, polygon);
        obtuses_faces = count_obtuse_triangles(simulated_cdt, polygon);
        cout<<"Number of obtuses after task 1: "<<obtuses_faces<<endl;
        cout<<"Sum of steiners after task 1: "<<count_vertices(simulated_cdt) - initial_vertexes<<endl;
        if(init_obtuse_faces > 0) success = ((double)obtuses_faces/(double)init_obtuse_faces)*100;
        cout<<100-success<<"%"<<" obtuse triangles reduction success after task 1"<<endl;
    }

    simulated_polygon = polygon;
    
    //Flips
    start_the_flips(simulated_cdt, simulated_polygon);
    //Local Search
    if(run_Local_Search){
        cout<<"Local Search is starting.."<<endl;
        local_search(simulated_cdt, simulated_polygon, L, instance_uid, randomization);
        cout <<"**Number of Obtuses after from Local Search: "<<count_obtuse_triangles(simulated_cdt, simulated_polygon)<<" **"<<endl;
    }

    //SA
    if(run_Simulated_Annealing){
        cout<<"Simulated Annealing is starting.. "<<endl;
        simulated_annealing(simulated_cdt, simulated_polygon, L, alpha, beta, batch_size, instance_uid);
        cout <<"**Number of Obtuses after from Simulated Annealing: "<<count_obtuse_triangles(simulated_cdt, simulated_polygon)<<" **"<<endl;
    }
    //Ant Colony
    if(run_Ant_Colony){
        cout<<"Ant Colony is starting.. "<<endl;
        ant_colony(simulated_cdt, simulated_polygon, alpha, beta, chi, psi, lamda , L, kappa, instance_uid);
        cout <<"**Number of Obtuses after from Ant Colony: "<<count_obtuse_triangles(simulated_cdt, simulated_polygon)<<" **"<<endl;
    }

    obtuses_faces = count_obtuse_triangles(simulated_cdt, simulated_polygon);
    cout<<"Final obtuses faces: "<<obtuses_faces<<endl;
    cout<<"Sum of steiners: "<<simulated_cdt.number_of_vertices() - initial_vertexes<<endl;
    cout<<"Final number of vertexes: "<<simulated_cdt.number_of_vertices()<<endl;
    if(init_obtuse_faces > 0) success = ((double)obtuses_faces/(double)init_obtuse_faces)*100;
    cout<<100-success<<"%"<<" obtuse triangles reduction success"<<endl;
    cout<<"Final form of Custom CDT "<<endl;
    //CGAL::draw(simulated_cdt);
    //print_polygon_edges(simulated_polygon);
    
    /*Print from Extra_Graphics*/
    //Calculate min_y and max_y
    double min_y = numeric_limits<double>::max();
    double max_y = numeric_limits<double>::lowest();

    for (const auto& point : points) {
        double y = CGAL::to_double(point.y());
        min_y = min(min_y, y);
        max_y = max(max_y, y);
    }
    QApplication app(argc, argv);
    CDTGraphicsView view(simulated_cdt, polygon);
    view.setRenderHint(QPainter::Antialiasing);
    view.setWindowTitle("Delaunay Triangulation with Point Coordinates");
    view.resize(1000, 1000);
    //Center and zoom the view
    view.fitInView(view.scene()->sceneRect(), Qt::KeepAspectRatio);
    view.scale(1.5, 1.5);

    double centerY = (min_y + max_y) / 2; //Calculate the center y position based on your points
    view.verticalScrollBar()->setValue(centerY);
    view.show();
    //////////// PHASE 3: JSON FILE OUTPUT //////////////////////////////

    output(jv, simulated_cdt, points, obtuses_faces, output_path, randomization);

    return app.exec();
    //return 0;
}
