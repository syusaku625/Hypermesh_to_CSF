#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include<list>
#include <boost/algorithm/string.hpp>
#include<boost/progress.hpp>
#include<set>
#include<cmath>
#include <iomanip>
using namespace std;

void export_element_part(string filename, vector<vector<int>> element)
{
    ofstream ofs(filename);
    for(int i=0; i<element.size(); i++){
        for(int j=0; j<element[i].size(); j++){
            ofs << element[i][j] << " ";
        }
        ofs << endl;
    }
    ofs.close();
}

void export_element_id(string filename, vector<vector<int>> element)
{
    boost::progress_display show_progress(element.size());
    ofstream ofs(filename);
    for(int i=0; i<element.size(); i++){
        ++show_progress;
        if(element[i].size()==4) ofs << 10 << endl;
        if(element[i].size()==6) ofs << 13 << endl;
    }
    ofs.close();
}

void input_surface_node_and_element_construct(string filename, vector<vector<double>> &x, vector<vector<int>> &element)
{
    ifstream ifs(filename);
    string str;
    for(int i=0; i<5; i++){
        getline(ifs,str);
    }
    int numOfnode, numOfelem;
    istringstream ss(str);
    for(int i=0; i<3; i++){
        getline(ss, str, ' ');
        if(i==2) numOfnode=stoi(str);
    }    
    for(int i=0; i<4; i++){
        getline(ifs,str);
    }
    istringstream ss2(str);
    for(int i=0; i<3; i++){
        getline(ss2, str, ' ');
        if(i==2) numOfelem=stoi(str);
    }
    for(int i=0; i<2; i++){
        getline(ifs,str);
    }
    for(int i=0; i<numOfnode; i++){
        getline(ifs,str);
        istringstream ss1(str);
        vector<double> tmp_x;
        for(int j=0; j<3; j++){
            getline(ss1, str, ' ');
            tmp_x.emplace_back(stod(str));
        }
        x.emplace_back(tmp_x);
    }
    for(int i=0; i<numOfelem; i++){
        getline(ifs,str);
        istringstream ss1(str);
        vector<int> tmp_x;
        getline(ss1,str,' ');
        int loop=stoi(str);
        for(int j=0; j<loop; j++){
            getline(ss1, str, ' ');
            tmp_x.emplace_back(stoi(str));
        }
        element.push_back(tmp_x);
    }
}

void export_boundary_surface_node_velocity(string filename, vector<vector<double>> x, vector<vector<double>> surface_x)
{
    boost::progress_display show_progress(x.size());
    ofstream ofs(filename);
    for(int i=0; i<x.size(); i++){
        ++show_progress;
        for(int j=0; j<surface_x.size(); j++){
            if(fabs(x[i][0]-surface_x[j][0])<0.0001 && fabs(x[i][1]-surface_x[j][1])<0.0001 && fabs(x[i][2]-surface_x[j][2])<0.0001){
                ofs << i << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
                break;
            }
        }
    }   
    ofs.close();
}

void export_boundary_surface_node_pressure(string filename, vector<vector<double>> x, vector<vector<double>> surface_x)
{
    boost::progress_display show_progress(x.size());
    ofstream ofs(filename);
    for(int i=0; i<x.size(); i++){
        ++show_progress;
        for(int j=0; j<surface_x.size(); j++){
            if(fabs(x[i][0]-surface_x[j][0])<0.0001 && fabs(x[i][1]-surface_x[j][1])<0.0001 && fabs(x[i][2]-surface_x[j][2])<0.0001){
                ofs << i << " " << 2 << " " << 0 << endl;
                break;
            }
        }
    }   
    ofs.close();
}

void export_boundary_surface_node_disp(string filename, vector<vector<double>> x, vector<vector<double>> surface_x)
{
    boost::progress_display show_progress(x.size());
    ofstream ofs(filename);
    for(int i=0; i<x.size(); i++){
        ++show_progress;
        for(int j=0; j<surface_x.size(); j++){
            if(fabs(x[i][0]-surface_x[j][0])<0.0001 && fabs(x[i][1]-surface_x[j][1])<0.0001 && fabs(x[i][2]-surface_x[j][2])<0.0001){
                ofs << i << " " << 0 << " " << 0 << " " << 0 << endl;
                break;
            }
        }
    }   
    ofs.close();
}

void export_boundary_surface_node_disp_lock(string filename, vector<vector<double>> x, vector<vector<double>> surface_x, vector<vector<double>> surface_x_full)
{
    set<int> surface_part_node;
    vector<int> surface_full_node;
    boost::progress_display show_progress(x.size());
    ofstream ofs(filename);
    for(int i=0; i<x.size(); i++){
        for(int j=0; j<surface_x.size(); j++){
            if(fabs(x[i][0]-surface_x[j][0])<0.0001 && fabs(x[i][1]-surface_x[j][1])<0.0001 && fabs(x[i][2]-surface_x[j][2])<0.0001){
                surface_part_node.insert(i);
                break;
            }
        }
    }
    for(int i=0; i<x.size(); i++){
        ++show_progress;
        for(int j=0; j<surface_x_full.size(); j++){
            if(fabs(x[i][0]-surface_x_full[j][0])<0.0001 && fabs(x[i][1]-surface_x_full[j][1])<0.0001 && fabs(x[i][2]-surface_x_full[j][2])<0.0001){
                surface_full_node.push_back(i);
                break;
            }
        }
    }

    for(int i=0; i<surface_full_node.size(); i++){
        if(surface_part_node.find(surface_full_node[i])==surface_part_node.end()){
            ofs << surface_full_node[i] << " " << 0 << " " << 0 << " " << 0 << endl;
        }
    }
    ofs.close();
}

vector<vector<double>> finite_element_surface_normal_vector(vector<vector<int>> surface_element_data, vector<vector<double>> node_coordinate_data)
{
    vector<vector<double>> normal_vector;
    normal_vector.resize(surface_element_data.size());

    #pragma omp parallel for
    for (int i = 0; i < normal_vector.size(); i++){
        normal_vector[i].resize(3);
    }
    cout << "calculate normal vector" << endl;
    boost::progress_display show_progress(surface_element_data.size());
    for (int i = 0; i < surface_element_data.size(); i++)
    {
        ++show_progress;
        vector<vector<double>> O(3, vector<double>(3));
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                O[j][k] = node_coordinate_data[surface_element_data[i][j]][k];
            }
        }
        vector<double> OA(3), OB(3);
        #pragma omp parallel for
        for (int j = 0; j < 3; j++){
            OA[j] = O[1][j] - O[0][j];
            OB[j] = O[2][j] - O[0][j];
        }

        normal_vector[i][0] = OA[1] * OB[2] - OA[2] * OB[1];
        normal_vector[i][1] = OA[2] * OB[0] - OA[0] * OB[2];
        normal_vector[i][2] = OA[0] * OB[1] - OA[1] * OB[0];

        double vec_size=sqrt(pow(normal_vector[i][0],2.0)+pow(normal_vector[i][1],2.0)+pow(normal_vector[i][2],2.0));
        #pragma omp parallel for
        for (int j = 0; j < 3; j++){
            normal_vector[i][j] /= -vec_size;
        }

    }
    return normal_vector;
}

void export_vector_double(string filename, vector<vector<double>> normal_vector)
{
    ofstream ofs(filename);
    for(int i=0; i<normal_vector.size(); i++){
        for(int j=0; j<3; j++){
            ofs << normal_vector[i][j] << " ";
        }
        ofs << endl;
    }
    ofs.close();
}

void export_vector_int(string filename, vector<vector<int>> normal_vector)
{
    ofstream ofs(filename);
    for(int i=0; i<normal_vector.size(); i++){
        for(int j=0; j<3; j++){
            ofs << normal_vector[i][j] << " ";
        }
        ofs << endl;
    }
    ofs.close();
}

void return_area_value(vector<vector<int>> element_node, vector<vector<double>> node_coordinate)
{
    ofstream ofs("area.dat");
    double area_sum = 0.0;
    for (int i = 0; i < element_node.size(); i++){
        vector<vector<double>> x(3, vector<double>(2));
        for (int j = 0; j < 3; j++){
            int node = element_node[i][j];
            x[j][0] = node_coordinate[node][0];
            x[j][1] = node_coordinate[node][1];
        }
        double delta2 = x[0][0] * (x[1][1] - x[2][1]) + x[1][0] * (x[2][1] - x[0][1]) + x[2][0] * (x[0][1] - x[1][1]);
        double delta = delta2 / 2.0;
        area_sum += fabs(delta);
    }
    ofs << area_sum << endl;
    ofs.close();
}

int main(int argc, char *argv[])
{
    string file = argv[1];
    string filename=file+".vtk";
    string str;
    ifstream ifs(filename);
    for(int i=0; i<5; i++){
        getline(ifs,str);
    }
    istringstream iss(str);
    string s;
    int count = 0;
    int numOfNode;
    while (iss >> s) {
        if(count==1){
            numOfNode=stod(s);
        }
        count++;
    }
    vector<vector<double>> x;

    for(int i=0; i<numOfNode/3; i++){
        getline(ifs,str);
        istringstream node_str(str);
        vector<double> node_buf;
        while (getline(node_str, str, ' ')){
            try
            {
                node_buf.push_back(stod(str));
            }
            catch(const std::exception& e)
            {
                continue;
            }
        }
        for(int j=0; j<node_buf.size()/3; j++){
            vector<double> x_tmp;
            x_tmp.push_back(node_buf[j*3]);
            x_tmp.push_back(node_buf[j*3+1]);
            x_tmp.push_back(node_buf[j*3+2]);
            x.push_back(x_tmp);
        }
    }

    while(getline(ifs,str)){
        list<string> list_string;
        boost::split(list_string, str, boost::is_any_of(" "));
        if(*(list_string.begin())=="CONNECTIVITY") break;
    }
    list<int> cell_node;
    
    while(getline(ifs,str)){
        bool flag = false;
        vector<std::string> result;
        istringstream stream(str);
        string tmp_cell_id;
        while (getline(stream, tmp_cell_id, ' ')){
            count = 0;
            if(count==0){
                if(tmp_cell_id == "CELL_TYPES"){
                    flag = true;
                    break;
                } 
            }
            count++;
            try
            {
                cell_node.push_back(stoi(tmp_cell_id));

            }
            catch(const std::exception& e)
            {
                continue;
            }
        }
        if(flag) break;
    }
    istringstream cell_num_read(str);
    int numOfCell;
    string tmp_cell_type;
    count = 0;
    while (getline(cell_num_read, tmp_cell_type, ' ')){
        if(count==1){
            numOfCell = stoi(tmp_cell_type);
        } 
        count++;
    }
    vector<int> cell_type;
    for(int i=0; i<numOfCell; i++){
        getline(ifs,str);
        cell_type.push_back(stoi(str));
    }
    vector<vector<int>> element;
    for(int i=0; i<cell_type.size(); i++){
        vector<int> tmp_element;
        int cell_construct_number;
        if(cell_type[i]==10) cell_construct_number = 4;
        else if(cell_type[i]==13) cell_construct_number = 6;
        for(int j=0; j<cell_construct_number; j++){
            try
            {
                tmp_element.push_back(*(cell_node.begin()));
                cell_node.pop_front();
            }
            catch(const std::exception& e)
            {
                continue;
            }
            
        }
        element.push_back(tmp_element);
    }
    ifs.close();

    
    string node_file = "node.dat";
    string element_file = "element.dat";
    string elementType_file = "elementType.dat";
    export_vector_double(node_file,x);
    export_vector_int(element_file,element);
    export_element_id(elementType_file,element);
    
    cout << "surface_process" << endl;

    vector<vector<double>> surface_x_velocity;  
    vector<vector<double>> surface_x_pressure;  
    vector<vector<double>> surface_x_disp;  
    vector<vector<double>> surface_x_full;
    vector<vector<int>> surface_element_velocity;
    vector<vector<int>> surface_element_pressure;
    vector<vector<int>> surface_element_disp;
    vector<vector<int>> surface_element_full;

    input_surface_node_and_element_construct("surface_boundary_v.ply",surface_x_velocity,surface_element_velocity);
    input_surface_node_and_element_construct("surface_boundary_p.ply",surface_x_pressure,surface_element_pressure);
    input_surface_node_and_element_construct("surface_boundary_disp.ply",surface_x_disp,surface_element_disp);
    input_surface_node_and_element_construct("surface_full.ply",surface_x_full,surface_element_full);

    cout << "export_boundary_v" << endl;
    export_boundary_surface_node_velocity("boundary_v2.dat", x, surface_x_velocity);
    cout << "export_boundary_p" << endl;
    export_boundary_surface_node_pressure("boundary_p.dat", x, surface_x_pressure);
    cout << "export_boundary_disp" << endl;
    export_boundary_surface_node_disp("boundary_disp.dat", x, surface_x_disp);
    export_boundary_surface_node_disp_lock("boundary_disp2.dat", x, surface_x_disp, surface_x_full);

    cout << "export_wall_vector" << endl;
    vector<vector<double>> normal_vector=finite_element_surface_normal_vector(surface_element_disp, surface_x_disp);
    export_vector_double("wall_vector.dat",normal_vector);

    return_area_value(surface_element_disp,surface_x_disp);
}