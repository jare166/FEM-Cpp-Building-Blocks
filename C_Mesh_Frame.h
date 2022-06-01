#ifndef C_MESH_FRAME_H
#define C_MESH_FRAME_H

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "f_MiscellaneousFunctions.h"
#include "C_Mesh.h"

// DECLARATIONS
template<typename T> std::vector<T> linspace(int num_in, T start_in, T end_in);

//! This class is used to store data describing a mesh of 1, 2, or 3 dimensions.
class C_Mesh_Frame : public C_Mesh{
    public:
        std::vector< std::vector<double> > nodes;
        std::vector< std::vector<int> >    elements;
        int num_NPE = 2; //! Nodes per element
        int num_Nd  = 0; //! Total number of nodes
        int num_El  = 0; //! Total number of elements
        int dim = 2;     //! Dimension of mesh

    C_Mesh_Frame() { 
        nodes    = { {0,0}, {1,1} };
        elements = { {0,1} };

        num_Nd = 2;
        num_El = num_Nd - 1;
        dim = 2;
    }

    C_Mesh_Frame(int dim_in){ dim = dim_in; }

    C_Mesh_Frame(double l, int num_Nd_in){
        // 1D Default Initializer
        num_Nd = num_Nd_in;
        num_El = num_Nd-1;

        dim = 1;

        nodes.resize(num_Nd, std::vector<double>(3, 0));
        elements.resize(num_Nd, std::vector<int>(2, 0));

        std::vector<double> x_start = {0, 0, 0};
        std::vector<double> x_end   = {l, 0, 0};

        construct_elems(x_start, x_end, num_Nd, 0);
    }

    //! I. ELEMENT CONSTRUCTION
    //!     Construct Individual 1D Element (Multiple sub-elements Between Principle Nodes)
    void construct_elems(std::vector<double> x_start, std::vector<double> x_end, int num_Nd_sub, int ii_s) {
        std::vector<double> sp_x, sp_y, sp_z;

        sp_x = linspace(num_Nd_sub, x_start[0], x_end[0]); 
        sp_y = linspace(num_Nd_sub, x_start[1], x_end[1]); 
        sp_z = linspace(num_Nd_sub, x_start[2], x_end[2]); 

        int ii = 0;
        for (ii = 0; ii <= (num_Nd_sub-1); ii++) {
            elements[ii_s+ii] = {ii, ii+1};

            nodes[ii_s+ii][0] = sp_x[ii]; 
            nodes[ii_s+ii][1] = sp_y[ii]; 
            nodes[ii_s+ii][2] = sp_z[ii]; 
        }

    }
    //!     Construct Frame from Nodal Data and Connectivity
    void construct_frame(std::vector< std::vector<double> > x_NODE, std::vector< std::vector<double> > x_CONN) {
        double len_elem;

        double x1_S, x2_S, x3_S; // Start Coordinates
        double x1_E, x2_E, x3_E; // End Coordinates

        int    jj_S, jj_E; // Start and End Indices, Per Nodal Pair
        int    N;          // Number of Elements, Per Nodal Pair

        for (int ii = 0; ii < x_CONN.size(); ii++) {
            
            // Start and End Indices
            jj_S = x_CONN[ii][0];
            jj_E = x_CONN[ii][1];

            // Number of Elements in Between Nodal Pair
            N  = x_CONN[ii][2];

            // Start and End Coordinates
            x1_S = x_NODE[jj_S][1];
            x2_S = x_NODE[jj_S][2];
            x3_S = x_NODE[jj_S][3];

            x1_E = x_NODE[jj_E][1];
            x2_E = x_NODE[jj_E][2];
            x3_E = x_NODE[jj_E][3];

            len_elem = norm( (x1_E-x1_S), (x2_E-x2_S), (x3_E-x3_S) );
            construct_elems({x1_S, x2_S, x3_S}, {x1_E, x2_E, x3_E}, N, jj_S);
        }
    }

    //! II. READ-INITIALIZE FROM FILE
    //!     Read Frame Data from File and Initialize Object
    void read_frame (std::string filePath_node, std::string filePath_conn) {
        
        if ( filePath_node.empty() | filePath_conn.empty() ) { return; }

        // I. READ NODAL FILE
        std::fstream fID_1;
        fID_1.open(filePath_node);

        std::string line, word;
        double temp_val_1, temp_val_2, temp_val_3, temp_val_4;
        int num_pr_nodes; // Number of Principal Nodes
        int num_connect;  // Total Number of Connections

        // Get number of lines: NODAL FILE
        std::getline(fID_1, line);
        std::stringstream str_1(line);

        std::getline(str_1, word, ',');
        str_1 >> num_pr_nodes;

        // Assign nodal storage vectors
        std::vector< std::vector<double> >  x_NODE( num_pr_nodes, std::vector<double>(4, 0.0) );
        //  Details: vector of num_pr_nodes, 4-column vectors, each w/ default value of 0.0.

        // Read Successive Lines for Data
        int iter_read = 0;
        while (std::getline(fID_1, line)) {
            // Read from next line
            std::stringstream str_1(line);
        
            // Read values ( Node, (x, y, z) )
            std::getline(str_1, word, ','); str_1 >> temp_val_1;
            std::getline(str_1, word, ','); str_1 >> temp_val_2;
            std::getline(str_1, word, ','); str_1 >> temp_val_3;
            std::getline(str_1, word, ','); str_1 >> temp_val_4;
            
            x_NODE[iter_read][0] = temp_val_1; // Node Number
            x_NODE[iter_read][1] = temp_val_2; // x-coord
            x_NODE[iter_read][2] = temp_val_3; // y-coord
            x_NODE[iter_read][3] = temp_val_4; // z-coord

            // Exit if requested number of lines has been read
            iter_read++;
            if (iter_read == (num_pr_nodes)) { break; }
        }
        fID_1.close();


        // II. READ CONNECTIVITY FILE
        std::fstream fID_2;
        fID_2.open(filePath_conn);

        // Get number of lines: CONNECTIVITY FILE
        std::getline(fID_2, line);
        std::stringstream str_2(line);

        std::getline(str_2, word, ',');
        str_2 >> num_connect;

        // Assign nodal storage vectors
        std::vector< std::vector<double> >  x_CONN( num_connect, std::vector<double>(3, 0.0) );
        //  Details: vector of num_connect, 3-column vectors, each w/ default value of 0.0.

        // Read Successive Lines for Data
        iter_read = 0;
        while (std::getline(fID_2, line)) {
            // Read from next line
            std::stringstream str_2(line);
        
            // Read values ( Node, (x, y, z) )
            std::getline(str_2, word, ','); str_2 >> temp_val_1;
            std::getline(str_2, word, ','); str_2 >> temp_val_2;
            std::getline(str_2, word, ','); str_2 >> temp_val_3;
            
            x_CONN[iter_read][0] = min(temp_val_1, temp_val_2); // Start Node Number
            x_CONN[iter_read][1] = max(temp_val_1, temp_val_2); // End   Node Number
            x_CONN[iter_read][2] = temp_val_3; // Number of Elements

            // Exit if requested number of lines has been read
            iter_read++;
            if (iter_read == (num_connect)) { break; }
        }
        fID_2.close();

        // USER SHOULD CHECK TO ENSURE THERE ARE NO REPEATS


        // III. POPULATE OBJECT WITH INPUT DATA
        // Determine total number of nodes and resize vectors
        num_Nd = num_pr_nodes;
        num_El = 0;
        for (int ii = 0; ii < num_connect; ii++) { 
            num_Nd += (x_CONN[ii][2]-1); 
            num_El +=  x_CONN[ii][2];
        }

        nodes.resize(num_Nd, std::vector<double>(3, 0));
        elements.resize(num_Nd, std::vector<int>(2, 0));

        // Add Data to Global Mesh Objects
        construct_frame(x_NODE, x_CONN);

        return;
    }
};

#endif