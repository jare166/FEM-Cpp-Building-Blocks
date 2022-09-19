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
    /* 
    INHERITED PROPERTIES:

    public:
        std::vector< std::vector<double> > nodes;
        std::vector< std::vector<int> >    elements;
        int num_NPE = 2; //! Nodes per element
        int num_Nd  = 0; //! Total number of nodes
        int num_El  = 0; //! Total number of elements
        int dim = 2;     //! Dimension of mesh
    */
    
    public:

    C_Mesh_Frame(int dim_in){ dim = dim_in; }

    //! I. ELEMENT CONSTRUCTION
    //!     Construct Individual 1D Element (Multiple sub-elements Between Principle Nodes)
    void construct_elems(std::vector<double> x_start, std::vector<double> x_end, std::vector<int> x_num, int num_El_sub, int& kk_NODE, int& kk_CONN) {
        /*!
        This function constructs 1D elements. They can be used singly or in a space frame.

        INPUT:
        \param x_start    Position and node number of first point used in construction of 1D element.   { (x1, y1, z1) }
        \param x_end      Position and node number of second point used in construction of 1D element.  { (x2, y2, z2) }
        \param x_num      Node numbers. {xN1, xN2}
        \param num_El_sub Number of elements to divide 1D section into.
        \param kk_NODE    Location to begin storing in global nodal vectors; note that this will also be the first new member for connectivity.
        \param kk_CONN    Location to begin storing in global connectivity vectors.
       
        EXAMPLE:
        Why are the end node numbers needed?
        Input: Two principal nodes, 1 and 2, with 3 sub-nodes between them.
            i.  If principal node numbers not provided
                Original: 1--2
                Result:   1--2--3--4--5
                    1 and 5 are /NEW/ principal node numbers.

            ii. If principal node numbers are provided 
                Original: 1--2
                Result:   1--3--4--5--2
                    Principal nodes can retain their numbers, making BC assignment easier.
        */
        
        std::vector<double> sp_x, sp_y, sp_z;

        sp_x = linspace(num_El_sub+1, x_start[0], x_end[0]); 
        sp_y = linspace(num_El_sub+1, x_start[1], x_end[1]); 
        sp_z = linspace(num_El_sub+1, x_start[2], x_end[2]); 

        // Principal node numbers
        int N1_num = x_num[0];
        int N2_num = x_num[1];

        int jj = 0;
        if (num_El_sub == 1) {
            // 1. Single Element: No additional nodes (other than principal nodes) added
            elements[kk_CONN] = {N1_num, N2_num};
            kk_CONN++;  
        }
        else {
            // 2. Multiple Interior Elements: Additional nodes added
            // NOTE: For each step, the first node in a segment is stored. Thus, in the 
            // case that the first node is a principal node, no new node is stored.
            for (int ii = 0; ii < num_El_sub; ii++) {
                if (ii == 0) { 
                    // i. Store Connectivity w/ First Principal Node
                    elements[kk_CONN] = {N1_num, kk_NODE}; 
                    }
                else if (ii == (num_El_sub-1)) { 
                    // ii. Store Connectivity w/ Second Principal Node
                    elements[kk_CONN] = {kk_NODE, N2_num};    

                    // Store FINAL Interior Node
                    nodes[kk_NODE][0] = sp_x[ii]; 
                    nodes[kk_NODE][1] = sp_y[ii]; 
                    nodes[kk_NODE][2] = sp_z[ii]; 
                    kk_NODE++;  
                }
                else { 
                    // iii. Interior Connectivity; No Principal Nodes Involved
                    elements[kk_CONN] = {kk_NODE, (kk_NODE+1)}; 

                    // Store Interior Node
                    nodes[kk_NODE][0] = sp_x[ii]; 
                    nodes[kk_NODE][1] = sp_y[ii]; 
                    nodes[kk_NODE][2] = sp_z[ii]; 
                    kk_NODE++;
                }
                kk_CONN++;
            }
        }

    }
    //!     Construct Frame from Nodal Data and Connectivity
    void construct_frame(std::vector< std::vector<double> > x_NODE, std::vector< std::vector<double> > x_CONN) {
        double len_elem;

        double x1_S, x2_S, x3_S; // Start Coordinates
        double x1_E, x2_E, x3_E; // End Coordinates

        int jj_S, jj_E;       // Start and End PRINCIPAL Indices, Per Nodal Pair
        int kk_NODE, kk_CONN; // Start local STORAGE indices
        int N;                // Number of Elements, Per Nodal Pair

        // Store Principal Nodes
        kk_NODE = x_NODE.size();
        kk_CONN = 0;
        for (int ii = 0; ii < kk_NODE; ii++) {
            nodes[ii][0] = x_NODE[ii][1]; 
            nodes[ii][1] = x_NODE[ii][2]; 
            nodes[ii][2] = x_NODE[ii][3]; 
        }

        // Store Interior Nodes and Connectivity
        for (int ii = 0; ii < x_CONN.size(); ii++) {
            
            // Start and End Indices
            jj_S = x_CONN[ii][0];
            jj_E = x_CONN[ii][1];

            // Number of Elements in Between Nodal Pair
            N = x_CONN[ii][2];

            // Start and End Coordinates
            x1_S = x_NODE[jj_S][1];
            x2_S = x_NODE[jj_S][2];
            x3_S = x_NODE[jj_S][3];

            x1_E = x_NODE[jj_E][1];
            x2_E = x_NODE[jj_E][2];
            x3_E = x_NODE[jj_E][3];

            construct_elems( {x1_S, x2_S, x3_S}, {x1_E, x2_E, x3_E}, {jj_S, jj_E}, N, kk_NODE, kk_CONN );
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
            str_1 >> temp_val_1;
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
        
            // Read values ( Start, End, Num_Elements )
            str_2 >> temp_val_1;
            std::getline(str_2, word, ','); str_2 >> temp_val_2;
            std::getline(str_2, word, ','); str_2 >> temp_val_3;
            
            x_CONN[iter_read][0] = min(temp_val_1, temp_val_2); // Start Node Number
            x_CONN[iter_read][1] = max(temp_val_1, temp_val_2); // End   Node Number
            x_CONN[iter_read][2] = temp_val_3;                  // Number of Elements

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
        elements.resize(num_El, std::vector<int>(2, 0));

        // Add Data to Global Mesh Objects
        construct_frame(x_NODE, x_CONN);

        return;
    }

};

#endif