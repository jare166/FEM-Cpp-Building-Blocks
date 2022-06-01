#ifndef F_MISCELLANEOUSFUNCTIONS_H
#define F_MISCELLANEOUSFUNCTIONS_H

#include <math.h>
#include <vector>

// Function Declarations
template <typename T> std::vector<T> linspace  (int num_in, T start_in, T end_in);
std::vector<int>                     intspace  (int start_in, int end_in);
template <typename T> T     signum    (T x);
template <typename T> T     Heaviside (T x);
template <typename T> T     norm      (T x, T y);
template <typename T> T     norm      (T x, T y, T z);
template <typename T> T     norm      (std::vector<T> x_v);
template <typename T> T     min       (T x, T y);
template <typename T> T     min       (std::vector<T> x_v);
template <typename T> T     min       (std::vector<T> x_v, int& ind_out);
template <typename T> T     max       (T x, T y);
template <typename T> T     max       (std::vector<T> x_v);
template <typename T> T     max       (std::vector<T> x_v, int& ind_out);

// i. Special Vector Initializations
template <typename T> std::vector<T> linspace(int num_in, T start_in, T end_in)
{
    //! Generate linear spaced STL-style vectors.
    std::vector<T> linspaced;

    T start = static_cast<T>(start_in);
    T end   = static_cast<T>(end_in);
    T num   = static_cast<T>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    T delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i) {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); 

    return linspaced;
}

std::vector<int> intspace(int start_in, int end_in)
{
    //! Generate integer spaced STL-style vectors.
    
    // Create a vector varying from start int to end int
    int num_ind = std::abs(end_in - start_in);
    int inc_sig = signum(end_in - start_in);

    // Initialize vector
    std::vector<int> intspaced(num_ind);
    for(int ii = 0; ii < num_ind; ii++) { intspaced[ii] = start_in + ii; }

    return intspaced;
}

// ii. General Mathematical Functions
template <typename T> T signum(T x) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> T Heaviside (T x) {
    if (x < T(0)) { return T(0); }
    else          { return x;    }
}

//  COMPUTE VECTOR NORM
template <typename T> T norm (T x, T y) {
    return sqrt( pow(x,2) + pow(y,2) );
}

template <typename T> T norm (T x, T y, T z) {
    return sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
}

template <typename T> T norm (std::vector<T> x_v) {
    T x_sum = T(0);

    for (int ii = 0; ii < x_v.size(); ii++) {
        x_sum += pow(x_v[ii], 2);
    }

    return sqrt(x_sum);
}

//  COMPUTE MINIMA
template <typename T> T min (T x, T y) {
    return x<y ? x : y;
}
//      ii. Test variable number of input values
template <typename T> T min (std::vector<T> x_v) {
    T x_out = x_v[0];

    for (int ii = 1; ii < x_v.size(); ii++) 
    { 
        x_out = min(x_out, x_v[ii]); 
    }

    return x_out;
}
//      iii. Return indices also
template <typename T> T min (std::vector<T> x_v, int& ind_out) {
    T x_out = x_v[0];

    for (int ii = 1; ii < x_v.size(); ii++) 
    { 
        x_out = min(x_out, x_v[ii]); 
        ind_out = ii;
    }
    

    return x_out;
}

//  COMPUTE MAXIMA
template <typename T> T max (T x, T y) {
    return x>y ? x : y;
}
//      ii. Test variable number of input values
template <typename T> T max (std::vector<T> x_v) {
    T x_out = x_v[0];

    for (int ii = 1; ii < x_v.size(); ii++) 
    { 
        x_out = max(x_out, x_v[ii]); 
    }
    return x_out;
}
//      iii. Return indices also
template <typename T> T max (std::vector<T> x_v, int& ind_out) {
    T x_out = x_v[0];

    for (int ii = 1; ii < x_v.size(); ii++) 
    { 
        x_out = max(x_out, x_v[ii]); 
        ind_out = ii;
    }
    

    return x_out;
}

#endif