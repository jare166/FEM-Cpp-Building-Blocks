#ifndef F_MISCELLANEOUSFUNCTIONS_H
#define F_MISCELLANEOUSFUNCTIONS_H

#include <math.h>
#include <vector>

// Function Declarations
template <typename T> std::vector<T> linspace  (int num_in, T start_in, T end_in);
std::vector<int>                     intspace  (int start_in, int end_in);
template <typename T> int            signum    (T x);
template <typename T> double         Heaviside (T x);


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
template <typename T> int signum(T x) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> double Heaviside (T x) {
    if (x < T(0)) { return T(0); }
    else          { return x;    }
}

#endif