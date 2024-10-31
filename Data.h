#ifndef DATA_H
#define DATA_H
#include <memory>

struct Data {

    int nx , ny ;
    std::shared_ptr<double[]> my_array; //This is one of the tools in the STD Library.
    // Think of it as a smart box that can hold something, like numbers, and it keeps track of how many people are using it.
    // When no one needs it anymore, it cleans up after itself.

    Data() = default; // it Helps to make an instance array from this class.
    Data(int Nx , int Ny) : nx(Nx) , ny(Ny), my_array(new double[nx*ny]){}; // n and my_array are the features of this class
    inline double& operator () (int i , int j) {return my_array[i*nx + j];} // to look at the element in position (i,j) and modify this
    //inline here says to computer that operator () with two input arguments(i,j) in this class means : my_array[i*n + j]
    //The double& before the operator() function indicates the return type of the function
    //The & after double indicates that the function returns a reference to a double.
    inline double operator () (int i , int j) const {return my_array[i*nx + j];} // by using " new " you're telling the computer to allocate memory for a
    // new array of doubles with a size of "n*n"

};
#endif // DAHA_H
