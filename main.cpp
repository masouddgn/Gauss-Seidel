#include <iostream>
#include <omp.h>
#include<fstream>
#include "Data.h"
#include "Timer.h"
#define pi  3.1416
#include "cmath"

int Nx , Ny;
int iteration_number;
int processor_number;
double hx , hy , alpha , Beta , gama;



using namespace std;

void set_b(Data& data){

    for (int i = 0 ; i < Nx ; i++){
        for (int j = 0 ; j < Ny ; j++){
            data(i,j) = 2*pi*pi*cos(pi*hx*i)*cos(pi*hy*j);
        }
    }

}

void set_BC_for_U(Data& data) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (i == 0) {
                data(0, j) = cos(pi * hy * j);
            } else if (i == Nx - 1) {
                data(Nx - 1, j) = -cos(pi * hy * j);
            } else if (j == 0) {
                data(i, 0) = cos(pi * hx * i);
            } else if (j == Ny - 1) {
                data(i, Ny - 1) = -cos(pi * hx * i);
            } else {
                data(i, j) = 0.0;
            }
        }
    }
}


void gauss_seidel_red_black(Data& u, Data& b, int itr) {
    //double c = -0.25 * h * h;

    for (int it = 0; it < itr; it++) {
        // Update Red cells
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < Nx-1; i++) {
            for (int j = 1; j < Ny-1; j++) {
                if ((i + j) % 2 == 0) { // Check for Red (even sum indices)
                    u(i,j) = (1.0/alpha) * (b(i,j) - Beta*(u(i+1 , j) + u(i-1 , j)) - gama* (u(i , j-1) + u(i , j+1)));
                }
            }
        }

        // Update Black cells
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < Nx-1; i++) {
            for (int j = 1; j < Ny-1; j++) {
                if ((i + j) % 2 == 1) { // Check for Black (odd sum indices)
                    u(i,j) = (1.0/alpha) * (b(i,j) - Beta*(u(i+1 , j) + u(i-1 , j)) - gama* (u(i , j-1) + u(i , j+1)));
                }
            }
        }
    }
}


void compute_r(Data& u , Data& b , Data& r){

    //double c = (1/(h*h));

    for (int i =0 ; i < Nx ; i++){
        for (int j = 0 ; j < Ny ; j++){

            if (i==0){
                r(i,j) = 0.0;
            }else if (i==Nx-1){
                r(i,j) = 0.0;
            } else if (j==0){
                r(i,j) = 0.0;
            }else if (j == Ny-1){
                r(i,j) = 0.0;
            }else {
                r(i,j) = b(i,j) - Beta*(u(i+1 , j) + u(i-1 , j)) - gama*(u(i , j-1) + u(i , j+1)) - alpha*u(i,j) ;
    
            }


        }
    }

}


double norm2(Data& data){

    double error = 0.0;

    for(int i=0 ; i < Nx-1 ; i++){
        for (int j = 0 ; j < Ny-1 ; j++){

            error += data(i,j)*data(i,j);

        }
    }

    error = sqrt( error / (Nx*Ny));
    return error ;

}

void printSolution(Data grid){

    ofstream solutionFile;
    solutionFile.open("solution.txt");
    solutionFile << "# x y grid(x,y) \n";

    for ( int i = 0 ; i < Nx ; i++)
    {
        for (int j = 0 ; j < Ny ; j++){
            solutionFile << (i*hx) <<" " << (j*hy) << " " << grid(i,j) << "\n" ;
        }
        solutionFile << "\n" ;
    }
    solutionFile.close();
    

}

int main(int argc, char* argv[]){

    siwir::Timer test_time;
    test_time.reset();

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <Nx> <Ny> <iteration_number>\n";
        return 1;
    }

    Nx = std::stoi(argv[1]); // Convert the first argument to an integer
    Ny = std::stoi(argv[2]); // Convert the first argument to an integer
    iteration_number = std::stoi(argv[3]); // Convert the second argument to an integer

    processor_number = 1;

    omp_set_num_threads(processor_number);    

    
    
    
    hx = 1.0 / (Nx-1) ;
    hy = 1.0 / (Ny-1) ;

    printf("hx is %f \n" , hx);
    printf("hx is %f \n" , hy);

    alpha = ((2.0/(hx*hx)) + (2.0/(hy*hy)));
    Beta = -(1.0/(hx*hx));
    gama = -(1.0/(hy*hy));

    
    double error;
    Data b(Nx , Ny) , u(Nx , Ny) , r(Nx , Ny);
    set_b(b);
    set_BC_for_U(u);
    compute_r( u , b , r );
    error = norm2(r) ;
    printf("residual value at before Gauss-Seidel is : %f \n" , error) ;

    

    double sequential_part = test_time.elapsed();
    test_time.reset();
    gauss_seidel_red_black(u , b , iteration_number);
    double parallel_time = test_time.elapsed();
    test_time.reset();
    compute_r( u , b , r );
    error = norm2(r);
    printf("residual value at end is : %f \n" , error) ;

    printf("the sequential_part is %f \n" , sequential_part);
    printf("the parallel_time is %f \n" , parallel_time);
    printf("the Wall clock time is %f \n" , parallel_time+sequential_part);

    printSolution(u) ;

    return 0;
}
