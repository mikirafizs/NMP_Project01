#include <iostream>
#include "numerics.hpp"
#include <cblas.h>
#include <fstream>
#include <iomanip>
#include <time.h>

int main(){
    const size_t n = 2000;
    const size_t ld = n;
    time_t t_start, t_end, t_statr_blas, t_end_blas;

    //30 times matrix vector multiplication for gamma = 2.0

    std::vector <double> x(n,1);
    std::vector <double> a = amat(n, ld, 2.0);
    std::vector <double> y = x;
    std::ofstream norm_gamma2;
    norm_gamma2.open("../data/norms_gamma2.csv");

    if(!norm_gamma2.is_open()){
        std::cerr<<"Cannot open file.\n";
        return -1;
    }

    norm_gamma2 << "iteration,norm2\n";
    
    
    for(int m = 0; m < 30; m++){
        y = matvec(a, y, ld, n);
        norm_gamma2 << m+1 <<","<< std::fixed << std::setprecision(6) << vnorm(y) << "\n";
    }
    norm_gamma2.close();

    //30 times matrix vector multiplication for gamma = 3.0
    
    a = amat(n, ld, 3.0);
    y = x;

    std::ofstream norm_gamma3;
    norm_gamma3.open("../data/norms_gamma3.csv");

    if(!norm_gamma3.is_open()){
        std::cerr<<"Cannot open file.\n";
        return -1;
    }
    norm_gamma3 << "iteration,norm2\n";
    
    for(int m = 0; m < 30; m++){
        y = matvec(a, y, ld, n);
        norm_gamma3 << m+1 <<","<< std::fixed << std::setprecision(6) << vnorm(y) << "\n";
    }
    
    norm_gamma3.close();

    // C = A*A Matrix-Matrix product(gamma = 3)
    
    t_start = time(NULL);
    std::vector <double> Cvec = matrix_multiply(a, a, ld, n);
    t_end = time(NULL);

    // C = A*A Matrix-Matrix product(cblas_dgemm)(gamma = 3)
    
    t_statr_blas = time(NULL);
    std::vector <double> Cvec_blas(n*ld, 0.0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, static_cast <int> (n), static_cast <int> (n), static_cast <int> (n), 1.0, &a[0], static_cast <int>(ld), &a[0], static_cast <int>(ld), 0.0, &Cvec_blas[0], static_cast <int>(ld));
    t_end_blas = time(NULL);

    std::cout << "Time for standard matrix multiplication: " << difftime(t_end, t_start) << "\n";
    std::cout << "Time for matrix multiplication using cblas_dgemm: " << difftime(t_end_blas, t_statr_blas) << "\n";
    return 0;
}