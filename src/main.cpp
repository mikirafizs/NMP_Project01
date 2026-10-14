#include <iostream>
#include "numerics.hpp"
#include <cblas.h>
#include <fstream>
#include <iomanip>
#include <chrono>

int main(){
    const size_t n = 4;
    const size_t ld = n;

    //30 times matrix vector multiplication for gamma = 2.0

    std::vector <double> x(n,1);
    std::vector <double> a = amat(n, ld, 2.0);
    std::vector <double> y = x;
    std::ofstream norm_gamma2;
    norm_gamma2.open("../data/norm_gamma2.txt");

    if(!norm_gamma2.is_open()){
        std::cerr<<"Cannot open file.\n";
        return -1;
    }

    norm_gamma2 << "# step norm\n";
    
    
    for(int m = 0; m < 30; m++){
        y = matvec(a, y, ld, n);
        norm_gamma2 << m+1 <<" "<< std::fixed << std::setprecision(6) << vnorm(y) << "\n";
    }
    norm_gamma2.close();

    //30 times matrix vector multiplication for gamma = 3.0
    
    a = amat(n, ld, 3.0);
    y = x;

    std::ofstream norm_gamma3;
    norm_gamma3.open("../data/norm_gamma3.txt");

    if(!norm_gamma3.is_open()){
        std::cerr<<"Cannot open file.\n";
        return -1;
    }
    norm_gamma3 << "# step norm\n";
    
    for(int m = 0; m < 30; m++){
        y = matvec(a, y, ld, n);
        norm_gamma3 << m+1 <<" "<< std::fixed << std::setprecision(6) << vnorm(y) << "\n";
    }
    
    norm_gamma3.close();

    // C = A*A Matrix-Matrix product(gamma = 3)

    std::vector <double> Cvec = matrix_multiply(a, a, ld, n);
    for(size_t i = 0; i < n; i++){
        std::cout << "[ ";
        for(size_t j = 0; j < ld; j++){
            std::cout << Cvec[i*ld+j]<<" ";
        }
        std::cout << "]\n";
    }
    std::cout<<"---------------------------------------------------------------------\n";
    // C = A*A Matrix-Matrix product(cblas_dgemm)(gamma = 3)
    std::vector <double> Cvec_blas(n*ld, 0.0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, static_cast <int> (n), static_cast <int> (n), static_cast <int> (n), 1.0, &a[0], static_cast <int>(ld), &a[0], static_cast <int>(ld), 0.0, &Cvec_blas[0], static_cast <int>(ld));
    for(size_t i = 0; i < n; i++){
        std::cout << "[ ";
        for(size_t j = 0; j < ld; j++){
            std::cout << Cvec_blas[i*ld+j]<<" ";
        }
        std::cout << "]\n";
    }
    return 0;
}