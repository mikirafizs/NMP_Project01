#ifndef NUMERICS_HPP
#define NUMERICS_HPP
#include <cmath>
#include <cstddef>
#include <vector>
#include <stdexcept>

inline double absdiff(std::size_t a, std::size_t b) {
    double i = static_cast <double>(a);
    double j = static_cast <double>(b);
    return abs(i-j);
}

inline std::vector <double> matrix_multiply(const std::vector <double>& a, const std::vector <double>& b, std::size_t ld, std::size_t n){
    if(ld<n){
        throw std::invalid_argument("Leading dimension is smaller than size of the matrix.");
    }
    if(a.size() != n*ld || b.size() != n*ld){
        throw std::invalid_argument("Cannot perform matrix-matrix multiplication (dimensions differ).");
    }
    std::vector <double> c (n*ld, 0.0);
    for(std::size_t i = 0; i < n; i++){
        for(std::size_t k = 0; k < n; k++){
            double s = 0.0;
            for(std::size_t j = 0; j < n; j++){
                s += a[i*ld+j]*b[j*ld+k];
            }
        c[i*ld+k] = s;
        }
    }
    return c;
    
}

inline std::vector <double> amat(std::size_t n, std::size_t ld, double gamma){
    if(ld<n){
        throw std::invalid_argument("Leading dimension is smaller than size of the matrix.");
    }
    std::vector <double> a(n*ld, 0.0);
    for(std::size_t i = 0; i<n; ++i) {
        for (std::size_t j=0; j<n; ++j) {
            double diff = absdiff(i,j);
            a[i * ld + j] = (1/(diff*diff*diff*diff+gamma));
        }
    }
    return a; 
}

inline std::vector<double> matvec(const std::vector <double>& a, const std::vector <double>& x, std::size_t ld, std::size_t n){
    if(ld<n){
        throw std::invalid_argument("Leading dimension is smaller than size of the matrix.");
    }
    if(a.size() != n*ld || x.size() != n){
        throw std::invalid_argument("Cannot perform matrix-vector multiplication (dimensions differ).");
    }
    std::vector <double> y(n, 0.0);
    for(std::size_t i = 0; i < n; i++){
        double sum = 0.0;
        for(std::size_t j = 0; j<n; j++){
            sum += a[i * ld + j]*x[j];
        }
        y[i] = sum;
    }
    return y;
}

inline double vnorm(const std::vector <double>& y){
    double norm = 0.0;
    for(std::size_t i = 0; i < y.size(); i++){
        norm += y[i]*y[i];
    }
    return std::sqrt(norm);
}
#endif