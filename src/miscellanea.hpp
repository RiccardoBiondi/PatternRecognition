#ifndef MISCELLANEA_HPP
#define MISCELLANEA_HPP

#include<complex>
#include<Eingen\Eigenvalues>




template<typename T, unsigned int N, unsigned int d>
void swap(Eigen::Matrix<T, N,1>& t_eigenval,Eigen::Matrix<T,N,N>& t_eigenvect
            unsigned int i, unsigned int j);

template<typename T, unsigned int N ,unsigned int d>
void sort(Eigen::Matrix<T,N,1>& t_eigenval, Eigen::Matrix<T,N,N>& t_eigenvec);


template<typename T, unsigned int N, unsigned int d>
void set_zero(Matrix<T,N,1> t_eigenval);
