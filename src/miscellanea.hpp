#ifndef MISCELLANEA_HPP
#define MISCELLANEA_HPP

#include<iostream>
#include<fstream>
#include<complex>
#include<cmath>
#include<string>
#include<cstdlib>

#include<eigen\Eigen\Eigenvalues>

typedef unsigned int uInt;

/**
*
* miscellanea.hpp provides templates function to work with Matrix
*@tparam T Type of entries, can be float or double
*
**/


template<typename T,uInt N>
void swap(Eigen::Matrix<T,N,1>& t_eigenval,
          Eigen::Matrix<T,N,N>& t_eigenvec,
          uInt i,uInt j);



template<typename T, uInt N>
void sort(Eigen::Matrix<T,N,1>& t_eigenval,
          Eigen::Matrix<T,N,N>& t_eigenvec);



template<typename T, uInt N>
void set_zero(Eigen::Matrix<T,N,1>& t_eigenval, uInt n);



template<typename T, uInt N>
void eigen_sqrt (Eigen::Matrix<T,N,1> & t_eigenval);


template<typename T, uInt R, uInt C>
Eigen::Matrix<T,R,C> cast_real (const Eigen::Matrix<std::complex<T>, R,C>
                                                            & t_M);


template<typename T, uInt R, uInt C>
T EuclideanNorm (const Eigen::Matrix<T,R,C>& t_M);


template<typename T, uInt R, uInt C>
void write (const Eigen::Matrix<T,R,C>& t_M,std::string filename);



//
//Starting implementation
//



template<typename T,uInt N>
void swap(Eigen::Matrix<T,N,1>& t_eigenval,
          Eigen::Matrix<T,N,N>& t_eigenvec,
          uInt i, uInt j){

  /**
  *Swap the ith and the jth eigenvalues and the
  *corresponding eigenvectors
  *
  *@param t_eigenval the vector that contains the eigenvalues to swap
  *@param t_eigenvec eigenvectors corresponding to the eingenvalues
  *@param i position of the first value to swap
  *@param j position of the first value to swap
  *
  *@note this function modiphy t_eigenval, t_eigenvec themselfs
  */

    T temp_eigenval = t_eigenval[i];
    Eigen::Matrix<T,N,1> temp_eigenvec = t_eigenvec.col(i);

    t_eigenval[i] = t_eigenval[j];
    t_eigenvec.col(i) = t_eigenvec.col(j);

    t_eigenval[j] = temp_eigenval;
    t_eigenvec.col(j) = temp_eigenvec;
  }



template<typename T, uInt N>
void sort(Eigen::Matrix<T,N,1>& t_eigenval,
          Eigen::Matrix<T,N,N>& t_eigenvec){

  /**
  *Sort the eigenvalues absolute value in decreasin order
  *
  *@param t_eigenval eingevalues to sort
  *@param t_eigenvec eigenvectors corresponding to the eigenvalues
  *
  *@note modiphy eigenvectors and eigenvalues themselfs
  **/
  unsigned int nSwap;
  do{
    nSwap = 0;
    for (int i=1; i< N ; ++i){

      if(std::fabs(t_eigenval[i]) > std::fabs(t_eigenval[i-1])){
        swap<T,N>(t_eigenval, t_eigenvec, i-1,i);
        nSwap ++;
      }

      else continue;
    }
  }while(nSwap > 0);
}



template<typename T, uInt N>
void set_zero(Eigen::Matrix<T,N,1>& t_eigenval, unsigned int n){

  /**
  * Set to zero all the N-n lowest values
  *
  *@tparam N, lenght of the vector(number of eigenvalues)
  *@param t_eigenval: vector that contains the eigenvalues
  *@param n : number of eigenvalues to not set to 0
  *
  *@note Moduphy t_eigenvec itself
  */

  for(unsigned int i=n; i< N ; ++i) t_eigenval[i]=0;
}



template<typename T, uInt N>
void eigen_sqrt (Eigen::Matrix<T,N,1> & t_eigenval){

  /**
  * Take the square root of t_eigenval
  *
  *@tparam N lenght of the vector(number of eigenvalues)
  *@param t_eigenval vector that contains the eigenvalues
  *
  *@note modiphy t_eigenval itself
  */
  for(unsigned int i=0; i<N; ++i)
    t_eigenval[i] = std::sqrt(t_eigenval[i]);
}



template<typename T, uInt R, uInt C>
Eigen::Matrix<T,R,C> cast_real (const Eigen::Matrix<std::complex<T>, R,C>&
                                                            t_M){
  /**
  * Converts a complex matrix to a real one
  *
  *@param t_M, complex matrix to cast
  *
  *@return Real matrix
  */
  Eigen::Matrix<T,R,C> Real;
  for(unsigned int i=0; i<R; ++i){
    for(unsigned int j=0; j<C;++j){
      Real(i,j) = t_M(i,j).real();
    }
  }
  return Real;
  }



  template<typename T, uInt R, uInt C>
  T EuclideanNorm (const Eigen::Matrix<T,R,C>& t_M){
    /**
    *Compute the eucliden norm of the matrix
    *
    *@tparam R number of rows of the input matrix
    *@tparam C  number of columns of the input matrix
    *
    *@param t_M matrix from which compute the euclidean norm
    *
    *@returns euclidean norm of the matrix
    **/

    return std::sqrt((t_M*t_M.transpose()).trace());
  }


  template<typename T, uInt R, uInt C>
  void write (const Eigen::Matrix<T,R,C>& t_M,std::string filename){


    std::ofstream out;
    out.open(filename);
    if(!out){
      std::cerr<<filename<<" could not be opened"<<std::endl;

    }

    for(unsigned int i=0; i<R; i++)
      for(unsigned int j=0; j<C; j++){
        out<<t_M(i,j)<<'\t';
      }
      out<<std::endl;
  }

#endif
