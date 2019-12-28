#ifndef MISCELLANEA_HPP
#define MISCELLANEA_HPP

#include<complex>
#include<cmath>
#include<Eigen\Eigenvalues>

//
//Staring declarations
//


template<typename T,unsigned int N>
void swap(Eigen::Matrix<T,N,1>& t_eigenval,Eigen::Matrix<T,N,N>& t_eigenvec,
            unsigned int i, unsigned int j);

template<typename T, unsigned int N>
void sort(Eigen::Matrix<T,N,1>& t_eigenval, Eigen::Matrix<T,N,N>& t_eigenvec);


template<typename T, unsigned int N>
void set_zero(Eigen::Matrix<T,N,1>& t_eigenval, unsigned int n);


template<typename T, unsigned int N>
void eigen_sqrt (Eigen::Matrix<T,N,1> & t_eigenval);


template<typename T, unsigned int R, unsigned int C>
Eigen::Matrix<T,R,C> cast_real (const Eigen::Matrix<std::complex<T>, R,C>&
                                                            t_eigenvec);








//
//Starting implementation
//


template<typename T,unsigned int N>
void swap(Eigen::Matrix<T,N,1>& t_eigenval,Eigen::Matrix<T,N,N>& t_eigenvec,
            unsigned int i, unsigned int j){

  /**
  * This function swap the ith and the jth eigenvalues and the corresponding eigenvectors
  *
  *@params t_eigenval, the vector that contains the eigenvalues to swap
  *@params t_eigenvec, tha matrix which contians the eigenvectors corresponding to the eingenvalues
  @params i,j , the position of the value t swap
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






template<typename T, unsigned int N>
void sort(Eigen::Matrix<T,N,1>& t_eigenval, Eigen::Matrix<T,N,N>& t_eigenvec){

  /**

  */
  unsigned int nSwap;
  do{
    nSwap = 0;
    for (int i=1; i< N ; ++i){

      if(t_eigenval[i] > t_eigenval[i-1]){
        swap<T,N>(t_eigenval, t_eigenvec, i-1,i);
        nSwap ++;
      }

      else continue;
    }
  }while(nSwap > 0);
}



template<typename T, unsigned int N>
void set_zero(Eigen::Matrix<T,N,1>& t_eigenval, unsigned int n){

  /**
  * Set to zero all the N-n lowest values
  *
  *@tparam T type of values of the vector. It can be double, float or int.
  *@tparam N lenght of the vector(number of eigenvalues)
  *@param t_eigenval: vector that contains the eigenvalues
  *@param n : number of eigenvalues to not set to 0
  *
  *@note Moduphy t_eigenvec itself
  */

  for(unsigned int i=n; i< N ; ++i) t_eigenval[i]=0;
}


template<typename T, unsigned int N>
void eigen_sqrt (Eigen::Matrix<T,N,1> & t_eigenval){

  /**
  * Take the square root of t_eigenval
  *
  *@tparam T type of values of the vector. It can be double, float or int.
  *@tparam N lenght of the vector(number of eigenvalues)
  *@param t_eigenval: vector that contains the eigenvalues
  *
  *@note modiphy t_eigenval itself
  */
  for(unsigned int i=0; i<N; ++i)
    t_eigenval[i] = std::sqrt(t_eigenval[i]);
}


template<typename T, unsigned int R, unsigned int C>
Eigen::Matrix<T,R,C> cast_real (const Eigen::Matrix<std::complex<T>, R,C>&
                                                            t_eigenvec){
  /**
  *
  */
  Eigen::Matrix<T,R,C> Real;
  for(unsigned int i=0; i<R; ++i){
    for(unsigned int j=0; j<C;++j){
      Real(i,j) = t_eigenvec(i,j).real();
    }
  }
  return Real;
  }






#endif
