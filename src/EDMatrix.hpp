#ifndef EDMATRIX_HPP
#define EDMATRIX_HPP

#include <Eigen\Eigenvalues>
#include <cmath>
#include <iostream>

#include "miscellanea.hpp"
/**
 *  EDMatrix v2.0.0
 *  Starting Generation: 2019-12-21
 *  ----------------------------------------------------------
 *  This file contains the implementations of the template class EDMatrix and some external methods usefull to work with Euclidean Distance Matrix
 *
 * Requirments:
 *   Eigen 3.3.7
 *I've performed test by using Check v2.11.0
 *
 */


template<typename T, unsigned int N, unsigned int d>
class EDMatrix{

/**
* \class EDMatrix, provides methods to build and works with Euclidean Distance Matrix(EDM)
* EDMatrix need Eigen 3.3.7 to work.
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* @tparam T type of entries of the EDM, can be double, float or int
* @tparam unsigned int N, number of points from which EDM is created. Is the dimension of the NxN EDM.
* @tparam unsigned ind dimension of the euclidean space in which points live
*/

  private:
    Eigen::Matrix<T,N,N> m_EDM;
    Eigen::Matrix<int,N,N> m_Mask;
    Eigen::Matrix<T,N,N> m_Noise;

  public:
    //Constructors
    EDMatrix(){};
    EDMatrix(Eigen::Matrix<T,N,N> t_EDM,
            Eigen::Matrix<int,N,N> t_Mask = Eigen::Matrix<int,N,N>::Ones()):
              m_EDM(t_EDM),m_Mask(t_Mask){};

    EDMatrix(const EDMatrix<T,N,d>& o_EDMatrix):m_EDM(o_EDMatrix.getEDM()),
                                                m_Mask(o_EDMatrix.getMask()){};
    //operator
    EDMatrix<T,N,d> operator  = (const EDMatrix<T,N,d>& t_EDMatrix);

    //sqrt
    Eigen::Matrix<T,N,N> sqrt() const;
    //getter and setter
    Eigen::Matrix<T,N,N> getEDM() const;
    Eigen::Matrix<int,N,N> getMask() const;
    Eigen::Matrix<T,N,N> getNoise() const;

    void setEDM(const Eigen::Matrix<T,N,N>& t_EDM);
    void setMask(const Eigen::Matrix<int,N,N>& t_Mask);
    void setNoise(const Eigen::Matrix<T,N,N>& t_Noise);

    //methods
    bool is_hollow();
    bool is_positive();
    bool is_symmetric();
    bool is_triang_inh();
    bool is_EDM();

    void make_hollow();
    void make_positive();
    void trim();//TODO

    Eigen::Matrix<T,N,N> gc_matrix() const;
    Eigen::Matrix<T,N,N> hadamard() const;
    Eigen::Matrix<T,N,N> gramm() const;

    T frobenius_norm();

};





//Declaratin of externam methods
template<typename T,unsigned int N, unsigned int d>
EDMatrix<T,N,d> EVTreshold(EDMatrix<T,N,d> t_EDM, unsigned int r);

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,d> ClassicalMDS(const EDMatrix<T,N,d>& t_EDM);

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,d> AlternatingDescend(const EDMatrix<T,N,d> & t_EDM);

template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> RankCompleteEDM(const EDMatrix<T,N,d>& t_EDM);


template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> OptSpace (const EDMatrix<T,N,d>& t_EDM);



//Starting definitions


//Copy assignement operator

template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> EDMatrix<T,N,d>:: operator  =
                                    (const EDMatrix<T,N,d>& t_EDMatrix){

      setEDM(t_EDMatrix.getEDM());
      setMask(t_EDMatrix.getMask());
      return *this;
                                    }

//
//sqrt operator
//
template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,N> EDMatrix<T,N,d>:: sqrt() const {

  Eigen::Matrix<T,N,N> R;
  for(unsigned int i=0; i< N ;i++){
    for(unsigned int j=i; j<N; j++){
      R(i,j) = std::sqrt(getEDM()(i,j));
      R(j,i) = std::sqrt(getEDM()(j,i));
    }
  }
  return R;

}





//
//Getter definition
//
template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,N> EDMatrix<T,N,d>:: getEDM() const{
  /**
  * getter. @returns m_EDM.
  *@see getMask
  */
  return m_EDM;
}


template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<int,N,N> EDMatrix<T,N,d>:: getMask() const{
  /**
  * getter. @returns m_Mask
  *@see getEDM
  */
  return m_Mask;
}


//
//Setter Definitions
//

template<typename T, unsigned int N, unsigned int d>
void EDMatrix<T,N,d>:: setEDM(const Eigen::Matrix<T,N,N>& t_EDM){
  /**
  * assign to m_EDM the t_EDM value
  * @param const Eigen::Matrix<T,N,N>& t_EDM
  */
  m_EDM = t_EDM;
}

template<typename T,unsigned int N, unsigned int d>
void EDMatrix<T,N,d>:: setMask(const Eigen::Matrix<int,N,N>& t_Mask){
  /**
  * set the value of m_Mask to t_Mask
  * @param const Eigen:: Matrix<int,N,N>& t_Mask
  */
  m_Mask = t_Mask;
}


//
//boolean methods
//

template<typename T, unsigned int N, unsigned int d>
bool EDMatrix<T,N,d>:: is_hollow(){
  /**
  * Check if the diagonal of the EDM is composed only by 0 values.
  *@returns true if the matrix is hollow
  *@returns false in the matrix is not hollow
  */

  int check = 0;
  for(unsigned int i=0; i<N; ++i){
    check += m_EDM(i,i) == 0 ? 0:1;
  }
  return check == 0;
}




template<typename T, unsigned int N, unsigned int d>
bool EDMatrix<T,N,d>:: is_positive(){
  /**
  * chek if all the entires of EDM are greater or equal to zero
  */
   int check = 0;
   for(unsigned int i=0; i<N;++i){
     for(unsigned int j= i; j<N; ++j )
        check += m_EDM(i,j) >= 0 ? 0:1;

   }
   return check == 0 ;
}



template<typename T, unsigned int N, unsigned int d>
bool EDMatrix<T,N,d>:: is_symmetric(){
  /**
  * Check if the matrix is symmetric
  */

  return m_EDM == m_EDM.transpose();
}



template<typename T, unsigned int N, unsigned int d>
bool EDMatrix<T,N,d>:: is_triang_inh(){
  /**
  * Check inf the triangular inhequality is verified
  */
  unsigned int check = 0;
  for(unsigned int i=0; i <N; ++i){
    for(unsigned int j=i; j< N ; ++j){
      for(unsigned int k=j; k< N ; ++k){

        check += std::sqrt(m_EDM(i,j)) <= std::sqrt(m_EDM(i,k))+std::sqrt(m_EDM(k,j))? 0:1;
      }
    }

  }
  return check == 0;
}


template<typename T,unsigned int N, unsigned int d>
bool EDMatrix<T,N,d>:: is_EDM(){
  /**
  * Check if all the proprieties of EDM are verified
  */
  return is_hollow()&&is_positive()&&is_symmetric()&&is_triang_inh();
}



//
//void methods
//


template<typename T,unsigned int N, unsigned int d>
void EDMatrix<T,N,d>:: make_hollow(){

  /**
  * Force the hollowness of the matrix by setting all the diagonal entires to zero.
  *@note Modiphy m_EDM itself
  */

  for(unsigned int i=0; i<N; ++i){
    m_EDM(i,i) = 0;
  }
}



template<typename T, unsigned int N,unsigned int d>
void EDMatrix<T,N,d>:: make_positive(){
  /**
  *Force the postivness of the matrix by setting to 0 all the negative entires.
  *@note Modiphy m_EDM itself
  */
  for(unsigned int i=0; i<N; ++i){
    for(unsigned int j=i; j<N;++j){
      m_EDM(i,j)= m_EDM(i,j)<0 ? 0: m_EDM(i,j);
      m_EDM(j,i)= m_EDM(j,i)<0 ? 0: m_EDM(j,i);

    }
  }
}





//
//Eigne::Matrix methods
//

template<typename T,unsigned int N, unsigned int d>
Eigen::Matrix<T,N,N> EDMatrix<T,N,d>:: gc_matrix() const{

  /**
  *@returns the geometric centring matrix given by:
  *@f$ J = I - \frac{1}{N}11^{T}@f$
  */
  return (Eigen::Matrix<T,N,N>::Identity()-(1./
        (double)N)*Eigen::Matrix<T,N,N>::Ones());

}




template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,N> EDMatrix<T,N,d>:: hadamard() const {

  /**
  *Perform the Hadamard product between EDM and mask
  *@return EDM with missing entires
  */

  Eigen::Matrix<T,N,N> Result;
  for(unsigned int i=0; i<N; ++i){
    for(unsigned int j=i; j<N;++j){
      Result(i,j) = m_EDM(i,j)*m_Mask(i,j);
      Result(j,i) = m_EDM(j,i)*m_Mask(j,i);
    }
  }
  return Result;
}




template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,N> EDMatrix<T,N,d>:: gramm() const{

  /**
  * @returns the Gramm matrix G, given by:
  *@f$G  = XX^{T} = -0.5JDJ@f$
  *Where:
  * X is the point set (dxN matrix)
  * J is the geometric centring matrix @see gc_matrix
  * D is the Eucliden Distance Matrix (m_EDM).
  */

  return -0.5*gc_matrix()*hadamard()*gc_matrix();



}




//
//T method
//



template<typename T, unsigned int N, unsigned int d>
T EDMatrix<T,N,d>:: frobenius_norm(){

  /**
  * Perform the frombenius norm of the EDM with missing entires
  */
  Eigen::Matrix<T,N,N> D = hadamard();
  T norm = 0;
  for(unsigned int i=0; i<N; ++i){
    for(unsigned int j=0; j<N; ++j){
      norm += std::pow(D(i,j),2);
    }
  }
  return std::sqrt(norm);
}



//
//
//
//
//NOW START THE IMPLEMENTATION OF THE EXTERNAL FUNCIONS
//
//
//
//


//
//
// INTERNAL FUNCTIONS
//


template<typename T,unsigned int N, unsigned int d>
EDMatrix<T,N,d> EVTreshold(EDMatrix<T,N,d> t_EDM, unsigned int r){

  /**
  *This function is usefull inside the Rank Complete Edm. Perform
  *and Eigenvalues decomposition and keeps only the r greater
  *eigenvalues and the corresponding eigenvalues.
  *@params t_EDM -> Euclidian Distance matrix in which we perform
  *the treshold
  *@params r -> Number of Eigenvalues to not kepp on zero
  *
  @returs An eucliden distance matrix.
  **/

  Eigen::EigenSolver<Eigen::Matrix<T,N,N>> es(t_EDM.hadamard());
  Eigen::Matrix<T,N,1> eigenval =
                          cast_real<T,N,1>(es.eigenvalues());
  Eigen::Matrix<T,N,N> eigenvec =
                          cast_real<T,N,N>(es.eigenvectors());

  sort<T,N>(eigenval, eigenvec);
  set_zero<T,N>(eigenval, r);
  Eigen::Matrix<T,N,N> D =
          eigenvec*eigenval.asDiagonal()*eigenvec.transpose();

  EDMatrix<T,N,d> R(D);
  return R;


}

//
//
//
//FUNCTIONS TO PERFORM THE MULTI DIMENSIONAL SCALING
//
//

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,d> ClassicalMDS(const EDMatrix<T,N,d>& t_EDM){

  /**
  * Perform the classical multi dimensional scaling
  *
  *@params t_EDM Euclidean distace matrix from which reconstruct the point set
  *
  é@returns Nxd matrix which is the reconstructed point set up to a translation and rotation.
  */


  Eigen::EigenSolver<Eigen::Matrix<T,N,N>> es(t_EDM.gramm());

  Eigen::Matrix<T,N,1> eigenval =
                      cast_real<T,N,1>(es.eigenvalues());
  Eigen::Matrix<T,N,N> eigenvec =
                      cast_real<T,N,N>(es.eigenvectors());

  sort<T,N>(eigenval, eigenvec);
  set_zero<T,N>(eigenval, d);
  eigen_sqrt<T,N>(eigenval);
  Eigen::Matrix<T,N,N>A =
                    eigenval.asDiagonal()*eigenvec.transpose();
  Eigen::Matrix<T,N,d> X;
  for(unsigned int i=0; i<d; ++i){
                  X.col(i) = A.transpose().col(i);}

  return X;

}
//
//
//FUNCTIONS TO PERFORM THE MATRIX COMPLETION
//
//

template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> RankCompleteEDM(const EDMatrix<T,N,d>& t_EDM, T MAX_TOL,
                                unsigned int MAX_ITER){

  /**

  **/

}



#endif
