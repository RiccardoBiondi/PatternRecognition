#include <iostream>
#include <cmath>
#include <eigen\Eigen\Eigenvalues>
#include <eigen\Eigen\SVD>

#include "miscellanea.hpp"

/**
 *  EDMatrix v2.0.0
 *  Starting Generation: 2020-03-20
 *  ----------------------------------------------------------
 *  This file contains the implementations of the template class Prcrustes
 *that contains tools to perform the Procrustes alignment of two set of *point in the euclidean space
 *
 * Requirments:
 *   Eigen 3.3.7, in particular the module Eigenvalues and the module SVD
 *I've performed test by using Catch v2.11.0
 *
 */

template<typename T, unsigned int N, unsigned int d>
class Procrustes{



  /**
  * \class Procrustes, provides methods to perform the prcrustes analysis
  *
  * Procrustes use Eigen 3.3.7 to work.
  * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  * @tparam T type of entries of the PointSet, can be double or float.
  * @tparam unsigned int N, number of points in the two set of points
  * @tparam unsigned int d, dimension of the euclidean space in
  *
  *
  */

private:

  //Reference
  Eigen::Matrix<T,d,N> m_Reference ;
  // Reference set of points
  Eigen::Matrix<T,d,N> m_Configuration;
  // Set of point to align to m_Reference
  Eigen::Matrix<T,d,N> m_NewConfiguration;
  //Set of poits aligned  with reference


  T m_Scale; //Scaling factor
  Eigen::Matrix<T,d,d> m_Rotation; //Rotation Matrix
public:



  Procrustes (const Eigen::Matrix<T,d,N>& t_Reference,
              const Eigen::Matrix<T,d,N>& t_Configuration);

  Procrustes (const Procrustes<T,N,d>& o_procrustes);

  Eigen::Matrix<T,d,N> getReference() const;
  Eigen::Matrix<T,d,N> getConfiguration() const;
  Eigen::Matrix<T,d,N> getNewConfiguration() const;
  Eigen::Matrix<T,d,d> Rotation() const;
  T Scale() const;


  void Centering();
  void Normalize();
  void FullOPA ();

};


//
//STARTING IMPLEMENTATION
//



//Parametric constructor
template<typename T,unsigned int N,unsigned int d>
Procrustes<T,N,d>:: Procrustes (const Eigen::Matrix<T,d,N>& t_Reference,
                                const Eigen::Matrix<T,d,N>& t_Configuration):
            m_Reference(t_Reference),
            m_Configuration(t_Configuration),
            m_NewConfiguration(t_Configuration),
            m_Rotation(Eigen::Matrix<T,d,d>:: Identity()),
            m_Scale(1.){
              /**
              Parametric Constructor
              @params t_Reference  dxN matrix which contains the reference
              *set of points
              @params t_COnfiguration dxN matrix which contains the set of
              *points to aligna to t_Reference
              **/
            }


//Copy Constructor
template<typename T, unsigned int N, unsigned int d>
Procrustes<T,N,d>:: Procrustes (const Procrustes<T,N,d>& o_procrustes):
                    m_Reference(o_procrustes.getReference()),
                    m_Configuration(o_procrustes.getConfiguration()),
                    m_Rotation(o_procrustes.Rotation()),
                    m_Scale(o_procrustes.Scale()){
                      /**
                      * Copy constructor
                      **/
                    }


//Getter
template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,d,N> Procrustes<T,N,d>:: getReference() const{
  /**
  @returns The reference set of points
  **/
  return m_Reference;
}

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,d,N> Procrustes<T,N,d>:: getConfiguration() const{
  /**
  *@returns the configuration set of points
  **/
  return m_Configuration;

}


template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,d,N> Procrustes<T,N,d>:: getNewConfiguration() const{
  /**
  *If procustes analysis has been performed
  *@returns the aligned set of points
  * If procrustes analysis hasn't been performed yet
  *@returns the original configuration
  **/
  return m_NewConfiguration;
}

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,d,d> Procrustes<T,N,d> :: Rotation() const {
  /**
  * @returns the rotation matrix used to align m_Configuration to m_Reference
  **/

  return m_Rotation;
}


template<typename T, unsigned int N, unsigned int d>
T Procrustes<T,N,d>:: Scale () const{
  /**
  *@returns the scaling factor
  **/

  return m_Scale;
}






template<typename T, unsigned int N, unsigned int d>
void Procrustes<T,N,d>:: Centering(){
  /**
  *Center the two set of points by using the centering matrix
  * C = I - 1/N11
  *@note Modiphy m_Reference and m_Configuration themself
  **/
  Eigen::Matrix<T,N,N> C = (Eigen::Matrix<T,N,N>::Identity()
                            -(1./(double)N)*Eigen::Matrix<T,N,N>::Ones());

  m_Reference = m_Reference*C;
  m_Configuration = m_Configuration*C;

}

template<typename T, unsigned int N, unsigned int d>
void Procrustes<T,N,d>:: Normalize(){
  /**
  *Normilize the two set of points by using the Euclidean Norm
  *@note Modiphy m_Configuration and m_Reference themself
  **/
  T N_1 = EuclideanNorm<T,d,N>(m_Configuration);
  T N_2 = EuclideanNorm<T,d,N>(m_Reference);
  m_Configuration = m_Configuration/N_1;
  m_Reference = m_Reference/N_2;
}



template<typename T, unsigned int N, unsigned int d>
void Procrustes<T,N,d>:: FullOPA(){
  /**
  *Perform the Full Ordinary Procrustes Analysis
  **/

  Centering();
  Normalize();

  Eigen::Matrix<T,N,d> X_1 = getReference().transpose();
  Eigen::Matrix<T,N,d> X_2 = getConfiguration().transpose();

  Eigen::Matrix<T,-1,-1> W;
  W = (X_1.transpose()*X_2);

  Eigen::BDCSVD<Eigen::Matrix<T,-1, -1>> S(W,22);
  m_Rotation = S.matrixU()*(S.matrixV().transpose());

  m_NewConfiguration = (X_2*m_Rotation.transpose()).transpose();

}
