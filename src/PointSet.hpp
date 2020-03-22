#ifndef POINT_SET_HPP
#define POINT_SET_HPP

#include <random>
#include <string>
#include "EDMatrix.hpp"

typedef unsigned int uInt;

/**
*
*This file contains the definition a the template class PointSet, that allow
*to  generates set of points in an euclidean space
**/
template<typename T, uInt N, uInt d>
class PointSet{

  /**
  *\class PointSet  allow the generation of set of Points in the *euclidean space
  *
  *@tparam N is the scalar type i.e. the type of coefficient. It can be
  *float or double.
  *@tparam N number of points of the PointSet
  *@tparam d dimension of the euclidean space
  **/

  private:

    Eigen::Matrix<N, d, N> m_point_set;
    EDMatrix<Scalar,N,d> m_EDM;

    std::random_device rd;//seed

    T EuclideanDistance(const Eigen::Matrix<T,d,1>& t_PointA,
                        const Eigen::Matrix<T,d,1>& t_PointB);

  public:

    PointSet ();
    PointSet(T h);
    PointSet(const PointSet& o_Point_Set);

    Eigen::Matrix<T,d,N> getPointSet() const;
    EDMatrix<T,N,d> getEDM() const;

    void AddNoise(T sdev);
    void AddMissEntires();
};

//
//Private Fucntions
//

template<typename T, uInt N, uInt d>
T PointSet<T,N,d>::EuclideanDistance(const Eigen::Matrix<T,d,1>& t_PA,
                                     const Eigen::Matrix<T,d,1>& t_PB){
  /**
  *Private function.
  *@param t_PA first point
  *@param t_PB second point
  *@returns the distance square between the poin PA and PB
  **/
  T distance = 0.;
  for(unsigned int i=0; i<d; i++){
    distance += std::pow((t_PA(i) - t_PB(i)),2);
  }
  return distance;
}



template<typename T, uInt N, uInt d>
PointSet<T,N,d>:: PointSet(){

  /**
  *Default constructor
  **/
}



template<typename T, uInt N, uInt d>
PointSet<T,N,d>:: PointSet(T h) {

  /**
  *Parametric constructor
  *
  *Generates a set of points uniformily distributed in and hypercube of
  *side h and compute the correspondig Euclidean Distance Matrix
  *
  *@param h side of the hypercube
  *
  **/

  std::mt19937 gen(rd());
  std::uniform_real_distribution<T> dis(-h/2, h/2);

  for(unsigned int i=0; i< d; i++){
    for(unsigned int j=0; j<N; j++){
    m_point_set(i,j) = dis(gen);
    }
  }

  Eigen::Matrix<T,N,N> D;
  for(unsigned int i=0; i <N; i++){
    for(unsigned int j=i; j<N; j++){
      D(i,j) = EuclideanDistance(getPointSet().col(i),
                                 getPointSet().col(j));
      D(j,i) = D(i,j);

    }
  }

  m_EDM.setEDM(D);
}



template<typename T, uInt N, uInt d>
PointSet<T,N,d> ::PointSet(const PointSet& o_PointSet):
  m_point_set(o_PointSet.getPointSet()),m_EDM(o_PointSet.getEDM()){
    /**
    *Copy Constructor
    **/
  }



template<typename T, uInt N, uInt d>
Eigen::Matrix<T,d,N> PointSet<T,N,d>:: getPointSet() const {

  /**
  *@returns m_point_set
  *
  *@see getEDM
  **/
  return m_point_set;
}





template<typename T, uInt N, uInt d>
EDMatrix<T,N,d> PointSet<T,N,d>:: getEDM() const {

  /**
  *@returns m_EDM
  *
  *@see getPointSet
  **/
  return m_EDM;
}



template<typename T, uInt N, uInt d>
void PointSet<T,N,d>:: AddMissEntires(){

  /**
  *Add a Mask to the EDM.
  *
  *@see AddNoise
  **/
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(0,1);

  Eigen::Matrix<int,N,N> M;
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=i; j<N;j++){
      M(i,j) = dis(gen);
      M(j,i) = M(i,j);
    }
  }
  m_EDM.setMask(M);
}



template<typename T, uInt N, uInt d>
void PointSet<T,N,d>:: AddNoise(T sdev){

  /**
  *Add a Gaussian Noise matrix to the EDM
  *
  *@param sdev standard deviation of Noise
  *
  *@see AddMissEntires
  **/
  std::mt19937 gen(rd());
  std::normal_distribution<T> dis(0,sdev);
  Eigen::Matrix<T,N,N> R;
  for(int i=0; i< N; i++){
    for(int j=i; j<N; j++){
      R(i,j) = dis(gen);
      R(j,i) = R(i,j);
    }
  }
  m_EDM.setNoise(R);
}



#endif
