#ifndef POINT_SET_HPP
#define POINT_SET_HPP

#include <random>
#include "EDMatrix.hpp"

template<typename T, unsigned int N, unsigned int d>
class PointSet{

  private:

    Eigen::Matrix<T, d, N> m_point_set;
    EDMatrix<T,N,d> m_EDM;

    //seed
    std::random_device rd;
    //set seed

    //private functions
    //compute the euclidean distance
    T EuclideanDistance(const Eigen::Matrix<T,d,1>& t_PointA,
                        const Eigen::Matrix<T,d,1>& t_PointB);

  public:

    //Constructors
    PointSet ();
    PointSet(T h);
    PointSet(const PointSet& o_Point_Set);

    //operator

    //getter and setter
    Eigen::Matrix<T,d,N> getPointSet() const;
    EDMatrix<T,N,d> getEDM() const;
    //Adder
    void AddNoise(T sdev);
    void AddMissEntires();



};

//
//Private Fucntions
//

//Euclidean EuclideanDistance

template<typename T, unsigned int N, unsigned int d>
T PointSet<T,N,d>::EuclideanDistance(const Eigen::Matrix<T,d,1>& t_PA,
                                     const Eigen::Matrix<T,d,1>& t_PB){
  /**
  *Private function.
  *@returns the distance square between the poin PA and PB
  **/
  T distance = 0.;
  for(unsigned int i=0; i<d; i++){
    distance += std::pow((t_PA(i) - t_PB(i)),2);
  }
  return distance;

                                     }

//Constructors

//Default Constructor
template<typename T, unsigned int N, unsigned int d>
PointSet<T,N,d>:: PointSet(){
  /**
  *Default constructor
  **/
}

//Parametric constructor
template<typename T, unsigned int N, unsigned int d>
PointSet<T,N,d>:: PointSet(T h) {
  /**
  *Parametric constructor
  **/

  //prepare the genertor
  std::mt19937 gen(rd());
  std::uniform_real_distribution<T> dis(-h/2, h/2);
  //genrate the point
  for(unsigned int i=0; i< d; i++){
    for(unsigned int j=0; j<N; j++){
    m_point_set(i,j) = dis(gen);
    }
  }

  //create the EDM
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

//Copy Constructor

template<typename T, unsigned int N, unsigned int d>
PointSet<T,N,d> ::PointSet(const PointSet& o_PointSet):
  m_point_set(o_PointSet.getPointSet()),m_EDM(o_PointSet.getEDM()){
    /**
    *Copy Constructor
    **/
  }
//Getter



template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,d,N> PointSet<T,N,d>:: getPointSet() const {

  /**
  **/
  return m_point_set;
}





template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> PointSet<T,N,d>:: getEDM() const {

  /**
  **/
  return m_EDM;

}


//Adding methods

template<typename T, unsigned int N, unsigned int d>
void PointSet<T,N,d>:: AddMissEntires(){
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



template<typename T, unsigned int N, unsigned int d>
void PointSet<T,N,d>:: AddNoise(T sdev){
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
