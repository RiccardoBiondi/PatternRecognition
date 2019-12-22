#ifndef EDMATRIX_HPP
#define EDMATRIX_HPP

#include <Eigen\Eigenvalues>
/**
 *  EDMatrix v2.0.0
 *  Starting Generation: 2019-12-21
 *  ----------------------------------------------------------
 *  This file contains the implementations of the template class EDMatrix and some external methods usefull to work with Euclidena Distance Matrix
 *
 * Requirments:
 *   Eigen 3.3.7
 *I've performed test by using Check v2.11.0
 *
 */


template<typename T, unsigned int N, unsigned int d>
class EDMatrix{

/**
* EDMatrix class, provides methods to build and works with Euclidean Distance Matrix(EDM)
* EDMatrix need Eigen 3.3.7 to work.
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* @tparam T type of entries of the EDM, can be double, float or int
* @tparam unsigned int N, number of points from which EDM is created. Is the dimension of the NxN EDM.
* @tparam unsigned ind dimension of the euclidean space in which points live
*/

  private:
    Eigen::Matrix<T,N,N> m_EDM;
    Eigen::Matrix<int,N,N> m_Mask;

  public:
    //Constructors
    EDMatrix(){};
    EDMatrix(Eigen::Matrix<T,N,N> t_EDM,
            Eigen::Matrix<int,N,N> t_Mask = Eigen::Matrix<int,N,N>::Ones()):
              m_EDM(t_EDM),m_Mask(t_Mask){};

    EDMatrix(const EDMatrix<T,N,d>& o_EDMatrix):m_EDM(o_EDMatrix.getEDM()),
                                                m_Mask(o_EDMatrix.getMask()){};
    //operator
    //EDMatrix<T,N,d> operator& = (const EDMatrix<T,N,d>& t_EDMatrix);


    //getter and setter
    Eigen::Matrix<T,N,N> getEDM() const;
    Eigen::Matrix<int,N,N> getMask() const;

    void setEDM(const Eigen::Matrix<T,N,N>& t_EDM);
    void setMask(const Eigen::Matrix<int,N,N>& t_Mask);

    //methods
    bool is_hollow();
    bool is_positive();
    bool is_symmetric();
    bool is_triang_inh();
    bool is_EDM();

    void make_hollow();
    void make_positive();
    void trim();

    Eigen::Matrix<T,N,N> gc_matrix();
    Eigen::Matrix<T,N,N> gramm();

    T frobenius_norm();





};





//Declaratin of externam methods
template<typename T,unsigned int N, unsigned int d>
EDMatrix<T,N,d> EVTreshold(EDMatrix<T,N,d> t_EDM);

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,N> ClassicalMDS(EDMatrix<T,N,d> t_EDM);






//Starting definitions

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

        check += sqrt(m_EDM(i,j)) <= sqrt(m_EDM(i,k))+sqrt(m_EDM(k,j))? 0:1;
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










#endif
