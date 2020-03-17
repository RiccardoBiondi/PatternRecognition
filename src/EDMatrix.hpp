#ifndef EDMATRIX_HPP
#define EDMATRIX_HPP

#include <eigen\Eigen\Eigenvalues>
#include <cmath>
#include <iostream>

#include "miscellanea.hpp"

//Typedef to improve code readability
typedef unsigned int uInt;

/**
 *  EDMatrix v2.0.0
 *  Starting Generation: 2019-12-21
 *  ----------------------------------------------------------
 *  This file contains the implementations of the template class EDMatrix and some external methods usefull to work with Euclidean Distance Matrix
 *
 * Requirments:
 *   Eigen 3.3.7
 *I've performed test by using Catch v2.11.0
 *
 */


template<typename T, uInt N, uInt d>
class EDMatrix{

/**
* \class EDMatrix, provides methods to build and works with Euclidean Distance Matrix(EDM)
* EDMatrix use Eigen 3.3.7 to work.
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
* @tparam T type of entries of the EDM, can be double, float.
* @tparam unsigned int N, number of points from which EDM is
*created. Is the dimension of the NxN EDM.
* @tparam unsigned int dimension of the euclidean space in
*which points
*live
*/

  private:
    Eigen::Matrix<T,N,N> m_EDM;
    Eigen::Matrix<int,N,N> m_Mask;
    Eigen::Matrix<T,N,N> m_Noise;

  public:
    //Constructors
    EDMatrix();

    //Parametric Constructor

    EDMatrix(Eigen::Matrix<T,N,N> t_EDM,
            Eigen::Matrix<int,N,N> t_Mask =
            Eigen::Matrix<int,N,N>::Ones(),
            Eigen::Matrix<T,N,N> t_Noise = Eigen::Matrix<T,N,N>::Zero()):

            m_EDM(t_EDM),m_Mask(t_Mask),m_Noise(t_Noise){
              /**
              *Parametric Constructor
              *@params t_EDM -> NxN matrix which contains the quare of
              *euclidean distances between the point of a given pointset
              *@params t_Mask -> NxN int matrix which contains only 0
              *and 1, 0 corresponds to an miss entire, 1 to a known
              *entries. As default is set as One matrix.
              *@params t_Noise -> NxN matrix which contains the error on
              *the estimates square distances. As default is set as a
              *matrix with only 0 entires.
              *
              **/
            };

    EDMatrix(const EDMatrix<T,N,d>& o_EDM);
    //operator
    EDMatrix<T,N,d> operator = (const EDMatrix<T,N,d>& t_EDMatrix);

    //Aritmetic operators
    EDMatrix<T,N,d> operator + (const EDMatrix<T,N,d>& t_EDM);
    EDMatrix<T,N,d> operator - (const EDMatrix<T,N,d>& t_EDM);
    //Access operator
    T& operator () (uInt Row , uInt Col);

    //sqrt
    EDMatrix<T,N,d> sqrt() const;
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
    void trim();

    Eigen::Matrix<T,N,N> gc_matrix() const;
    //Eigen::Matrix<T,N,N> hadamard() const;
    Eigen::Matrix<T,N,N> gramm() const;

    T frobenius_norm();

};

//Some usefull tools to work with matrices
template<typename T,typename X, uInt Row, uInt Col>
Eigen::Matrix<T,Row,Col> hadamard(const Eigen::Matrix<T,Row,Col>& t_A,
                                  const Eigen::Matrix<X,Row,Col>& t_B);

template<typename T, uInt Row, uInt Col>
T frobenius_norm (const Eigen::Matrix<T,Row,Col>& t_M);

//Declaratin of externam methods
template<typename T, uInt N, uInt d>
EDMatrix<T,N,d> EVTreshold(EDMatrix<T,N,d> t_EDM, unsigned int r);

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,d,N> ClassicalMDS(const EDMatrix<T,N,d>& t_EDM);

template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<T,N,d> AlternatingDescend(const EDMatrix<T,N,d> & t_EDM);

template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> RankCompleteEDM(const EDMatrix<T,N,d>& t_EDM);


template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> OptSpace (const EDMatrix<T,N,d>& t_EDM);



//Starting definitions

//Constructors


//Default Constructor
template<typename T, uInt N, uInt d>
EDMatrix<T,N,d>:: EDMatrix() :
                             m_EDM(Eigen::Matrix<T,N,N>::Zero()),
                             m_Mask(Eigen::Matrix<int,N,N>:: Ones()),
                             m_Noise(Eigen::Matrix<T,N,N>::Zero()){
  /**
  *Dafult constructor: Initzilize the EDM as a Zero matrix, the mask as
  *1 matrix and the noise matrix as Zero matrix.
  *
  **/
}



//Copy Constructor
template<typename T, uInt N, uInt d>
EDMatrix<T,N,d>:: EDMatrix(const EDMatrix<T,N,d> & o_EDM):
                  m_EDM(o_EDM.getEDM()),
                  m_Mask(o_EDM.getMask()),
                  m_Noise(o_EDM.getNoise()){

            /**
            *Copy constructor
            **/
                  }
//Copy assignement operator

template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> EDMatrix<T,N,d>:: operator  =
                                    (const EDMatrix<T,N,d>& t_EDMatrix){

      setEDM(t_EDMatrix.getEDM());
      setMask(t_EDMatrix.getMask());
      setNoise(t_EDMatrix.getNoise());
      return *this;
                                    }


//Aritmetic operators
template<typename T, uInt N, uInt d>
EDMatrix<T,N,d> EDMatrix<T,N,d>:: operator + (const EDMatrix<T,N,d>&
                                                    t_EDM){
  /**
  *Perform the sum of two euclidean distance matrix
  *@perams this-> First Adder
  *@params t_EDM -> Second Adder
  *@return the sum of two euclidean distance matrix
  *
  **/
  EDMatrix<T,N,d> Sum (getEDM()+t_EDM.getEDM());
  return Sum;
}

template<typename T, uInt N, uInt d>
EDMatrix<T,N,d> EDMatrix<T,N,d>:: operator - (const EDMatrix<T,N,d>&
                                                    t_EDM){
  /**
  *Perform the sum of two euclidean distance matrix
  *@perams this-> First Adder
  *@params t_EDM -> Second Adder
  *@return the sum of two euclidean distance matrix
  *
  **/
  EDMatrix<T,N,d> Diff (getEDM()-t_EDM.getEDM());
  return Diff;
}



//Access operator
template<typename T, uInt N, uInt d>
T& EDMatrix<T,N,d>:: operator () (uInt Row , uInt Col)
{
  /**
  @return the element of the euclidean distance matrix in row Row and
  *column Col
  **/
  return m_EDM(Row,Col);

}



//
//sqrt operator
//
template<typename T, unsigned int N, unsigned int d>
EDMatrix<T,N,d> EDMatrix<T,N,d>:: sqrt() const {

  EDMatrix<T,N,d> R(getEDM(), getMask());
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
  *@see getNoise
  */
  return m_EDM;
}


template<typename T, unsigned int N, unsigned int d>
Eigen::Matrix<int,N,N> EDMatrix<T,N,d>:: getMask() const{
  /**
  * getter. @returns m_Mask
  *@see getEDM
  *@see getNoise
  */
  return m_Mask;
}


template<typename T, uInt N, uInt d>
Eigen::Matrix<T,N,N> EDMatrix<T,N,d>:: getNoise() const{
  /**
  *Getter
  *@returns m_Noise
  *@see getEDM
  *@see getMask
  **/
  return m_Noise;
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



template<typename T, uInt N, uInt d>
void EDMatrix<T,N,d>:: setNoise(const Eigen::Matrix<T,N,N>& t_Noise){

  /**
  * set the value of m_Noise to t_Noise
  * @param const Eigen:: Matrix<T,N,N>& t_Noise
  */
  m_Noise = t_Noise;
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




template<typename T, unsigned int N, unsigned int d>
void EDMatrix<T,N,d> :: trim(){

  //define some useful quantity
  unsigned int nEntires = 0;
  unsigned int eRow [N] = {0};
  //counts the number of entires
  //counts also the entires foe ach row and column

  for(unsigned int i=0 ; i<N; i++){
    for(unsigned int j=0; j<N; j++){
      nEntires += getMask()(i,j) == 1 ? 1:0;
      eRow[i]  += getMask()(i,j) == 1 ? 1:0;
    }
  }

  //define the threshold quantity
  double Tresh = 2*(double)nEntires / (double)N;



  //suppres over rapresented rows

  for(unsigned int i=0; i<N; i++){
    if(eRow[i]> Tresh ) {
      m_Mask.row(i) = Eigen::Matrix<int,1,N>::Zero();
      m_Mask.col(i) = Eigen::Matrix<int,N,1>::Zero();}
    else continue;
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
Eigen::Matrix<T,N,N> EDMatrix<T,N,d>:: gramm() const{

  /**
  * @returns the Gramm matrix G, given by:
  *@f$G  = XX^{T} = -0.5JDJ@f$
  *Where:
  * X is the point set (dxN matrix)
  * J is the geometric centring matrix @see gc_matrix
  * D is the Eucliden Distance Matrix (m_EDM).
  */

  Eigen::Matrix<T,N,N> D;
  D = hadamard<T,int,N,N>((getEDM()+getNoise()), getMask());
  return -0.5*gc_matrix()*D*gc_matrix();



}




//
//T method
//



template<typename T, unsigned int N, unsigned int d>
T EDMatrix<T,N,d>:: frobenius_norm(){

  /**
  * Perform the frombenius norm of the EDM with missing entires
  */
  Eigen::Matrix<T,N,N> D = hadamard<T,int,N,N>(getEDM()+getNoise(),getMask());
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

template<typename T,typename X, uInt Row, uInt Col>
Eigen::Matrix<T,Row,Col> hadamard(const Eigen::Matrix<T,Row,Col>& t_A,
                                  const Eigen::Matrix<X,Row,Col>& t_B)
{
/**
@return the hadamart product between matrix t_A and matrix t_B.
**/
  Eigen::Matrix<T,Row,Col> Result;
  for(unsigned int i=0; i<Row; ++i){
    for(unsigned int j=0; j<Col;++j){
      Result(i,j) = t_A(i,j)*t_B(i,j);
    }
  }
  return Result;

}


template<typename T, uInt Row, uInt Col>
T frobenius_norm(const Eigen::Matrix<T,Row,Col>& t_M){
/**
Compute the frobenius norm of the matrix t_M
**/
  T norm = 0;
  for(unsigned int i=0; i<Row; ++i){
    for(unsigned int j=0; j<Col; ++j){
      norm += std::pow(t_M(i,j),2);
    }
  }
  return std::sqrt(norm);
}


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

  Eigen::Matrix<T,N,N> EDM;
  EDM = hadamard<T,int,N,N>(t_EDM.getEDM()+t_EDM.getNoise(), t_EDM.getMask());

  Eigen::EigenSolver<Eigen::Matrix<T,N,N>> es(EDM);
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
Eigen::Matrix<T,d,N> ClassicalMDS(const EDMatrix<T,N,d>& t_EDM){

  /**
  * Perform the classical multi dimensional scaling
  *
  *@params t_EDM Euclidean distace matrix from which reconstruct the point set
  *
  @returns Nxd matrix which is the reconstructed point set up to a translation and rotation.
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

  return X.transpose();

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
  *Perform the matrix completion and the noise reduction by using the
  *rank proprieties of EDM.
  *@params t_EDM -> EDM to Complete
  *@params MAX_TOL -> treshold to establish the convergence of the
  *logoritm.
  *@params MAX_ITER -> max number of iter to do before stop the
  *algroitm, even if it doesn't converges.
  **/
  //quantity to use inside the cicle
  uInt count = 0;
  T conv;
  EDMatrix<T,N,d> D_old;
  Eigen::Matrix<T,N,N> D;
  D = hadamard<T,int,N,N>(t_EDM.getEDM()+t_EDM.getNoise(), t_EDM.getMask());
  //Compute the mean
  T mu = D.mean();
  //Create the matrix to work with
  Eigen::Matrix<T,N,N> A ;
  for(uInt R = 0; R< N ; R++){
    for(uInt C = R; C<N; C++ ){
      A(R,C) = t_EDM.getMask()(R,C) == 1 ? D(R,C) : mu;
      A(C,R) = A(R,C);
    }
  }

  EDMatrix<T,N,d> EDM(A);
  //start the cicle
  do{
    //Store the old value
    D_old = EDM;
    //Upgrade the value: EVTreshold
    EDM = EVTreshold(EDM,d+2);
    //Enforce known entires
    for(uInt i = 0; i<N; i++){
      for(uInt j=i; j<N; j++){
        EDM(i,j) = t_EDM.getMask()(i,j) == 1 ? D(i,j) : EDM(i,j) ;

        EDM(j,i) = EDM(i,j);
      }
    }
    //make hollow
    EDM.make_hollow();
    //make positive
    EDM.make_positive();

    //Check convergence
    conv = (D_old-EDM).frobenius_norm();
    count++;
  }while(count < MAX_ITER && conv > MAX_TOL);

  return EDM;
}



#endif
