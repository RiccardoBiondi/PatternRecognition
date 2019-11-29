#ifndef MDS_HPP
#define MDS_HPP

/**
*This hpp file contains the implementation of a template class that
* contains methods to perform the multi dimensional scaling(), in order to
*reconstruct the reciprocal position of some point in space starting from an
*euclidean distances matrix(EDM).
*/


namespace MDS{


template<unsigned int N, unsigned int d>
class MultiDimensionalScaling{
  /**
  This template class allow to perform the Multididimensional scaling of an
  *euclidean distance matrix, in order to obtain a matrix with the coordinate
  *of the points.
  *Template @params:
  *      unsigned int d -> Dimension of the euclidean space from which the
  *                         distances are mesured.
  *      unsigned int N -> Dimension of the NxN Euclidean Distance Matrix
  */

private:
  Eigen::MatrixXd m_EDM;
  //! The euclidean distance matrix
  void swap(Eigen::VectorXcd& t_eigenvalues,
            Eigen::MatrixXcd& t_eigenvactors,
            unsigned int i, unsigned int j);
          /**< This func swap  the i-th and j-th eigenvalues and the
           *corresponding eigenvectors. This function
          * is private and built to make more readable the code.
          *@param Eigen::VectorXcd& t_eigenvlues: A complex vector that
          * contains the eigenvalues to swap.
          * Eigen::MatrixXcd& t_eigenvactors: The corresponding eigenvactor
          *matrix.
          * unsigned int i, unsigned int j: The position of the two eigenvalues *to swap:
          */

  void sort(Eigen::VectorXcd& t_eigenvalues,Eigen::MatrixXcd& t_eigenvectors);
          /*
          This method sort the eigenvalues in defcreasing order
          */
  void set_zero(Eigen::VectorXcd& t_eigenvalues);
          /*
          This method set to 0 all the n-d smallest eigenvalues
          */
  void eigen_sqrt(Eigen::VectorXcd& t_eigenvalues);
          /*
           *This method is used to take the square root of the eigenvalues
           *belonghing form the MDS algorithm
           */

  public:
    MultiDimensionalScaling(); /*default costructor*/
    MultiDimensionalScaling(Eigen::MatrixXd t_EDM); /*Parametric costructor*/
    MultiDimensionalScaling(const MultiDimensionalScaling & t_original);  /*Copy constructor*/
   ~MultiDimensionalScaling();

    //getter
    Eigen::MatrixXd getEDM () const;

    /*Start method declaration*/


    Eigen:: MatrixXcd classicalMDS();

  //  Eigen::MatrixXd AlternatingDescent(const Eigen::MatrixXd &W); //TODO

    //Semdefinite relaxation

};



//constructor definitions
template<unsigned int N, unsigned int d>
Eigen::MatrixXd MultiDimensionalScaling<N,d> :: getEDM() const{
  return m_EDM;

}


template<unsigned int N, unsigned int d>
MultiDimensionalScaling<N,d> :: MultiDimensionalScaling(){}


template< unsigned int N, unsigned int d>
MultiDimensionalScaling<N,d> :: MultiDimensionalScaling(Eigen::MatrixXd t_EDM):m_EDM(t_EDM){}


template<unsigned int N, unsigned int d>
MultiDimensionalScaling<N,d>:: MultiDimensionalScaling
(const MultiDimensionalScaling & t_original){
        m_EDM = t_original.getEDM();}

template<unsigned int N, unsigned int d>
MultiDimensionalScaling<N,d> :: ~MultiDimensionalScaling(){};
//private functions definitions

template<unsigned int N, unsigned int d>
void MultiDimensionalScaling<N,d>:: swap(Eigen::VectorXcd& t_eigenvalues,
                    Eigen::MatrixXcd& t_eigenvectors,
                    unsigned int i,unsigned int j)
                    {
              Eigen::VectorXcd vector_temporary = t_eigenvectors.col(i);
              std::complex<double> value_temporary = t_eigenvalues[i];
              t_eigenvectors.col(i) = t_eigenvectors.col(j);
              t_eigenvalues[i] = t_eigenvalues[j];
              t_eigenvectors.col(j) = vector_temporary;
              t_eigenvalues[j] = value_temporary;
            }


template<unsigned int N, unsigned int d>
void MultiDimensionalScaling<N,d> :: sort(Eigen::VectorXcd& t_eigenvalues,
                      Eigen::MatrixXcd& t_eigenvectors){
  unsigned int nSwap;
  do{
     nSwap = 0;
    for(unsigned int i=1; i<N; ++i){
      if(std::abs(t_eigenvalues[i-1]) < std::abs(t_eigenvalues[i])){
        swap(t_eigenvalues, t_eigenvectors,i-1,i);
        nSwap+=1;
      }
      else continue;
    }


  }while(nSwap >0);

}





template<unsigned int N, unsigned int d>
void MultiDimensionalScaling<N,d> :: set_zero(Eigen::VectorXcd& t_eigenvalues){
  for(unsigned int i=d; i<N; ++i) t_eigenvalues[i]=0;
}



template<unsigned int N, unsigned int d>
void MultiDimensionalScaling<N,d> :: eigen_sqrt
                    (Eigen::VectorXcd& t_eigenvalues){

            for(unsigned int i=0; i<N; ++i)
                          t_eigenvalues[i] = std::sqrt(t_eigenvalues[i]);

}
//public metrods definitiocns


template<unsigned int N ,unsigned int d>
Eigen::MatrixXcd MultiDimensionalScaling<N,d> :: classicalMDS(){
  /** < This function implements the
  *classical multidimensional scaling.
  *@return \f$\mathbb{C}^{d\timesn}\f$
  */

  //initzialize identity
  Eigen::MatrixXd I= Eigen::MatrixXd::Identity(N,N);
  //define a matrix with all 1
  Eigen::MatrixXd U = Eigen::MatrixXd::Ones(N,N);
  //find the gram MatrixXd
  Eigen::MatrixXd G = -0.5*(I-(double)1/N*U)*getEDM()*(I-(double)1/N*U);
  //diagonalize the gram Matrix
  Eigen::EigenSolver<Eigen::MatrixXd> es(G);
  Eigen::VectorXcd eigenval = es.eigenvalues();
  Eigen::MatrixXcd eigenvec = es.eigenvectors();
  //sort the eigenvalues and the corresponding eigenvectors
  sort(eigenval, eigenvec);
  //set to zero all the N-d lowest Eigenvalues
  set_zero(eigenval);
  //take the square root of the remaining Eigenvalues
  eigen_sqrt(eigenval);

  Eigen::MatrixXcd A = eigenvec*eigenval.asDiagonal();
  Eigen::MatrixXcd X (N,d);
  for(unsigned int i=0; i<d; ++i)
                  X.col(i) = A.col(i);

  //return the resulting matrix
  return X;



}

}




#endif
