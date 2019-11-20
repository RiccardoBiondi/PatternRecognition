#ifndef MDS_HPP
#define MDS_HPP

/**
*This hpp file contains the implementation of a template class that
* contains methods to perform the multi dimensional scaling(), in order to
*reconstruct the reciprocal position of some point in space starting from an
*euclidean distances matrix(EDM).
*/





template<unsigned int N, unsigned int d>
class MDS{

private:
  Eigen::MatrixXd m_EDM;
  /**
  < NxN matrix of eculeidean distances. The entires type is T */
  void swap(Eigen::VectorXcd& eigenvalueis,
            Eigen::MatrixXcd& eigenvactors,
            unsigned int i, unsigned int j);
          /**This func swap  the i-th and j-th eigenvalues.
          */

  void sort(Eigen::VectorXcd& eigenvalues, Eigen::MatrixXcd& eigenvectors);
          /*
          This method sort the eigenvalues in defcreasing order
          */
  void set_zero(Eigen::VectorXcd& sorted_eigenvalues);
          /*
          This method set to 0 all the n-d smallest eigenvalues
          */
  void eigen_radq(Eigen::VectorXcd& eigenvalues);
          /*
           *This method is used to take the square root of the eigenvalues
           *belonghing form the MDS algorithm
           */

  public:
    MDS(); /*default costructor*/
    MDS(Eigen::MatrixXd t_EDM); /*Parametric
    *costructor*/
    MDS(const MDS & original);  /*Copy constructor*/
  //  ~MDS();

    //getter
    Eigen::MatrixXd getEDM () const;

    /*Start method declaration*/

    /** < This function implements the classical multidimensional scaling.
    *Will return an dxn matrix where d-> is the euclidean space dimension
    *(usually 2 or 3) and n is the number of points to determine.
    */
    Eigen:: MatrixXd classicalMDS();

};



//constructor definitions
template<unsigned int N, unsigned int d>
Eigen::MatrixXd MDS<N,d> :: getEDM() const{
  return m_EDM;

}


template<unsigned int N, unsigned int d>
MDS<N,d> :: MDS(){}


template< unsigned int N, unsigned int d>
MDS<N,d> :: MDS(Eigen::MatrixXd t_EDM):m_EDM(t_EDM){}


template<unsigned int N, unsigned int d>
MDS<N,d>:: MDS(const MDS & original){
        m_EDM = original.getEDM();

}


//private functions definitions

template<unsigned int N, unsigned int d>
void MDS<N,d>:: swap(Eigen::VectorXcd& eigenvalues,
                    Eigen::MatrixXcd& eigenvectors,
                    unsigned int i,unsigned int j){
              Eigen::VectorXcd vector_temporary = eigenvectors.col(i);
              std::complex<double> value_temporary = eigenvalues[i];
              eigenvectors.col(i) = eigenvectors.col(j);
              eigenvalues[i] = eigenvalues[j];
              eigenvectors.col(j) = vector_temporary;
              eigenvalues[j] = value_temporary;
            }


template<unsigned int N, unsigned int d>
void MDS<N,d> :: sort(Eigen::VectorXcd& eigenvalues,
                      Eigen::MatrixXcd& eigenvectors){
  unsigned int nSwap;
  do{
     nSwap = 0;
    for(unsigned int i=1; i<N; ++i){
      if(std::abs(eigenvalues[i-1]) < std::abs(eigenvalues[i])){
        swap(eigenvalues, eigenvectors,i-1,i);
        nSwap+=1;
      }
      else continue;
    }


  }while(nSwap >0);

}





template<unsigned int N, unsigned int d>
void MDS<N,d> :: set_zero(Eigen::VectorXcd& sorted_eigenvalues){
  for(unsigned int i=d; i<N; ++i) sorted_eigenvalues[i]=0;
}



template<unsigned int N, unsigned int d>
void MDS<N,d> :: eigen_radq(Eigen::VectorXcd& eigenvalues){

            for(unsigned int i=0; i<N; ++i)
                          eigenvalues[i] = std::sqrt(eigenvalues[i]);

}
//public metrods definitiocns


template<unsigned int N ,unsigned int d>
Eigen::MatrixXd MDS<N,d> :: classicalMDS(){

  //initzialize identity
  Eigen::MatrixXd I= Eigen::MatrixXd::Identity(N,N);

  //define a matrix with all 1
  Eigen::MatrixXd U = Eigen::MatrixXd::Ones(N,N);
  //define the centring Matrix
  Eigen::MatrixXd J = I-(double)1/N*U;
  //find the gram MatrixXd
  Eigen::MatrixXd G = -0.5*J*getEDM()*J;
  //diagonalize the gram Matrix
  Eigen::EigenSolver<Eigen::MatrixXd> es(G);
  //sort eigenvalues and corresponding eigenvectors
  Eigen::VectorXcd eigenval = es.eigenvalues();
  Eigen::MatrixXd eigenvec = es.eigenvectors().cast<double>();
  //sort(eigenval, eigenvec);
  //set to zero all the N-d lowest Eigenvalues
  set_zero(eigenval);
  //take the square root of the remaining Eigenvalues
  eigen_radq(eigenval);
  //eigenvalues to diagonal Matrix
  Eigen::MatrixXcd D = eigenval.asDiagonal();
  //return the resulting matrix
  return (double)eigenvec;



}






#endif
