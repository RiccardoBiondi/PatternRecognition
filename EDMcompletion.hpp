#ifndef EDMCOMPLETION_HPP
#define EDMCOMPLETION_HPP



template<unsigned int N, unsigned int d>
class Completion{
  /*
    Insert here class description

  */
  private:
    Eigen:: MatrixXd m_EDM;
    Eigen:: MatrixXd m_Mask;
    //some private functions
    Eigen::MatrixXd evThreshold(const Eigen::MatrixXd &D,unsigned int r);

    void swap(Eigen::VectorXcd& t_eigenvalues,
              Eigen::MatrixXcd& t_eigenvactors,
              unsigned int i, unsigned int j);

    void sort(Eigen::VectorXcd& t_eigenvalues,
              Eigen::MatrixXcd& t_eigenvectors);

    void set_zero(Eigen::VectorXcd& t_eigenvalues,unsigned int r);

    void eigen_sqrt(Eigen::VectorXcd& t_eigenvalues);


  public:
    Completion();
    Completion(Eigen::MatrixXd t_EDM, Eigen::MatrixXd t_Mask);
    Completion(const Completion<N,d> & original);
    ~Completion();

    Eigen::MatrixXd getEDM() const;
    Eigen::MatrixXd getMask() const;


};


template<unsigned int N, unsigned int d>
void Completion<N,d> :: swap(Eigen::VectorXcd& t_eigenvalues,
                                          Eigen::MatrixXcd& t_eigenvectors,
                                          unsigned int i, unsigned int j){
      Eigen::VectorXcd vector_temporary=t_eigenvectors.col(i);
      std::complex<double> value_temporary = t_eigenvalues[i];
      t_eigenvectors.col(i) = t_eigenvectors.col(j);
      t_eigenvalues[i] = t_eigenvalues[j];
      t_eigenvectors.col(j) = vector_temporary;
      t_eigenvalues[j] = value_temporary;}



template<unsigned int N, unsigned int d>
void Completion<N,d>:: sort(Eigen::VectorXcd& t_eigenvalues,
                                          Eigen::MatrixXcd& t_eigenvectors){

            unsigned int nSwap;
            do{
                nSwap = 0;
                for(unsigned int i=1; i<N; ++i){
                    if(std::abs(t_eigenvalues[i-1]) <
                      std::abs(t_eigenvalues[i])){
                        swap(t_eigenvalues, t_eigenvectors,i-1,i);
                        nSwap+=1;
                        }
                    else continue;
                          }
            }while(nSwap >0);
                                          }


template<unsigned int N, unsigned int d>
void Completion<N,d>:: set_zero(Eigen::VectorXcd&t_eigenvalues,
                                              unsigned int r){
      for(unsigned int i=r; i< N ; ++i) t_eigenvalues[i]=0;
}



template<unsigned int N, unsigned int d>
void Completion<N,d>::eigen_sqrt(Eigen::VectorXcd&t_eigenvalues){
  for(unsigned int i=0; i<N; ++i)
    t_eigenvalues[i] = std::sqrt(t_eigenvalues[i]);
}

template<unsigned int N, unsigned int d>
Completion<N,d>:: Completion(){};

template<unsigned int N, unsigned int d>
Completion<N,d>:: Completion(Eigen::MatrixXd t_EDM, Eigen::t_Mask):
                              m_EDM(t_EDM), m_Mask(t_Mask){};

template<unsigned int N, unsigned int d>
Completion<N,d>:: Completion(const Completion<N,d>& original ){
  m_EDM = original.getEDM();
  t_Mask = original.getMask();
}

template<unsigned int N ,unsigned int d>
Completion<N,d>:: ~Completion(){}



//getter
template<unsigned int N, unsigned int d>
Eigen::MatrixXd Completion<N,d>:: getEDM() const {
  /*
  Getter: return the Euclidean Distance Matrix. Is const declared.
  */
  return m_EDM;
}

template<unsigned int N, unsigned int d>
Eigen::MatrixXd Completion<N,d>:: getMask() const{ return m_Mask;}
// definiton of private function

template<unsigned int N, unsigned int d>
Eigen::MatrixXd Completion<N,d>::evThreshold(const Eigen::MatrixXd &D,
                                                unsigned int r){

      Eigen::EigenSolver<Eigen::MatrixXd> es(getEDM());
      Eigen::VectorXcd eigenval = es.eigenvalues();
      Eigen::MatrixXcd eigenvec = es.eigenvectors();

      sort(eigenval, eigenvec);
      set_zero(eigenval,r);
      return eigenvec*eigenval.asDiagonal()*eigenvec.transpose();
                                                }






#endif
