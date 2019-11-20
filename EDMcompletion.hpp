#ifndef EDMCOMPLETION_HPP
#define EDMCOMPLETION_HPP

template<unsigned int N, unsigned int d>
class Completion{

  private:
    Eigen:: MatrixXd m_EDM;

    //some private functions
    Eigen::MatrixXd evThreshold(const Eigen::MatrixXd &D,unsigned int r);

  public:
    Completion();
    Completion(Eigen::MatrixXd t_EDM);
    Completion(const Completion<N,d> & original);
    ~Completion();

    Eigen::MatrixXd getEDM() const;


};


template<unsigned int N, unsigned int d>
Eigen::MatrixXd Completion<N,d> :: getEDM() const {return m_EDM;}


template<unsigned int N, unsigned int d>
Completion<N,d>:: Completion(){};

template<unsigned int N, unsigned int d>
Completion<N,d>:: Completion(Eigen::MatrixXd t_EDM): m_EDM(t_EDM){};

template<unsigned int N, unsigned int d>
Completion<N,d>:: Completion(const Completion<N,d>& original ){
  m_EDM = oriinal.getEDM();
}

template<unsigned int N ,unsigned int d>
Completion<N,d>:: ~Completion(){};




#endif
