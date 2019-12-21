#ifndef EDMATRIX_HPP
#define EDMATRIX_HPP



template<unsigned int N, unsigned int d>
class EDMatrix : Eigen::MatrixXd {

  private:
    Eigen::MatrixXd m_EDM;
    EIgen::MatrixXd m_Mask;

  public:
    //Constructors
    EDMatrix();
    EDMatrix(t_EDM, t_Mask = Ones(N,N));
    EDMatrix(const EDMatrix& o_Matrix);

    //operator

    //getters, setter
    MatrixXd getEDM() const;
    MatrixXd getMask() const;

    void setMask();



    //Member functions

    bool is_hollow();
    bool is_symmetric();
    bool is_positive();
    bool is_triang_inh();
    bool is_EDM();

    void make_hollow();
    void make_positive();
    void trim();

    MatrixXd gc_matrix();
    MatrixXd gramm();

    double frobenius_norm();





};





//Declaratin of externam methods
template<unsigned int N, unsigned int d>
EDMatrix<N,d> EVTreshold(EDMatrix<N,d> t_EDM);

template<unsigned int N, unsigned int d>
MatrixXd ClassicalMDS(EDMAtrix<N,d> t_EDM);





#endif
