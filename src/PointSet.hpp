#ifndef POINT_SET_HPP
#define POINT_SET_HPP

#include <random>

template<unsingned int N, unsigned int d>
class PointSet{

  private:
    Eigen::MatrixXd m_point_set;
    //add seed

  public:

    //Constructors
    PointSet ();
    PointSet(Eigen::MatrixXd t_pont_set);
    PointSet(const PointSet& o_Point_Set);

    //operator

    //getter and setter
    Eigen::MatrixXd get_point_set() const
    EDMatrix<N,d> get_EDM() const;
    Eigen::MatrixXd get_Noise() const ;
    Eigen::MatrixXd get_miss_entires() const;

    //Add set seed
    void setPoint(Eigen::MatrixXd t_point);

}



#endif
