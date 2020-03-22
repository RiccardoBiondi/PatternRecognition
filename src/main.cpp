#include <iostream>
#include <fstream>
#include <cstdlib>
#include <eigen\Eigen\Eigenvalues>

#include "EDMatrix.hpp"
#include "PointSet.hpp"
#include "Procrustes.hpp"

void write(const Eigen::Matrix<double,-1,-1>& t_M, std::string filename){


  std::ofstream o;
  o.open(filename);
  if(!o){
    std::cerr<<"file not opened"<<std::endl;
    std::exit(1);
  }
  for(unsigned int i = 0; i< t_M.rows(); i++){
    for(unsigned int j=0; j< t_M.cols(); j++){
      o<<t_M(i,j)<<'\t';
    }
    o<<'\n';
  }
  o.close();
}




int main (){

const unsigned int nPoints = 70; // number of points
const unsigned int dim = 2; // dimension of euclidean space
double h = 20. ; // side of hypercube in which points are distributes
double sigma = 10. ; // standard deviation of noise
double TOL = 1e-09 ;//MAX convergence tollerance
unsigned int ITER = 1000; //MAX number of iterations


EDMatrix<double,nPoints,dim> Complete  ;
Eigen::Matrix<double,dim, nPoints> REC;

std::clog<<"Generation of the pointset "<<std::endl;
PointSet<double, nPoints, dim> P (h);

std::clog<<"Adding noise and miss entires"<<std::endl;
P.AddMissEntires();
P.AddNoise(sigma);

std::clog<<"Performing matrix completion"<<std::endl;
Complete = RankCompleteEDM<double,nPoints,dim>(P.getEDM(), TOL, ITER);


std::clog<<"Performing the MultiDimensional Scaling"<<std::endl;
REC = ClassicalMDS(Complete);

std::clog<<"Starting procrutes analysis"<<std::endl;

Procrustes<double,nPoints,dim> OPA(P.getPointSet(), REC);

OPA.FullOPA();

std::clog<< "Writing data"<<std::endl;

write(OPA.getNewConfiguration().transpose(), "C:\\Users\\Riccardo\\github\\PatternRecognition\\file_1.dat");
write(OPA.getReference().transpose(), "C:\\Users\\Riccardo\\github\\PatternRecognition\\file_2.dat");


write(P.getPointSet().transpose(),"C:\\Users\\Riccardo\\github\\PatternRecognition\\file_3.dat");

write(REC.transpose(),"C:\\Users\\Riccardo\\github\\PatternRecognition\\file_4.dat");
std::clog<<"Finish"<<std::endl;











return 0;
}
