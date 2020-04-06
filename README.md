# Euclidean Distance Matrix

In this repository I've implemented some useful tools to work with EDM. Here I've provide toll to generate a N point set in a d dimensional Euclidean space, find the EDM , complete the EDM, perform Multi Dimensional Scaling(MDS) and the Procrustes analysis.

Name | Contents
-----|---------
EDMatrix.hpp | Class to work with EDM, perform matrix completion and the multi-dimensional scaling
PointSet.hpp | Class that allow to generate a random point configuration
Procrustes.hpp | Class that allow to perform the full ordinary Procrustes analysis in order to align two points configuration
utilities.hpp | Implements some functions useful to work with EDM


## Requirements

To use this tools is necessary to include Eigen/Eigenvalues from eigen 3.7.0 (http://eigen.tuxfamily.org/dox/GettingStarted.html)
To perform test it is necessary t use Catch2(https://github.com/catchorg/Catch2)

## Usage

```c++
#include <Eigen\Eigenvalues>

#include "EDMatrix.hpp"
#include "PointSet.hpp"
#include "Procrustes.hpp"
#include "utilities.hpp"

int main(){

const unsigned int nPoints = 79; // number of points
const unsigned int dim = 2; // dimension of euclidean space
double h = 20. ; // side of hypercube in which points are distributes
double sigma = 10. ; // standard deviation of noise
double TOL = 1e-09 ;//MAX convergence tollerance
unsigned int ITER = 1000; //MAX number of iterations

EDMatrix<double,nPoints,dim> Complete  ;
Eigen::Matrix<double,dim, nPoints> REC;

//Generation of the point configuration
PointSet<double, nPoints, dim> P (h);

//Adding noise and miss entries
P.AddMissEntries();
P.AddNoise(sigma);

//Performing matrix completion
Complete = RankCompleteEDM<double,nPoints,dim>(P.getEDM(), TOL, ITER);

//Performing the MultiDimensional Scaling
REC = ClassicalMDS(Complete);

//Starting procrutes analysis
Procrustes<double,nPoints,dim> OPA(P.getPointSet(), REC);
OPA.FullOPA();


return 0;
}


```
