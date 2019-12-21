
#include "C:\Users\Riccardo\github\PatternRecognition\src\EDMatrix.hpp"
/*
 *  EDMAtrix_test
 *
 *  ----------------------------------------------------------
 * This filecontains all the test performed to check EDMatrix.hpp functions and classes
 *  Requiriments:
 *
 * Catch v2.11.0
 */

#define CATCH_CONFIG_MAIN
#include "catch.hpp"



SCENARIO("getter, setter, costructors","[one]"){
  GIVEN("two different EDM A, B and a mask W"){
    Eigen:: MatrixXd A(4,4);
      A << 0. , 8649., 6724.,17689.,
        8649. , 0.   , 2704., 3600.,
        6724. , 2704.,0.    ,12321.,
        17689., 3600.,12321., 0.;

    Eigen:: MatrixXd B(4,4);
      B << 0. , 16., 17.,18.,
          16. , 0. , 19., 20.,
          17. , 19., 0. , 21.,
          18. , 20., 21., 0.;


    Eigen:: MatrixXi W(4,4);

      W<< 0,1,1,0,
          1,1,1,1,
          1,1,0,0,
          0,1,0,0;


    WHEN("Create an EDMatrix with A"){

            EDMatrix<double,4,4> D (A);
      THEN("Check the getter"){
        REQUIRE(D.getEDM()==A);
        REQUIRE(D.getMask()==Eigen::MatrixXi::Ones(4,4));
      }
    }
    WHEN("Create EDMamtrix by default constructor and set the mask and the EDM"){
      EDMatrix<double,4,4> D;
      D.setEDM(B);
      D.setMask(W);
      THEN("check if setter woerks"){
        REQUIRE(D.getEDM()==B);
        REQUIRE(D.getMask()==W);
      }
    }
  }




}
