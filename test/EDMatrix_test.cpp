#define CATCH_CONFIG_MAIN
#include "catch.hpp"

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






SCENARIO("getter, setter, costructors, operators","[one]"){
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
    WHEN("Create EDMatrix by default constructor and set the mask and the EDM"){
      EDMatrix<double,4,4> D;
      D.setEDM(B);
      D.setMask(W);
      THEN("check if setter works"){
        REQUIRE(D.getEDM()==B);
        REQUIRE(D.getMask()==W);
      }
    }
    WHEN("Create an EDMatrix by using the copy assignement operator"){
      EDMatrix<double,4,4> D (A);
      EDMatrix<double,4,4> E = D;

      THEN("The two EDM are equal"){
        REQUIRE(E.getEDM()==D.getEDM());
        REQUIRE(E.getMask() == D.getMask());
      }
    }
  }




}









SCENARIO("test the boolen methods", "[two]"){

  GIVEN("An EDM and a not EDM matrix"){

    Eigen:: MatrixXd A(4,4);
      A << 0. , 8649., 6724.,17689.,
        8649. , 0.   , 2704., 3600.,
        6724. , 2704.,0.    ,12321.,
        17689., 3600.,12321., 0.;

    Eigen:: MatrixXd B(4,4);
      B << 0. , 16., 17.,18.,
          16. , 3. , 19., 20.,
          7. , 49., -1 , 21.,
          13. , 25., 21., 0.;

    WHEN("Check the boolen methods on the EDM"){
      EDMatrix<double,4,4> D (A);
      THEN("All methods returns true"){
        REQUIRE(D.is_hollow());
        REQUIRE(D.is_positive());
        REQUIRE(D.is_symmetric());
        REQUIRE(D.is_triang_inh());
        REQUIRE(D.is_EDM());
      }


    }

    WHEN("Check boolean methods on the not EDM"){
      EDMatrix<double,4,4> D (B);
      THEN("All methods return false"){

        REQUIRE_FALSE(D.is_hollow());
        REQUIRE_FALSE(D.is_positive());
        REQUIRE_FALSE(D.is_symmetric());
        REQUIRE_FALSE(D.is_triang_inh());
        REQUIRE_FALSE(D.is_EDM());

      }
  }



}
}






SCENARIO("make_hollow and make_positive tests","[three]"){
  GIVEN("A not hollow and not positive matrix"){
    Eigen:: MatrixXd B(4,4);
      B << 0. , 16., 17.,18.,
          16. , 3. , 19., 20.,
          17. , 49., -1 , 21.,
          13. , 25., 21., 0.;
      EDMatrix<double,4,2> D(B);

      WHEN("use make hollow and make positive methods"){
        D.make_hollow();
        D.make_positive();

        THEN("D becomes hollow and positive"){
          REQUIRE(D.is_hollow());
          REQUIRE(D.is_positive());
        }
      }
  }
}






TEST_CASE("gc_matrix, gramm, hadamard and frobenius norm ", "[four]"){

//provide a starting matrix and the solution
  Eigen:: Matrix<double,4,4> A;
  Eigen:: Matrix<int,4,4> W;
  Eigen:: Matrix<double,4,4> S;
  A << 0.   , 8649., 6724.,17689.,
      8649. , 0.  , 2704., 3600.,
      6724. , 2704.,0.    ,12321.,
      17689., 3600.,12321., 0.;

  W << 0,1,1,0,
      1,1,1,1,
      1,1,0,0,
      0,1,0,0;

  S << 0   , 8649., 6724., 0    ,
      8649., 0.   , 2704., 3600.,
      6724., 2704., 0.   , 0.   ,
      0.   , 3600., 0.   , 0.   ;

    EDMatrix<double,4,2> D(A);

    //solutions

    Eigen::Matrix<double,4,4> J;
    Eigen::Matrix<double,4,4> G;
    double norm = 0;

    J<<   .75,-0.25,-0.25,-0.25,
        -0.25,.75  ,-0.25,-0.25,
        -0.25,-0.25,.75 ,-0.25,
        -0.25,-0.25,-0.25,.75;

    G = -0.5*J*A*J;


    //so now check te methods

    REQUIRE(D.gc_matrix() == J);
    REQUIRE(D.gramm() == G);

    //now set the mask

    D.setMask(W);
    REQUIRE(D.hadamard() == S);

    //compute the check value for the frobenius frobenius_norm
    Eigen::Matrix<double,4,4> R;
    R = S*S.transpose();
    norm = std::sqrt(R.trace());
    REQUIRE(D.frobenius_norm() == norm);



}
