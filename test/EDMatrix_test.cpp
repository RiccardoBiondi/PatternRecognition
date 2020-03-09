



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

    Eigen::Matrix<double,4,4> S;
    S<< 0  , 93, 82 , 133,
        93 , 0 , 52 , 60 ,
        82 , 52, 0  , 111,
        133, 60, 111, 0  ;



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
    WHEN("Square root of the distances square is required"){
      EDMatrix<double,4,2> D(A);
      THEN("The operator computes it correctrly"){
        REQUIRE(S ==D.sqrt());
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



SCENARIO("gc_matrix, gramm, hadamard and frobenius norm ", "[four]"){

  GIVEN("An EDM"){
    Eigen:: Matrix<double,4,4> A;
    A << 0.   , 8649., 6724.,17689.,
        8649. , 0.  , 2704., 3600.,
        6724. , 2704.,0.    ,12321.,
        17689., 3600.,12321., 0.;
    EDMatrix<double,4,2> D(A);

    WHEN("Given the correct results"){
      Eigen::Matrix<double,4,4> J;
      Eigen::Matrix<double,4,4> G;
      J<<   .75,-0.25,-0.25,-0.25,
          -0.25,.75  ,-0.25,-0.25,
          -0.25,-0.25,.75 ,-0.25,
          -0.25,-0.25,-0.25,.75;

      G = -0.5*J*A*J;


      THEN("The results coincide with the correct one"){

        REQUIRE(D.gc_matrix() == J);
        REQUIRE(D.gramm() == G);
      }

    WHEN("Set the mask, and give a control matrix an norm"){
      Eigen:: Matrix<int,4,4> W;
      Eigen:: Matrix<double,4,4> S;
      W << 0,1,1,0,
          1,1,1,1,
          1,1,0,0,
          0,1,0,0;
      S << 0   , 8649., 6724., 0    ,
          8649., 0.   , 2704., 3600.,
          6724., 2704., 0.   , 0.   ,
          0.   , 3600., 0.   , 0.   ;

      D.setMask(W);

      Eigen::Matrix<double,4,4> R;
      R = S*S.transpose();
      double norm = std::sqrt(R.trace());
          THEN("The method return the correct result"){
            REQUIRE(D.hadamard() == S);
            REQUIRE(D.frobenius_norm() == norm);
          }
    }}
  }
}




SCENARIO("classicalMDS","[five]"){

  GIVEN("An input EDM, resulting point set and a control matrix"){
    Eigen:: MatrixXd A(4,4);
      A << 0. , 8649., 6724.,17689.,
        8649. , 0.   , 2704., 3600.,
        6724. , 2704.,0.    ,12321.,
        17689., 3600.,12321., 0.;
    EDMatrix<double,4,2> D(A);

    Eigen::Matrix<double,4,2> Sol;

    Sol<< -62.8311,  32.9745,
           18.4029, -12.027 ,
          -24.9602, -39.7109,
           69.3884,  18.7634;




    WHEN("Perform the classicalMDS"){
      Eigen::Matrix<double, 4,2> Res;
      Res = ClassicalMDS<double,4,2>(D);

      THEN("We obtain the reconstructed point set, compared trhough frobenius norm"){
      Eigen::Matrix<double,4,2> Diff;
      Diff = Res + (-Sol);
      Eigen::Matrix<double,4,4> NS;
      NS = Diff*Diff.transpose();
      double norm = std::sqrt(NS.trace());
        //define the comparing functions
        REQUIRE(norm <1e-04);
      }
    }
  }
}



SCENARIO("EVTreshold test", "[six]"){

  GIVEN("An input EDM and a Control matrix"){
    Eigen:: MatrixXd A(4,4);
      A << 0. , 8649., 6724.,17689.,
        8649. , 0.   , 2704., 3600.,
        6724. , 2704.,0.    ,12321.,
        17689., 3600.,12321., 0.;
    EDMatrix<double,4,2> D(A);

    Eigen::Matrix<double,4,4> Sol;

    Sol<< 1703.59119256248,7291.03635562772,4073.23453154243,18643.5361281355,
7291.03635562771,1882.63548916683,4644.5785168113,2563.85641733326,
4073.23453154242,4644.5785168113,4161.69797519517,10895.0591521313,
18643.5361281355,2563.85641733326,10895.0591521313,629.527178256831;

    WHEN("Perform the EVTreshold"){


      EDMatrix<double,4,2> Res = EVTreshold(D,2);

      THEN("the frobenius norm of the diffrence ins lover then 1e-10"){

        Eigen::Matrix<double,4,4> Diff;
        Diff = Res.getEDM() + (-Sol);
        Eigen::Matrix<double,4,4> NS;
        NS = Diff*Diff.transpose();
        double norm = std::sqrt(NS.trace());
          //define the comparing functions
          REQUIRE(norm <1e-10);

      }

    }
  }

}
