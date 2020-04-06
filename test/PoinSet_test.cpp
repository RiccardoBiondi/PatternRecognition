#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include<cstdlib>

#include "C:\Users\Riccardo\github\PatternRecognition\src\PointSet.hpp"

TEMPLATE_TEST_CASE("Testing Constructors","[one][template]",double){

  //Generate the size of the hypercube
  srand (time(NULL)); //set seed
  TestType h = rand() % 10 +1; // side between one and ten

  //create the PointSet object
  PointSet<TestType, 100, 2> P(h);
  //if the pointset is correctly created, the pointset matrix has 2 rows
  REQUIRE(P.getPointSet().cols()==100);
  REQUIRE(P.getPointSet().rows() == 2);
  //

  //if the EDM is correctly created, it is hollow, positive and respect
  //the tring inequality

  REQUIRE(P.getEDM().is_EDM());
  //the EDM has the default noise matrix and Mask
  REQUIRE(P.getEDM().getMask() == Eigen::Matrix<int,100,100>::Ones());
  REQUIRE(P.getEDM().getNoise() ==
                              Eigen::Matrix<TestType,100,100>::Zero());

  //
  //If I add the noise and mask
  //
  P.AddMissEntries();
  P.AddNoise(0.5);
  REQUIRE(P.getEDM().getMask() != Eigen::Matrix<int,100,100>::Ones());
  REQUIRE(P.getEDM().getNoise() !=
                              Eigen::Matrix<TestType,100,100>::Zero());
  REQUIRE(P.getEDM().getMask() == P.getEDM().getMask().transpose());
  REQUIRE(P.getEDM().getNoise() == P.getEDM().getNoise().transpose());
  //If I create a copy of the point set by using the copy Constructor
  //
  PointSet<TestType, 100,2> C(P);

  //EDM and point set are equal
  REQUIRE(C.getPointSet() == P.getPointSet());
  REQUIRE(C.getEDM().getEDM() == P.getEDM().getEDM());
}
