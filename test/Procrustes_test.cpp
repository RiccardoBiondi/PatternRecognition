#define CATCH_CONFIG_MAIN
#include "catch.hpp"


#include "C:\Users\Riccardo\github\PatternRecognition\src\Procrustes.hpp"
#include "C:\Users\Riccardo\github\PatternRecognition\src\PointSet.hpp"


#include<iostream>
#include<typeinfo>
#include<complex>

TEMPLATE_TEST_CASE("tsting basic features", "[one][template]", double){
  //Generates two rando set of points
  PointSet<TestType,10,2> P(10.);
  PointSet<TestType,10,2> Q(10.);

  Procrustes<TestType, 10,2> OPA(P.getPointSet(), P.getPointSet());

  //Check if the constructor works well
  //Check also the getter
  REQUIRE(OPA.getReference() == P.getPointSet());
  REQUIRE(OPA.getConfiguration() == P.getPointSet());
  REQUIRE(OPA.getNewConfiguration()== P.getPointSet());
  REQUIRE(OPA.Rotation() == Eigen::Matrix<TestType,2,2>::Identity());
  REQUIRE(OPA.Scale() == 1.);

  
  //Perform the FullOPA()
  OPA.FullOPA();

  REQUIRE(typeid(OPA.Scale()).name() == typeid(TestType).name());
  REQUIRE(OPA.Rotation().rows() == 2);
  REQUIRE(OPA.Rotation().cols() == 2);
  REQUIRE(OPA.Rotation().determinant() == 1.);
  REQUIRE(OPA.Rotation().transpose() == OPA.Rotation().inverse());
  REQUIRE(OPA.getConfiguration().rows() == 2);
  REQUIRE(OPA.getConfiguration().cols() ==10);





}
