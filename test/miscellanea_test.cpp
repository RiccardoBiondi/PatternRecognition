
#include "C:\Users\Riccardo\github\PatternRecognition\src\miscellanea.hpp"

#include<typeinfo>
#include<complex>
/*
 *  miscellanea_test
 *
 *  ----------------------------------------------------------
 * This filecontains all the test performed to check EDMatrix.hpp functions and classes
 *  Requiriments:
 *
 * Catch v2.11.0
 */

#define CATCH_CONFIG_MAIN
#include "catch.hpp"




TEMPLATE_TEST_CASE("testing sorting algorithms","[one][template]", float){


  // decalration of the vector and matrix to sort
  Eigen::Matrix<TestType, 4,1> eigenval;
  Eigen::Matrix<TestType, 4,4> eigenvec;
  eigenval<< 4,3,2,1;
  eigenvec<< 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15;


  Eigen::Matrix<TestType, 4,1> val_result;
  Eigen::Matrix<TestType, 4,4> vec_result;

  val_result<< 4,2,3,1;
  vec_result<< 0,2,1,3,4,6,5,7,8,10,9,11,12,14,13,15;

  Eigen::Matrix<TestType, 4,1> sort_val_result = eigenval;
  Eigen::Matrix<TestType, 4,4> sort_vec_result = eigenvec;




  //if use swap

  swap<TestType,4>(eigenval, eigenvec, 1,2);

      REQUIRE(eigenval==val_result);
      REQUIRE(eigenvec==vec_result);

  //if use sort
  sort<TestType, 4> (eigenval, eigenvec);

      REQUIRE(eigenval == sort_val_result);
      REQUIRE(eigenvec == sort_vec_result);






}





TEMPLATE_TEST_CASE("test set zero and sqrt function", "[two][template]", double){


  //initial definitions
  Eigen::Matrix<TestType, 7,1> vec;
  Eigen::Matrix<TestType, 7,1> control;
  vec << 1,4,9,16, 4,5,6;
  control << 1,2,3,4,0,0,0;

  set_zero<TestType, 7>(vec, 4);
  REQUIRE(vec[4]==0);
  REQUIRE(vec[5]==0);
  REQUIRE(vec[6]==0);

  eigen_sqrt<TestType,7>(vec);

  REQUIRE(vec == control);


}






TEMPLATE_TEST_CASE("casting test", "[three][template]", double){

  Eigen::Matrix<std::complex<TestType>, 4,1> vec;
  vec << 1,2,3,4;

  //if no cast
  REQUIRE(typeid(vec).name() == typeid(Eigen::Matrix<std::complex<TestType>,4,1>).name());

  //if cast
  REQUIRE(typeid(cast_real<TestType,4,1>(vec)).name()== typeid(Eigen::Matrix<TestType,4,1>).name());

}
