#ifndef SORTING_HPP
#define SORTING_HPP

namespace srt{

  void swap(Eigen::VectorXcd& t_eigenvalues,
            Eigen::MatrixXcd& t_eigenvactors,
            unsigned int i, unsigned int j);
        /**< This func swap  the i-th and j-th eigenvalues and the
         *corresponding eigenvectors. This function
        * is private and built to make more readable the code.
        *@param Eigen::VectorXcd& t_eigenvlues: A complex vector that
        * contains the eigenvalues to swap.
        * Eigen::MatrixXcd& t_eigenvactors: The corresponding eigenvactor
        *matrix.
        * unsigned int i, unsigned int j: The position of the two eigenvalues
        *to swap:
        */
  void sort(Eigen::VectorXcd& t_eigenvalues,
            Eigen::MatrixXcd& t_eigenvectors,
          unsigned int N);
                /*
                This method sort the eigenvalues in defcreasing order
                */
  void set_zero(Eigen::VectorXcd& t_eigenvalues,
                unsigned int r, unsigned int N);
                /*
                This method set to 0 all the n-d smallest eigenvalues
                */
  void eigen_sqrt(Eigen::VectorXcd& t_eigenvalues, unsigned int N);
                /*
                 *This method is used to take the square root of the eigenvalues
                 *belonghing form the MDS algorithm
                 */


  void swap(Eigen::VectorXcd& t_eigenvalues,
            Eigen::MatrixXcd& t_eigenvectors,
            unsigned int i,unsigned int j){

    Eigen::VectorXcd vector_temporary = t_eigenvectors.col(i);
    std::complex<double> value_temporary = t_eigenvalues[i];
    t_eigenvectors.col(i) = t_eigenvectors.col(j);
    t_eigenvalues[i] = t_eigenvalues[j];
    t_eigenvectors.col(j) = vector_temporary;
    t_eigenvalues[j] = value_temporary;
                            }



  void  sort(Eigen::VectorXcd& t_eigenvalues,
            Eigen::MatrixXcd& t_eigenvectors,
            unsigned int N){
            unsigned int nSwap;
            do{
                 nSwap = 0;
                for(unsigned int i=1; i<N; ++i){
                    if(std::abs(t_eigenvalues[i-1]) <
                        std::abs(t_eigenvalues[i])){
                              swap(t_eigenvalues, t_eigenvectors,i-1,i);
                              nSwap+=1;
                              }
                    else continue;
                                }


                }while(nSwap >0);

                            }





  void set_zero(Eigen::VectorXcd& t_eigenvalues,
                unsigned int r, unsigned int N){

            for(unsigned int i=r; i< N ; ++i) t_eigenvalues[i]=0;
                            }






  void eigen_sqrt (Eigen::VectorXcd& t_eigenvalues, unsigned int N){

          for(unsigned int i=0; i<N; ++i)
            t_eigenvalues[i] = std::sqrt(t_eigenvalues[i]);

                            }









}
#endif
