#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <armadillo>

#define MAT(i) geometry::get_mat(i)
#define EPS(i) geometry::get_eps(i)

class geometry
{
    public:

        // METHODS TO READ PRIVATE
        arma::cx_mat get_mat(int i);   // returns matrix mat[i]
        int get_eps(int i);   // returns eps[i]
        int get_p();
        int get_q();
        int get_dim();
        double get_g(int i);    // gets i-th coupling constant
        int get_npow();
        int get_nH();
        int get_nL();
        int get_n();
        int get_dim_gamma();


        // ACTION FUNCTIONS
        // (they depend on npow)
        double calculate_S();                       // using H and L decomposition
        double calculate_S_fromD();                 // using whole Dirac operator
        arma::cx_double calculate_S();              // using H and L decomposition (both real and imaginary part)
        arma::cx_double calculate_S_fromD();        // using whole Dirac operator (both real and imaginary part)



    private:

        // ********** MATRICES
        
        // H and L matrices (all hermitian)
        arma::cx_mat* mat;   
        
        // epsilon: +1 for H, -1 for L
        int* eps;

        // ********** END MATRICES
        
        
        
        
        // ********** INPUT PARAMETERS
        
        // (p,q) numbers
        int p;
        int q;
        
        // size of H and L matrices
        int dim;

        // coupling constants array (g), and dimension of array (npow)
        double* g;
        int npow;

        // ********** END INPUT PARAMETERS

        


        // ********** DERIVED PARAMETERS

        // number of H and L matrices (nH and nL) and total number of matrices (n)
        int nH;
        int nL;
        int n;

        // size of gamma matrices
        int dim_gamma;

        // ********** END DERIVED PARAMETERS
}













#endif
