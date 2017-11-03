#ifndef WORKSPACE_H
#define WORKSPACE_H

#include <Eigen/Dense>


class workspace
{
    public:
        
        // MATRICES
        Eigen::MatrixXcd* H;
        Eigen::MatrixXcd* L;

        Eigen::MatrixXcd* gammaH;
        Eigen::MatrixXcd* gammaL;
        Eigen::MatrixXcd D;

        // PROPOSED MOVE AND ASSOCIATED INDEX
        Eigen::MatrixXcd MOVE;
        int move_idx;
        
        
        // methods to read private
        int get_p();
        int get_q();
        int get_dim();
        double get_g(int i);    // gets i-th coupling constant
        int get_nH();
        int get_nL();
        int get_n();
        int get_dim_gamma();
        double get_S();


        // ACTION FUNCTIONS
        // the powers of D that have to be computed are given by
        // the non-vanishing elements of the coupling constants array
        double calculate_S(); // H and L decomposition, calculate the whole action
        double calculate_delta_S(int i); // H and L decomposition, calculate difference when i-th matrix is changed
        double calculate_S_bruteforce();  // calculate S from the Dirac operator


        // METROPOLIS MOVE
        double propose_move(double scale);
        void accept_move(double delta_S);




    private:

        // INPUT PARAMETERS
        
        // (p,q) numbers
        int p;
        int q;
        
        // size of H and L matrices
        int dim;

        // coupling constants array
        double* g;

        // END INPUT PARAMETERS

        


        // DERIVED PARAMETERS

        // number of H and L matrices (nH and nL) and total number of matrices (n)
        int nH;
        int nL;
        int n;

        // size of gamma matrices
        int dim_gamma;

        // END DERIVED PARAMETERS
        
        
        
        // SCALAR OBSERVABLES
        double S;
        // END SCALAR OBSERVABLES
}













#endif
