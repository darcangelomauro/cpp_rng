#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <armadillo>

#define MAT(i) get_mat(i)
#define EPS(i) get_eps(i)
#define GAMMA(i) get_gamma(i)

class Geom24
{
    public:
        

        // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR
        Geom24(int p, int q, int dim, double g2);
        Geom24(std::istream& in);
        Geom24(const Geom24& G);
        Geom24& operator=(const Geom24& G);
        ~Geom24();
        // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR



        // ============== GET METHODS
        int get_p() const {return p;}
        int get_q() const {return q;}
        int get_dim() const {return dim;}
        int get_nH() const {return nH;}
        int get_nL() const {return nL;}
        int get_nHL() const {return nHL;}
        int get_dim_omega() const {return dim_omega;}
        int get_eps(int i) const {return eps[i];}
        double get_g2() const {return g2;}
        arma::cx_mat get_mat(int i) const {return mat[i];}
        arma::cx_mat get_omega(int i) const {return omega[i];}
        // ============== GET METHODS

        

        // ============== ACTION METHODS
        double dirac2() const;
        double dirac4() const;
        double calculate_S() const; // using H and L decomposition
        double calculate_S_fromD() const; // using whole Dirac operator
        // ============== ACTION METHODS
        
        
        
        void derived_parameters();
        void init_gamma();
        void shuffle();
        std::istream& read_parameters(std::istream& in);



    protected:

        // ============== MATRICES
        // H and L matrices (all hermitian)
        arma::cx_mat* mat;   
        
        // epsilon: +1 for H, -1 for L
        int* eps;
        // ============== MATRICES
        
        
        
        
        // ============== BASE PARAMETERS
        // (p,q) numbers
        int p;
        int q;
        
        // size of H and L matrices
        int dim;

        // coupling constant
        double g2;
        // ============== BASE PARAMETERS

        


        // ============== DERIVED PARAMETERS
        // number of H and L matrices (nH and nL) and total number of matrices (nHL)
        int nH;
        int nL;
        int nHL;

        // size of gamma matrices
        int dim_omega;

        // omega matrices (all hermitian)
        arma::cx_mat* omega;   
        // ============== DERIVED PARAMETERS
        
};


std::ostream& operator<<(std::ostream& out, const Geom24& C);

#endif
