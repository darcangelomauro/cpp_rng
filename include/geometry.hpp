#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <armadillo>
#include <gsl/gsl_rng.h>

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
        arma::cx_mat get_mom(int i) const {return mom[i];}
        arma::cx_mat get_omega(int i) const {return omega[i];}
        arma::cx_double get_omega_table_4(int i) const {return omega_table_4[i];}
        // ============== GET METHODS

        

        // ============== ACTION METHODS
        arma::cx_mat build_dirac() const;
        double dirac2() const;
        double dirac4() const;
        double dirac4_old() const; // DEPRECATED
        double compute_A4(const int&, const int&, const int&, const int&) const;
        double compute_A2(const int&, const int&) const;
        double compute_A(const int&) const;
        double calculate_S() const; // using H and L decomposition
        double calculate_S_old() const; // DEPRECATED using H and L decomposition
        double calculate_S_from_dirac() const; // using whole Dirac operator
        // ============== ACTION METHODS
        

        // ============== DERIVATIVE METHODS
        arma::cx_mat compute_B4(const int&, const int&, const int&, const int&, const double&, const bool&) const;
        arma::cx_mat compute_B2(const int&, const int&) const;
        arma::cx_mat compute_B(const int&) const;
        arma::cx_mat der_dirac4(const int&, const bool&) const;
        arma::cx_mat der_dirac2(const int&) const;
        arma::cx_mat der_dirac24(const int&, const bool&) const;
        // ============== DERIVATIVE METHODS

        // ============== HAMILTONIAN METHODS
        void sample_mom(gsl_rng*);
        double calculate_K() const;
        double calculate_H() const;
        void leapfrog(const int&, const double&);
        double HMC_core(const int&, const double&, gsl_rng*, double*, double*);
        double HMC(const int&, double&, const int&, gsl_rng*, std::ostream&);
        double HMC(const int&, double&, const int&, gsl_rng*);
        double HMC(const int&, const double&, const int&, const int&, gsl_rng*, std::ostream&, std::ostream&);
        // ============== HAMILTONIAN METHODS
        
        // ============== METROPOLIS METHODS
        double delta2(const int&, const int&, const int&, const arma::cx_double&);
        double delta4(const int&, const int&, const int&, const arma::cx_double&);
        double delta24(const int&, const int&, const int&, const arma::cx_double&);
        double MMC_core(const double&, gsl_rng*, double*, double*);
        double MMC(double&, const int&, const bool&, gsl_rng*, std::ostream&, std::ostream&);
        void delta24_debug(const double&, const int&, gsl_rng*, std::ostream&);
        // ============== METROPOLIS METHODS
        
        void derived_parameters();
        void shuffle(gsl_rng*);
        std::istream& read_mat(std::istream& in);
        void reverse_mom();
        void init_omega_table_4();
        void print_omega_table_4() const;
        std::istream& read_parameters(std::istream& in);



    protected:

        // ============== MATRICES
        // H and L matrices (all hermitian)
        arma::cx_mat* mat;   
        
        // conjugate momenta
        arma::cx_mat* mom;
        
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

        // omega 4-product table
        arma::cx_double* omega_table_4;
        // ============== DERIVED PARAMETERS
        
};

std::vector<int> base_conversion(int dec, const int& base, const int& max);
std::ostream& operator<<(std::ostream& out, const Geom24& C);


#endif
