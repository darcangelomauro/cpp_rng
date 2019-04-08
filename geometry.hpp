#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

class Core
{
    public:
        
        Core(int p, int q, int dim);
        Core(std::istream& in);
        Core(const Core& C);
        Core& operator=(const Core& C);
        ~Core(){};


        std::istream& read_parameters(std::istream& in);

        int get_p() const {return p;}
        int get_q() const {return q;}
        int get_dim() const {return dim;}
        int get_nH() const {return nH;}
        int get_nL() const {return nL;}
        int get_nHL() const {return nHL;}
        int get_dim_gamma() const {return dim_gamma;}


    protected:
        
        // size of H and L matrices
        int dim;
    
        // number of H and L matrices (nH and nL) and total number of matrices (nHL)
        int nH;
        int nL;
        int nHL;


    private:

        // (p,q) numbers
        int p;
        int q;
        
        
        // size of gamma matrices
        int dim_gamma;
};


std::ostream& operator<<(std::ostream& out, const Core& C);

/*
class Geometry
{
    public:

        // CONSTRUCTOR / DESTRUCTOR
        Geometry(std::istream& in);
        ~Geometry();


        // METHODS TO INITIALIZE PARAMETERS
        std::istream& read_parameters(std::istream& in);
        void derived_parameters(int p, int q);


        // METHODS TO READ PRIVATE
        arma::cx_mat get_mat(int i) const;   // returns matrix mat[i]
        int get_eps(int i) const;   // returns eps[i]
        int get_p() const;
        int get_q() const;
        int get_dim() const;
        double get_g(int i) const;    // gets i-th coupling constant
        int get_ncc() const;
        int get_nH() const;
        int get_nL() const;
        int get_n() const;
        int get_dim_gamma() const;


        // ACTION FUNCTIONS
        // (they depend on ncc)
        double calculate_S() const;                       // using H and L decomposition
        double calculate_S_fromD() const;                 // using whole Dirac operator
        
        * not sure this makes sense (the formulas are already manifestly real)
         * also: functions cannot be overloaded by return type, so this will clash
         * with previous calculate_S functions
         *
        arma::cx_double calculate_S() const;              // using H and L decomposition (both real and imaginary part)
        arma::cx_double calculate_S_fromD() const;        // using whole Dirac operator (both real and imaginary part)
        *


    private:

        // ********** MATRICES
        
        // H and L matrices (all hermitian)
        arma::cx_mat* mat;   
        
        // epsilon: +1 for H, -1 for L
        int* eps;

        // ********** END MATRICES
        
        
        
        
        // ********** FUNDAMENTAL PARAMETERS
        
        // (p,q) numbers
        int p;
        int q;
        
        // size of H and L matrices
        int dim;

        // coupling constants array (g), and dimension of array (ncc)
        // NOTE: these are even powers only. odd powers not implemented
        //
        double* g;
        int ncc;

        // ********** END FUNDAMENTAL PARAMETERS

        


        // ********** DERIVED PARAMETERS

        // number of H and L matrices (nH and nL) and total number of matrices (nHL)
        int nH;
        int nL;
        int nHL;

        // size of gamma matrices
        int dim_gamma;

        // ********** END DERIVED PARAMETERS
}






Geometry::Geometry(std::istream& in)
try
{
    // read parameters from input stream
    read_parameters(in);

    // find derived parameters
    derived_parameters(p, q);

    // initialize H and L matrices to identity
    mat = new arma::cx_mat [nHL];
    eps = new int [nHL];
    for(int i=0; i<nHL; i++)
    {
        if(i<nH)
            eps[i] = 1;
        else
            eps[i] = -1;

        mat[i].eye(dim, dim); 
    }
}
catch(const std::domain_error& e)
{
    std::cout << e.what();
}


Geometry::~Geometry()
{
    delete [] g;
    delete [] mat;
    delete [] eps;
}



*/







#endif
