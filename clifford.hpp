#ifndef CLIFFORD_HPP
#define CLIFFORD_HPP

#include <armadillo>

#define GAM(i) get_gamma(i)

class Cliff
{
    public:
        

        // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR
        
        Cliff(int p, int q);
        Cliff(std::istream& in);
        Cliff(const Cliff& G);
        Cliff& operator=(const Cliff& G);
        ~Cliff();

        // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR

        Cliff& operator*=(const Cliff& C);
        friend Cliff operator*(Cliff C1, const Cliff& C2){ C1*=C2; return C1; }


        // ============== GET METHODS
        
        int get_p() const { return p; }
        int get_q() const { return q; }
        arma::cx_mat get_gamma(int i) const { return gamma[i]; }

        // ============== GET METHODS

        

        std::istream& read_parameters(std::istream& in);


    private:

        int p;
        int q;

        arma::cx_mat* gamma;   
        arma::cx_mat chiral;   

};


std::ostream& operator<<(std::ostream& out, const Cliff& C);



#endif
