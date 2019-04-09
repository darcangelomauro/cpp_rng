#include <iostream>
#include "geometry.hpp"

using std::istream;
using std::ostream;
using std::clog;
using std::cerr;
using std::endl;

// Constructors

Geom24::Geom24(int p_, int q_, int dim_, double g2_)
    : p(p_), q(q_), dim(dim_), g2(g2_)
{

    // initialize derived parameters
    derived_parameters();

    // allocate and initialize H and L matrices to identity
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


    clog << "Geometry initialized with the following parameters:" << endl;
    clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;
}

Geom24::Geom24(istream& in)
{
    read_parameters(in);
    
    // initialize derived parameters
    derived_parameters();

    // allocate and initialize H and L matrices to identity
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
    
    clog << "Geometry initialized with the following parameters:" << endl;
    clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;

}

// Copy constructor
Geom24::Geom24(const Geom24& G)
{
    // copy parameters
    dim = G.get_dim();
    p = G.get_p();
    q = G.get_q();
    g2 = G.get_g2();
    nH = G.get_nH();
    nL = G.get_nL();
    nHL = G.get_nHL();
    dim_gamma = G.get_dim_gamma();


    // allocate and copy matrices
    mat = new arma::cx_mat [nHL];
    eps = new int [nHL];
    for(int i=0; i<nHL; i++)
    {
        mat[i] = G.MAT(i);
        eps[i] = G.EPS(i);
    }
    

    clog << "Geometry initialized with the following parameters:" << endl;
    clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;
}

// Operator =
Geom24& Geom24::operator=(const Geom24& G)
{
    dim = G.get_dim();
    p = G.get_p();
    q = G.get_q();
    g2 = G.get_g2();
    nH = G.get_nH();
    nL = G.get_nL();
    nHL = G.get_nHL();
    dim_gamma = G.get_dim_gamma();

    // delete, reallocate and copy matrices
    delete [] mat;
    delete [] eps;
    mat = new arma::cx_mat [nHL];
    eps = new int [nHL];
    for(int i=0; i<nHL; i++)
    {
        mat[i] = G.MAT(i);
        eps[i] = G.EPS(i);
    }
    
    
    clog << "Geometry overwritten with the following parameters:" << endl;
    clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;
    
    return *this;
}

// Destructor
Geom24::~Geom24()
{
    delete [] mat;
    delete [] eps;
}


// Read parameters from istream
istream& Geom24::read_parameters(istream& in)
{
    if(in)
    {
        // read basic parameters
        in >> p >> q >> dim >> g2;
         
        // clear input stream state
        in.clear();
    }

    return in;
}


void Geom24::derived_parameters()
{
    nH = 1;
    nL = 0;
    nHL = nH + nL;
    dim_gamma = 0;
}

ostream& operator<<(ostream& out, const Geom24& G)
{
    out << "Geometry (p, q, dim, g2) = (" << G.get_p() << ", " << G.get_q() << ", " << G.get_dim() << ", " << G.get_g2() << ") ";

    return out;
}


