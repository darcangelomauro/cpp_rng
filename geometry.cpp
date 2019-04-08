#include <iostream>
#include "geometry.hpp"

using std::istream;
using std::ostream;
using std::clog;
using std::cerr;
using std::endl;

// Core constructors

Core::Core(int p_, int q_, int dim_)
    : p(p_), q(q_), dim(dim_)
{
    clog << "Geometry initialized with the following parameters:" << endl;
    clog << "(p, q, dim) = (" << p << ", " << q << ", " << dim << ")" << endl;
}

Core::Core(istream& in)
{
    read_parameters(in);
}

// Core copy constructor
Core::Core(const Core& C)
{
    dim = C.get_dim();
    p = C.get_p();
    q = C.get_q();
    nH = C.get_nH();
    nL = C.get_nL();
    nHL = C.get_nHL();
    dim_gamma = C.get_dim_gamma();
    
    clog << "Geometry initialized with the following parameters:" << endl;
    clog << "(p, q, dim) = (" << p << ", " << q << ", " << dim << ")" << endl;
}

// Core operator =
Core& Core::operator=(const Core& C)
{
    dim = C.get_dim();
    p = C.get_p();
    q = C.get_q();
    nH = C.get_nH();
    nL = C.get_nL();
    nHL = C.get_nHL();
    dim_gamma = C.get_dim_gamma();

    clog << "Geometry overwritten with the following parameters:" << endl;
    clog << "(p, q, dim) = (" << p << ", " << q << ", " << dim << ")" << endl;
    
    return *this;
}

// Core read parameters from istream
istream& Core::read_parameters(istream& in)
{
    if(in)
    {
        // read basic parameters
        in >> p >> q >> dim;
         
        // clear input stream state
        in.clear();
    }

    clog << "Geometry initialized with the following parameters:" << endl;
    clog << "(p, q, dim) = (" << p << ", " << q << ", " << dim << ")" << endl;

    return in;
}

ostream& operator<<(ostream& out, const Core& C)
{
    out << "Geometry (p, q, dim) = (" << C.get_p() << ", " << C.get_q() << ", " << C.get_dim() << ") ";

    return out;
}


