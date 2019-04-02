#include "geometry.hpp"

using std::istream;
using std::clog;
using std::endl;
using namespace Geometry;

istream& read_parameters(istream& in)
{
    p = q = dim = ncc = 0;

    if(in)
    {
        // read basic parameters
        in >> p >> q >> dim >> ncc;
         
        if(ncc < 1)
            throw domain_error("there must be at least one coupling constant");

        // initialize coupling constant array to zero
        g = new double [ncc];
        for(int i=0; i<ncc; i++)
            g[i] = 0; 
        
        // fill coupling constant array from input stream
        for(int i=0; i<ncc; i++)
            in >> g[i]; 

        // check that at least one coupling constant is non-zero
        bool check = 0
        for(int i=0; i<ncc; i++)
            if(g[i]) check = 1; 
        if(!check)
            throw domain_error("at least one coupling constant must be non-zero");

        // clear input stream state
        in.clear();
    }

    if((p < 0) || (q < 0))
        throw domain_error("p and q must be positive integers");
    if(!p && !q)
        throw domain_error("(p,q) = (0,0) is not a valid geometry");
    if(dim < 1)
        throw domain_error("matrix dimension must be positive non-zero integer");

    
    clog << "Geometry initialized with the following parameters" << endl;
    clog << "(p,q) = (" << p << "," << q << ")" << endl;
    clog << "dim = " << dim << endl;
    for(int i=0; i<ncc; i++)
        clog << "g" << 2*(i+1) << " = " << g[i] << endl;


    return in;
}

