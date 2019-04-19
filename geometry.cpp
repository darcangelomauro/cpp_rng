#include <iostream>
#include <cmath>
#include <vector>
#include "geometry.hpp"

using std::vector;
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
    int n = p+q;

    int* gamma = new int [n];
    for(int i=0; i<n; i++)
        gamma[i] = i+1;

    nH = 0;
    nL = 0;

	int  count = pow(2, n);
	// The outer for loop will run 2^n times (the number of all possible subsets).
	// Here variable i will act as a binary counter
	for (int i=0; i<count; i++)
	{
        vector<int> vec;
		// The inner for loop will run n times, As the maximum number of elements a set can have is n
		// This loop will generate a subset
		for (int j=0; j<n; j++)
		{
			// This if condition will check if jth bit in binary representation of i is set or not
			// if the value of (i & (1 << j)) is greater than 0, include arr[j] in the current subset
			// otherwise exclude arr[j]
			if ((i & (1 << j)) > 0)
                vec.push_back(gamma[j]);
		}
        
        // Now print subset if it has odd number of elements
        int k = vec.size();
        if(k % 2)
        {
            vector<int>::const_iterator end(vec.end());
            vector<int>::const_iterator begin(vec.begin());
            int first = 0;
            int second = 0;
            for(vector<int>::const_iterator l = begin; l != end; ++l)
            {
                if((*l) > p)
                {
                    first = k - (l-begin);
                    break;
                }
            }
            second = k*(k-1)/2;

            if ((first+second) % 2)
                ++nL;
            else
                ++nH;
        }
	}

    delete [] gamma;

    nHL = nH+nL;

    
    if(n % 2)
        dim_gamma = pow(2, (n-1)/2);
    else
        dim_gamma = pow(2, n/2);

}


ostream& operator<<(ostream& out, const Geom24& G)
{
    out << "Geometry (p, q, dim, g2) = (" << G.get_p() << ", " << G.get_q() << ", " << G.get_dim() << ", " << G.get_g2() << ") ";

    return out;
}


