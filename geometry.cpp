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
using arma::cx_mat;
using arma::cx_double;

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

void Geom24::shuffle()
{
    for(int i=0; i<nHL; i++)
    {
        arma::cx_mat temp(dim, dim, arma::fill::randu);
        mat[i] = temp + temp.t();
    }
}


double Geom24::dirac2() const
{
    double res = 0.;
    for(int i=0; i<nHL; i++)
    {
        double trMM = trace(mat[i]*mat[i]);
        double trM = trace_herm(mat[i]);

        res += (dim*trMM + eps[i]*trM*trM);
    }

    return 2.*dim_gamma*res;
}



double Geom24::dirac4() const
{
    double res = 0.;
    int* i = new int [4];


    for(i[3]=0; i[3]<nHL; i[3]++)
    {
        for(i[2]=0; i[2]<nHL; i[2]++)
        {
            for(i[1]=0; i[1]<=i[2]; i[1]++)
            {
                for(i[0]=0; i[0]<=i[3]; i[0]++)
                {
                    if(i[0]==i[1] && i[1]==i[2] && i[2]==i[3])
                    {
                        // compute matrix products
                        cx_mat M0M0 = mat[i[0]]*mat[i[0]];
                        cx_mat M0M0M0 = M0M0*mat[i[0]];
                        cx_mat M0M0M0M0 = M0M0M0*mat[i[0]];

                        // compute traces
                        double trM0 = trace(mat[i[0]]).real();
                        double trM0M0 = trace(M0M0).real();
                        double trM0M0M0 = trace(M0M0M0).real();
                        double trM0M0M0M0 = trace(M0M0M0M0).real();
                        
                        // add to total
                        double temp = 0;

                        // tr4 term
                        temp += dim*2.*trM0M0M0M0;

                        // tr3tr1 term
                        temp += 8.*eps[i[0]]*trM0M0M0*trM0;

                        // tr2tr2 term
                        temp += 6.*trM0M0*trM0M0;
                        
                        res += dim_gamma*temp;
                    }


                    
                    else if(i[0] == i[3] && i[1] == i[2])
                    {
                        // alloc matrix products (trace needed)
                        cx_mat M0M0 = mat[i[0]]*mat[i[0]];
                        cx_mat M1M1 = mat[i[1]]*mat[i[1]];
                        cx_mat M0M1 = mat[i[0]]*mat[i[1]];
                        cx_mat M1M0M0M1 = mat[i[1]]*M0M0*mat[i[1]];
                        cx_mat M0M1M0 = M0M1*mat[i[0]];
                        cx_mat M1M0M1 = mat[i[1]]*M0M1;
                            
                        // compute traces
                        double trM0 = trace(mat[i[0]]).real();
                        double trM1 = trace(mat[i[1]]).real();
                        double trM0M0 = trace(M0M0).real();
                        double trM1M1 = trace(M1M1).real();
                        double trM0M1 = trace(M0M1).real();
                        double trM0M1M0 = trace(M0M1M0).real();
                        double trM1M0M1 = trace(M1M0M1).real();
                        double trM1M0M0M1 = trace(M1M0M0M1).real();
                        
                        // add to total
                        double temp = 0;

                        // tr4 term
                        temp += dim*2.*trM1M0M0M1;

                        // tr3tr1 term
                        temp += 4.*eps[i[0]]*trM1M0M1*trM0;
                        temp += 4.*eps[i[1]]*trM0M1M0*trM1;

                        // tr2tr2 term
                        temp += 4.*eps[i[0]]*eps[i[1]]*trM0M1*trM0M1;
                        temp += 2.*trM0M0*trM1M1;
                        
                        res += dim_gamma*temp;
                    }
                    

                    else
                    {
                        cx_double cliff = gamma_table[4][i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                        
                        if(cliff.real() != 0. || cliff.imag() != 0.)
                        {
                            // compute matrix products
                            cx_mat M0M1 = mat[i[0]]*mat[i[1]];
                            cx_mat M0M2 = mat[i[0]]*mat[i[2]];
                            cx_mat M0M3 = mat[i[0]]*mat[i[3]];
                            cx_mat M1M2 = mat[i[1]]*mat[i[2]];
                            cx_mat M1M3 = mat[i[1]]*mat[i[3]];
                            cx_mat M2M3 = mat[i[2]]*mat[i[3]];
                            cx_mat M0M1M2 = M0M1*mat[i[2]];
                            cx_mat M0M1M3 = M0M1*mat[i[3]];
                            cx_mat M0M2M3 = M0M2*mat[i[3]];
                            cx_mat M1M2M3 = mat[i[1]]*M2M3;
                            cx_mat M0M1M2M3 = mat[i[0]]*M1M2M3;

                            // compute traces
                            cx_double trM0M1M2M3 = trace(M0M1M2M3);
                            cx_double trM0M1M2 = trace(M0M1M2);
                            cx_double trM0M1M3 = trace(M0M1M3);
                            cx_double trM0M2M3 = trace(M0M2M3);
                            cx_double trM1M2M3 = trace(M1M2M3);
                            double trM0M1 = trace(M0M1).real();
                            double trM0M2 = trace(M0M2).rel();
                            double trM0M3 = trace(M0M3).real();
                            double trM1M2 = trace(M1M2).real();
                            double trM1M3 = trace(M1M3).real();
                            double trM2M3 = trace(M2M3).real();
                            double trM0 = trace(mat[i[0]]).real();
                            double trM1 = trace(mat[i[1]]).real();
                            double trM2 = trace(mat[i[2]]).real();
                            double trM3 = trace(mat[i[3]]).real();
                            

                            // add to total

                            // tr4 terms
                            cx_double T1 = trM0M1M2M3 + eps[i[0]]*eps[i[1]]*eps[i[2]]*eps[i[3]]*conj(trM0M1M2M3);
                            res += 2.*dim*(cliff*T1).real();

                            // tr3tr1 terms

                            cx_double T3 = eps[i[3]]*trM0M1M2 + eps[i[0]]*eps[i[1]]*eps[i[2]]*conj(trM0M1M2);
                            T3 = trM3*T3;

                            cx_double T4 = eps[i[2]]*trM0M1M3 + eps[i[0]]*eps[i[1]]*eps[i[3]]*conj(trM0M1M3);
                            T3 += trM2*T4;

                            cx_double T5 = eps[i[1]]*trM0M2M3 + eps[i[0]]*eps[i[2]]*eps[i[3]]*conj(trM0M2M3);
                            T3 += trM1*T5;

                            cx_double T6 = eps[i[0]]*trM1M2M3 + eps[i[1]]*eps[i[2]]*eps[i[3]]*conj(trM1M2M3);
                            T3 += trM0*T6;

                            res += 2.*(cliff*T3).real();

                            // tr2tr2 terms
                            double T7 = trM0M1*trM2M3*(eps[i[0]]*eps[i[1]] + eps[i[2]]*eps[i[3]]);
                            double T8 = trM0M2*trM1M3*(eps[i[0]]*eps[i[2]] + eps[i[1]]*eps[i[3]]);
                            double T9 = trM0M3*trM1M2*(eps[i[0]]*eps[i[3]] + eps[i[1]]*eps[i[2]]);

                            res += 2.*cliff.real()*(T7+T8+T9);
                        }
                    }
                }
            }
        }
    }
    
    for(i[3]=0; i[3]<nHL; i[3]++)
    {
        for(i[1]=0; i[1]<nHL; i[1]++)
        {
            for(i[2]=0; i[2]<i[1]; i[2]++)
            {
                for(i[0]=0; i[0]<i[3]; i[0]++)
                {
                    cx_double cliff = gamma_table[4][i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                    
                    if(cliff.real() != 0. || cliff.imag() != 0.)
                    {
                        // compute matrix products
                        cx_mat M0M1 = mat[i[0]]*mat[i[1]];
                        cx_mat M0M2 = mat[i[0]]*mat[i[2]];
                        cx_mat M0M3 = mat[i[0]]*mat[i[3]];
                        cx_mat M1M2 = mat[i[1]]*mat[i[2]];
                        cx_mat M1M3 = mat[i[1]]*mat[i[3]];
                        cx_mat M2M3 = mat[i[2]]*mat[i[3]];
                        cx_mat M0M1M2 = M0M1*mat[i[2]];
                        cx_mat M0M1M3 = M0M1*mat[i[3]];
                        cx_mat M0M2M3 = M0M2*mat[i[3]];
                        cx_mat M1M2M3 = M1M2*mat[i[3]];
                        cx_mat M0M1M2M3 = mat[i[0]]*M1M2M3;

                        // compute traces
                        cx_double trM0M1M2M3 = trace(M0M1M2M3);
                        cx_double trM0M1M2 = trace(M0M1M2);
                        cx_double trM0M1M3 = trace(M0M1M3);
                        cx_double trM0M2M3 = trace(M0M2M3);
                        cx_double trM1M2M3 = trace(M1M2M3);
                        double trM0M1 = trace(M0M1).real();
                        double trM0M2 = trace(M0M2).real();
                        double trM0M3 = trace(M0M3).real();
                        double trM1M2 = trace(M1M2).real();
                        double trM1M3 = trace(M1M3).real();
                        double trM2M3 = trace(M2M3).real();
                        double trM0 = trace(mat[i[0]]).real();
                        double trM1 = trace(mat[i[1]]).real();
                        double trM2 = trace(mat[i[2]]).real();
                        double trM3 = trace(mat[i[3]]).real();
                        
                        // add to total

                        // tr4 terms
                        cx_double T1 = trM0M1M2M3 + eps[i[0]]*eps[i[1]]*eps[i[2]]*eps[i[3]]*conj(trM0M1M2M3);
                        res += dim*2.*(cliff*T1).real();

                        // tr3tr1 terms
                        cx_double T3 = eps[i[3]]*trM0M1M2 + eps[i[0]]*eps[i[1]]*eps[i[2]]*conj(trM0M1M2);
                        T3 = trM3*T3;

                        cx_double T4 = eps[i[2]]*trM0M1M3 + eps[i[0]]*eps[i[1]]*eps[i[3]]*conj(trM0M1M3);
                        T3 += trM2*T4;

                        cx_double T5 = eps[i[1]]*trM0M2M3 + eps[i[0]]*eps[i[2]]*eps[i[3]]*conj(trM0M2M3);
                        T3 += trM1*T5;

                        cx_double T6 = eps[i[0]]*trM1M2M3 + eps[i[1]]*eps[i[2]]*eps[i[3]]*conj(trM1M2M3);
                        T3 += trM0*T6;

                        res += 2.*(cliff*T3).real();

                        // tr2tr2 terms
                        double T7 = trM0M1*trM2M3*(eps[i[0]]*eps[i[1]] + eps[i[2]]*eps[i[3]]);
                        double T8 = trM0M2*trM1M3*(eps[i[0]]*eps[i[2]] + eps[i[1]]*eps[i[3]]);
                        double T9 = trM0M3*trM1M2*(eps[i[0]]*eps[i[3]] + eps[i[1]]*eps[i[2]]);

                        res += 2.*cliff.real()*(T7+T8+T9);
                    }
                }
            }
        }
    }

    delete [] i;
    return res;
}




double Geom24::calculate_S() const
{
    return g2*dirac2() + dirac4();
}
