#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "clifford.hpp"

using namespace std;
using namespace arma;

// Constructors

Cliff::Cliff(int mode)
{
    //(1,0)
    if(mode==3)
    {
        p=1;
        q=0;
        dim_gamma = 1;

        gamma = new cx_mat [1];
        gamma[0] = cx_mat(1, 1, fill::eye);
        chiral = cx_mat(1, 1, fill::eye);
    }
    //(0,1)
    else if(mode==4)
    {
        p=0;
        q=1;
        dim_gamma = 1;
        
        cx_double z(0.,-1.);
        gamma = new cx_mat [1];
        gamma[0] = cx_mat(1, 1);
        chiral = cx_mat(1, 1, fill::eye);
        gamma[0](0,0) = z;
    }
    //(2,0)
    else if(mode==0)
    {
        p=2;
        q=0;
        dim_gamma = 2;
        
        cx_double z(0.,1.);

        gamma = new cx_mat [2];
        
        gamma[0] = cx_mat(2, 2, fill::zeros);
        gamma[0](0,0) = 1.;
        gamma[0](1,1) = -1.;
        
        gamma[1] = cx_mat(2, 2, fill::zeros);
        gamma[1](0,1) = 1.;
        gamma[1](1,0) = 1.;

        chiral = cx_mat(2, 2, fill::zeros);
        chiral(0,1) = z;
        chiral(1,0) = -z;
    }
    //(1,1)
    else if(mode==2)
    {
        p=1;
        q=1;
        dim_gamma = 2;
        
        gamma = new cx_mat [2];
        
        gamma[0] = cx_mat(2, 2, fill::zeros);
        gamma[0](0,0) = 1.;
        gamma[0](1,1) = -1.;
        
        gamma[1] = cx_mat(2, 2, fill::zeros);
        gamma[1](0,1) = 1.;
        gamma[1](1,0) = -1.;

        chiral = cx_mat(2, 2, fill::zeros);
        chiral(0,1) = 1;
        chiral(1,0) = 1;
    }
    //(0,2)
    else if(mode==1)
    {
        p=0;
        q=2;
        dim_gamma = 2;
        
        cx_double z(0.,1.);
        
        gamma = new cx_mat [2];
        
        gamma[0] = cx_mat(2, 2, fill::zeros);
        gamma[0](0,0) = z;
        gamma[0](1,1) = -z;
        
        gamma[1] = cx_mat(2, 2, fill::zeros);
        gamma[1](0,1) = 1.;
        gamma[1](1,0) = -1.;

        chiral = cx_mat(2, 2, fill::zeros);
        chiral(0,1) = 1.;
        chiral(1,0) = 1.;
    }
}

Cliff::Cliff(int p_, int q_)
    : p(p_), q(q_)
{
    //(1,0)
    if(p==1 && q==0)
    {
        dim_gamma = 1;

        gamma = new cx_mat [1];
        gamma[0] = cx_mat(1, 1, fill::eye);
        chiral = cx_mat(1, 1, fill::eye);
    }
    //(0,1)
    if(p==0 && q==1)
    {
        dim_gamma = 1;
        
        cx_double z(0.,-1.);
        gamma = new cx_mat [1];
        gamma[0] = cx_mat(1, 1);
        chiral = cx_mat(1, 1, fill::eye);
        gamma[0](0,0) = z;
    }
    //(2,0)
    if(p==2 && q==0)
    {
        dim_gamma = 2;
        
        cx_double z(0.,1.);

        gamma = new cx_mat [2];
        
        gamma[0] = cx_mat(2, 2, fill::zeros);
        gamma[0](0,0) = 1.;
        gamma[0](1,1) = -1.;
        
        gamma[1] = cx_mat(2, 2, fill::zeros);
        gamma[1](0,1) = 1.;
        gamma[1](1,0) = 1.;

        chiral = cx_mat(2, 2, fill::zeros);
        chiral(0,1) = z;
        chiral(1,0) = -z;
    }
    //(1,1)
    if(p==1 && q==1)
    {
        dim_gamma = 2;
        
        gamma = new cx_mat [2];
        
        gamma[0] = cx_mat(2, 2, fill::zeros);
        gamma[0](0,0) = 1.;
        gamma[0](1,1) = -1.;
        
        gamma[1] = cx_mat(2, 2, fill::zeros);
        gamma[1](0,1) = 1.;
        gamma[1](1,0) = -1.;

        chiral = cx_mat(2, 2, fill::zeros);
        chiral(0,1) = 1;
        chiral(1,0) = 1;
    }
    //(0,2)
    if(p==0 && q==2)
    {
        dim_gamma = 2;
        
        cx_double z(0.,1.);
        
        gamma = new cx_mat [2];
        
        gamma[0] = cx_mat(2, 2, fill::zeros);
        gamma[0](0,0) = z;
        gamma[0](1,1) = -z;
        
        gamma[1] = cx_mat(2, 2, fill::zeros);
        gamma[1](0,1) = 1.;
        gamma[1](1,0) = -1.;

        chiral = cx_mat(2, 2, fill::zeros);
        chiral(0,1) = 1.;
        chiral(1,0) = 1.;
    }
    //any other case
    else
        init_gamma();
}


// Copy constructor
Cliff::Cliff(const Cliff& C)
{
    // copy parameters
    p = C.get_p();
    q = C.get_q();
    dim_gamma = C.get_dim_gamma();


    // allocate and copy matrices
    gamma = new arma::cx_mat [p+q];
    for(int i=0; i<p+q; i++)
        gamma[i] = C.get_gamma(i);

    chiral = C.get_chiral();
}

// Operator =
Cliff& Cliff::operator=(const Cliff& C)
{
    p = C.get_p();
    q = C.get_q();
    dim_gamma = C.get_dim_gamma();

    // delete, reallocate and copy matrices
    delete [] gamma;
    gamma = new arma::cx_mat [p+q];
    for(int i=0; i<p+q; i++)
        gamma[i] = C.get_gamma(i);
    
    chiral = C.get_chiral();
    
    return *this;
}

// Destructor
Cliff::~Cliff()
{
    delete [] gamma;
}


void decomp(int p, int q, int* dec)
{
    if(p)
    {
        if(!(p%2))
        {
            dec[0] = p/2;
            dec[3] = 0;
        }
        else
        {
            dec[0] = (p-1)/2;
            dec[3] = 1;
        }
    }
    else
    {
        dec[0] = 0;
        dec[3] = 0;
    }

    if(q)
    {
        if(!(q%2))
        {
            dec[1] = q/2;
            dec[4] = 0;
        }
        else
        {
            dec[1] = (q-1)/2;
            dec[4] = 1;
        }
    }
    else
    {
        dec[1] = 0;
        dec[4] = 0;
    }

    if(dec[3] && dec[4])
    {
        dec[3] = 0;
        dec[4] = 0;
        dec[2] = 1;
    }
    else
        dec[2] = 0;
}


// init_gamma gets called only if p+q > 2
void Cliff::init_gamma()
{
    int dec[5];
    decomp(p, q, dec);

    vector<Cliff> vec;
    for(int i=0; i<5; ++i)
    {
        for(int j=0; j<dec[i]; ++j)
            vec.push_back( Cliff(i) );
    }
    //std::reverse(vec.begin(), vec.end());

    vector<Cliff>::const_iterator begin = vec.begin();
    vector<Cliff>::const_iterator end = vec.end();
    
    for(vector<Cliff>::const_iterator iter = begin; iter != end; ++iter)
        cout << (*iter) << endl;
    
    Cliff C1 = (*begin);

    for(vector<Cliff>::const_iterator iter = begin+1; iter != end; ++iter)
        C1 *= (*iter);

    dim_gamma = C1.get_dim_gamma(); 
    gamma = new cx_mat [p+q];
    for(int i=0; i<p+q; i++)
        gamma[i] = C1.get_gamma(i);
    chiral = C1.get_chiral();

    sort_gamma();
}

Cliff& Cliff::operator*=(const Cliff& C2)
{
    // store C2 frequently used variables
    int p2, q2, dim2;
    p2 = C2.get_p();
    q2 = C2.get_q();
    dim2 = C2.get_dim_gamma();
    

    // temporary variables to avoid overwriting on (*this)
    int p_, q_, dim_gamma_;
    cx_mat* gamma_;
    
    p_ = p + p2;
    q_ = q + q2;
    dim_gamma_ = dim_gamma * dim2;
    gamma_ = new cx_mat [p_+q_];
    

    // start computing product    
    cx_mat id2(dim2, dim2, fill::eye);

    for(int i=0; i<p+q; ++i)
        gamma_[i] = kron( gamma[i], id2 );
    for(int i=0; i<p2+q2; ++i)
        gamma_[p+q+i] = kron( chiral, C2.get_gamma(i) );


    // compute chirality
    int s = (q-p) % 8;
    int s2 = (q2-p2) % 8;
    if((s2 % 8) % 2)
    {
        cx_mat id1(dim_gamma, dim_gamma, fill::eye);
        chiral = kron(id1, C2.get_chiral());
    }
    else
        chiral = kron(chiral, C2.get_chiral());

    // overwrite on (*this)
    p = p_;
    q = q_;
    dim_gamma = dim_gamma_;
    delete [] gamma;
    gamma = new cx_mat [p+q];
    for(int i=0; i<p+q; ++i)
        gamma[i] = gamma_[i];


    // cleanup
    delete [] gamma_;

    return (*this);
}


    
ostream& operator<<(ostream& out, const Cliff& C)
{
    out << "Clifford (p, q) = (" << C.get_p() << ", " << C.get_q() << ") ";

    return out;
}

bool hermiticity(const cx_mat& M1, const cx_mat& M2)
{
    return(M2.is_hermitian() && !(M1.is_hermitian()));
}


void Cliff::sort_gamma()
{
    vector<cx_mat> vec;

    for(int i=0; i<p+q; ++i)
        vec.push_back(gamma[i]);

    sort(vec.begin(), vec.end(), hermiticity);

    vector<cx_mat>::const_iterator begin = vec.begin();
    vector<cx_mat>::const_iterator end = vec.end();
    for(vector<cx_mat>::const_iterator iter = begin; iter != end; ++iter)
        gamma[iter-begin] = (*iter);
}
