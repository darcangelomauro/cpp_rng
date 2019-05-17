#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "geometry.hpp"
#include "clifford.hpp"

using namespace std;
using namespace arma;

// Constructors

Geom24::Geom24(int p_, int q_, int dim_, double g2_)
    : p(p_), q(q_), dim(dim_), g2(g2_)
{
    // initialize derived parameters
    derived_parameters();

    // allocate and initialize H and L matrices to identity
    mat = new arma::cx_mat [nHL];
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    for(int i=0; i<nHL; i++)
    {
        if(i<nH)
            eps[i] = 1;
        else
            eps[i] = -1;

        mat[i].eye(dim, dim); 
        mom[i].eye(dim, dim); 
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
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    for(int i=0; i<nHL; i++)
    {
        if(i<nH)
            eps[i] = 1;
        else
            eps[i] = -1;

        mat[i].eye(dim, dim); 
        mom[i].eye(dim, dim); 
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
    dim_omega = G.get_dim_omega();


    // allocate and copy matrices
    mat = new arma::cx_mat [nHL];
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    omega = new arma::cx_mat [nHL];
    for(int i=0; i<nHL; i++)
    {
        mat[i] = G.get_mat(i);
        mom[i] = G.get_mom(i);
        eps[i] = G.get_eps(i);
        omega[i] = G.get_omega(i);
    }
    
    int nHL4 = pow(nHL, 4);
    omega_table_4 = new cx_double [nHL4];
    for(int i=0; i<nHL4; i++)
        omega_table_4[i] = G.get_omega_table_4(i);

    
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
    dim_omega = G.get_dim_omega();

    // delete, reallocate and copy matrices
    delete [] mat;
    delete [] mom;
    delete [] omega;
    delete [] eps;
    mat = new arma::cx_mat [nHL];
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    omega = new arma::cx_mat [nHL];
    for(int i=0; i<nHL; i++)
    {
        mat[i] = G.get_mat(i);
        mom[i] = G.get_mat(i);
        eps[i] = G.get_eps(i);
        omega[i] = G.get_omega(i);
    }
    
    delete [] omega_table_4;
    int nHL4 = pow(nHL, 4);
    omega_table_4 = new cx_double [nHL4];
    for(int i=0; i<nHL4; i++)
        omega_table_4[i] = G.get_omega_table_4(i);
    
    clog << "Geometry overwritten with the following parameters:" << endl;
    clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;
    
    return *this;
}

// Destructor
Geom24::~Geom24()
{
    delete [] mat;
    delete [] mom;
    delete [] eps;
    delete [] omega;
    delete [] omega_table_4;
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

    // create a type (p, q) clifford module
    Cliff C(p, q);
    vector<cx_mat> gamma = C.get_gamma();
    dim_omega = C.get_dim_gamma();
    
    vector<cx_mat> herm;
    vector<cx_mat> anti;

    for(int i=0; i<p; ++i)
        herm.push_back(gamma[i]);
    for(int i=0; i<q; ++i)
        anti.push_back(cx_double(0, 1)*gamma[p+i]);

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
                vec.push_back(j);
		}
        
        // Now calculate and push product if it has odd number of gammas
        int k = vec.size();
        if(k % 2 && k != 1)
        {
            vector<int>::const_iterator begin(vec.begin());
            vector<int>::const_iterator end(vec.end());
            cx_mat M = gamma.at(*begin);

            for(vector<int>::const_iterator iter = vec.begin() + 1; iter != end; ++iter)
                M *= gamma.at((*iter));

            if(M.is_hermitian())
                herm.push_back(M);
            else
                anti.push_back(cx_double(0, 1)*M);
        }
	}

    nH = herm.size();
    nL = anti.size();
    nHL = nH+nL;

    omega = new cx_mat [nHL];
    for(int i=0; i<nH; ++i)
        omega[i] = herm[i];
    for(int i=0; i<nL; ++i)
        omega[nH+i] = anti[i];

    init_omega_table_4();
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
        arma::cx_mat temp(dim, dim, arma::fill::randn);
        mat[i] = (temp + temp.t())/(2*dim*dim);
    }
}

istream& Geom24::read_mat(istream& in)
{
    if(in)
    {
        // loop on matrices
        for(int i=0; i<nHL; ++i)
        {
            // loop on indices
            for(int j=0; j<dim; ++j)
            {
                for(int k=0; k<dim; ++k)
                {
                    double x, y;
                    in >> x >> y;
                    mat[i](j,k) = cx_double(x, y);
                }
            }
        }
         
        // clear input stream state
        in.clear();
    }

    return in;
}


void Geom24::reverse_mom()
{
    for(int i=0; i<nHL; ++i)
        mom[i] *= -1.;
}

void Geom24::init_omega_table_4()
{
    omega_table_4 = new cx_double [nHL*nHL*nHL*nHL];

    for(int i=0; i<nHL; ++i)
    {
        for(int j=0; j<nHL; ++j)
        {
            for(int k=0; k<nHL; ++k)
            {
                for(int l=0; l<nHL; ++l)
                    omega_table_4[l + nHL*(k + nHL*(j + nHL*i))] = trace(omega[i]*omega[j]*omega[k]*omega[l]);
            }
        }
    }
}

vector<int> base_conversion(int dec, const int& base, const int& max)
{
    vector<int> rem;

    while(dec)
    {
        rem.push_back(dec % base);
        dec /= base;
    }

    for(int i=rem.size(); i<max; i++)
        rem.push_back(0);

    reverse(rem.begin(), rem.end());

    return rem;
}


void Geom24::print_omega_table_4() const
{
    const int n = pow(nHL, 4);

    for(int i=0; i<n; ++i)
    {
        cx_double z = omega_table_4[i];
        if(z != cx_double(0., 0.))
        {
            int e = 1;
            vector<int> prod = base_conversion(i, nHL, 4);
            vector<int>::const_iterator end(prod.end());
            for(vector<int>::const_iterator iter = prod.begin(); iter != end; ++iter)
            {
                cout << (*iter);
                e *= eps[(*iter)];
            }
            cout << " " << omega_table_4[i] << e << endl;
        }
    }
}



cx_mat Geom24::build_dirac() const
{
    // initialize dirac op to zero
    int dim_dirac = dim*dim*dim_omega;
    cx_mat dirac(dim_dirac, dim_dirac, fill::zeros);

    static cx_mat id(dim, dim, fill::eye);
    for(int i=0; i<nHL; ++i)
    {
        cx_mat bracket = kron(mat[i], id) + eps[i]*kron(id, mat[i].st());
        dirac += kron(omega[i], bracket);
    }

    return dirac;
}

double Geom24::calculate_S_from_dirac() const
{
    cx_mat dirac = build_dirac();
    cx_mat dirac2 = dirac*dirac;
    double trdirac2 = trace(dirac2).real();
    double trdirac4 = trace(dirac2*dirac2).real();
    return g2*trdirac2 + trdirac4;
}

double Geom24::dirac2() const
{
    double res = 0.;
    for(int i=0; i<nHL; ++i)
    {
        double tr2 = trace(mat[i]*mat[i]).real();
        double tr1 = trace(mat[i]).real();

        res += (dim*tr2 + eps[i]*tr1*tr1);
    }

    return 2.*dim_omega*res;
}




double Geom24::calculate_S() const
{
    return g2*dirac2() + dirac4();
}


double Geom24::compute_A4(const int& i1, const int& i2, const int& i3, const int& i4) const
{
    // epsilon factor
    int e = eps[i1]*eps[i2]*eps[i3]*eps[i4];

    // if e=-1, then [1+*e] becomes 2i*imag
    // and the clifford part is guaranteed to
    // be pure imaginary
    if(e<0)
    {
        // clifford product
        double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].imag(); 

        if(cliff != 0.)
        {
            // base matrix products
            cx_mat M1M2 = mat[i1]*mat[i2];
            cx_mat M1M3 = mat[i1]*mat[i3];
            cx_mat M1M4 = mat[i1]*mat[i4];
            cx_mat M2M3 = mat[i2]*mat[i3];
            cx_mat M2M4 = mat[i2]*mat[i4];
            cx_mat M3M4 = mat[i3]*mat[i4];

            // traces
            double tr1234 = trace(M1M2*M3M4).imag();
            double tr234 = trace(M2M3*mat[i4]).imag();
            double tr134 = trace(M1M3*mat[i4]).imag();
            double tr124 = trace(M1M2*mat[i4]).imag();
            double tr123 = trace(M1M2*mat[i3]).imag();
            double tr1 = trace(mat[i1]).real();
            double tr2 = trace(mat[i2]).real();
            double tr3 = trace(mat[i3]).real();
            double tr4 = trace(mat[i4]).real();
            
            // compute sum
            double res = dim*tr1234;
            res += eps[i1]*tr1*tr234;
            res += eps[i2]*tr2*tr134;
            res += eps[i3]*tr3*tr124;
            res += eps[i4]*tr4*tr123;

            return -2*cliff*res;
            // NOTE: this minus here comes from the i in cliff
            // and the i coming from 2i*imag
        }
        else
            return 0.;
    }
    else
    {
        // clifford product
        double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].real(); 

        if(cliff != 0.)
        {
            // base matrix products
            cx_mat M1M2 = mat[i1]*mat[i2];
            cx_mat M1M3 = mat[i1]*mat[i3];
            cx_mat M1M4 = mat[i1]*mat[i4];
            cx_mat M2M3 = mat[i2]*mat[i3];
            cx_mat M2M4 = mat[i2]*mat[i4];
            cx_mat M3M4 = mat[i3]*mat[i4];

            // traces
            double tr1234 = trace(M1M2*M3M4).real();
            double tr234 = trace(M2M3*mat[i4]).real();
            double tr134 = trace(M1M3*mat[i4]).real();
            double tr124 = trace(M1M2*mat[i4]).real();
            double tr123 = trace(M1M2*mat[i3]).real();
            double tr12 = trace(M1M2).real();
            double tr34 = trace(M3M4).real();
            double tr13 = trace(M1M3).real();
            double tr24 = trace(M2M4).real();
            double tr14 = trace(M1M4).real();
            double tr23 = trace(M2M3).real();
            double tr1 = trace(mat[i1]).real();
            double tr2 = trace(mat[i2]).real();
            double tr3 = trace(mat[i3]).real();
            double tr4 = trace(mat[i4]).real();


            double res = dim*tr1234;
            res += eps[i1]*tr1*tr234;
            res += eps[i2]*tr2*tr134;
            res += eps[i3]*tr3*tr124;
            res += eps[i4]*tr4*tr123;
            res += eps[i1]*eps[i2]*tr12*tr34;
            res += eps[i1]*eps[i3]*tr13*tr24;
            res += eps[i1]*eps[i4]*tr14*tr23;

            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].real(); 

            return 2*cliff*res;
        }
        else
            return 0.;
    }
}

double Geom24::compute_A2(const int& i1, const int& i2) const
{
    // clifford product
    double cliff = omega_table_4[i2 + nHL*(i1 + nHL*(i2 + nHL*i1))].real();

    // base matrix products
    cx_mat M1M1 = mat[i1]*mat[i1];
    cx_mat M2M2 = mat[i2]*mat[i2];
    cx_mat M1M2 = mat[i1]*mat[i2];

    // traces
    double tr1122 = trace(M1M1*M2M2).real();
    double tr1212 = trace(M1M2*M1M2).real();
    double tr122 = trace(M1M2*mat[i2]).real();
    double tr112 = trace(M1M1*mat[i2]).real();
    double tr11 = trace(M1M1).real();
    double tr22 = trace(M2M2).real();
    double tr12 = trace(M1M2).real();
    double tr1 = trace(mat[i1]).real();
    double tr2 = trace(mat[i2]).real();
    
    
    if(cliff < 0)
    {
        // compute sum
        double res = dim*(2*tr1122 - tr1212);
        res += 2*eps[i1]*tr1*tr122;
        res += 2*eps[i2]*tr2*tr112;
        res += tr11*tr22;
        res += 2*eps[i1]*eps[i2]*tr12*tr12;

        return 2*dim_omega*res;
    }
    else
    {
        // compute sum
        double res = dim*(2*tr1122 + tr1212);
        res += 6*eps[i1]*tr1*tr122;
        res += 6*eps[i2]*tr2*tr112;
        res += 3*tr11*tr22;
        res += 6*eps[i1]*eps[i2]*tr12*tr12;

        return 2*dim_omega*res;
    }
}

double Geom24::compute_A(const int& i) const
{
    // base matrix products
    cx_mat M2 = mat[i]*mat[i];
    cx_mat M3 = M2*mat[i];

    // traces
    double tr1 = trace(mat[i]).real();
    double tr2 = trace(M2).real();
    double tr3 = trace(M3).real();
    double tr4 = trace(M3*mat[i]).real();

    double res = dim*tr4;
    res += 4*eps[i]*tr1*tr3;
    res += 3*tr2*tr2;

    return 2*dim_omega*res;
}

double Geom24::dirac4() const
{
    double res = 0.;

    // four distinct indices
    for(int i=0; i<nHL; ++i)
    {
        for(int j=i+1; j<nHL; ++j)
        {
            for(int k=j+1; k<nHL; ++k)
            {
                for(int l=k+1; l<nHL; ++l)
                    res += 8*(compute_A4(i,j,k,l)+compute_A4(i,j,l,k)+compute_A4(i,k,j,l));
            }
        }
    }

    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        for(int j=i+1; j<nHL; ++j)
            res += 2*compute_A2(i,j);
    }

    // all indices equal
    for(int i=0; i<nHL; ++i)
        res += compute_A(i);

    return res;
}


void Geom24::sample_mom(gsl_rng* engine)
{
    for(int i=0; i<nHL; ++i)
    {
        cx_mat temp(dim, dim);
        temp.imbue( [&engine](){ return cx_double(gsl_ran_gaussian(engine, 1.), gsl_ran_gaussian(engine, 1.)); } );

        mom[i] = (temp+temp.t())/2.;

        /*
        for(int j=0; j<dim; ++j)
        {
            for(int k=j+1; k<dim; ++k)
            {
                double x = gsl_ran_gaussian(engine, 1.)/sqrt(2.);
                double y = gsl_ran_gaussian(engine, 1.)/sqrt(2.);

                mom[i](j,k) = cx_double(x,y);
                mom[i](k,j) = cx_double(x,-y);
            }
        }

        for(int j=0; j<dim; ++j)
        {
            double x = gsl_ran_gaussian(engine, 1.);
            mom[i](j,j) = cx_double(x, 0.);
        }
        */
    }
}

double Geom24::calculate_K() const
{
    double res = 0;

    for(int i=0; i<nHL; ++i)
        res += trace(mom[i]*mom[i]).real();

    return res/2.;
}

double Geom24::calculate_H() const
{
    return calculate_S() + calculate_K();
}

void Geom24::leapfrog(const int& Nt, const double& dt)
{
    for(int i=0; i<nHL; ++i)
        mat[i] += (dt/2.)*mom[i].st();

    for(int j=0; j<Nt-1; ++j)
    {
        for(int i=0; i<nHL; ++i)
        {
            mom[i] += -dt*der_dirac24(i, true);
            mat[i] += dt*mom[i].st();
        }
    }
    
    for(int i=0; i<nHL; ++i)
    {
        mom[i] += -dt*der_dirac24(i, true);
        mat[i] += (dt/2.)*mom[i].st();
    }
}


double Geom24::HMC(const int& Nt, const double& dt, const int& iter, gsl_rng* engine, ostream& out_s, ostream& out_hl)
{
    double ar = 0;

    double Si, Ki, Hi;
    double Sf, Kf, Hf;

    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        sample_mom(engine);

        // store previous configuration
        cx_mat* mat_bk = new cx_mat [nHL];
        for(int j=0; j<nHL; j++)
            mat_bk[j] = mat[j];

        // set potential to previous final value,
        // unless it's the first iteration
        if(i)
            Si = Sf;
        else
            Si = calculate_S();

        Ki = calculate_K();
        Hi = Si+Ki;

        // leapfrog
        leapfrog(Nt, dt);

        // final hamiltonian
        Sf = calculate_S();
        Kf = calculate_K();
        Hf = Sf+Kf;

        // metropolis test
        if(Hf > Hi)
        {
            double r = gsl_rng_uniform(engine);
            double e = exp(Hi-Hf);

            if(r > e)
            {
                // restore old configuration
                for(int j=0; j<nHL; ++j)
                {
                    mat[j] = mat_bk[j];
                    Sf = Si;
                    Kf = Ki;
                    Hf = Hi;
                }
            }
            else ++ar;
        }
        else ++ar;
    
        // print S, K, H
        out_s << Sf << " " << Kf << " " << Hf << endl; 

        // print mat
        for(int j=0; j<nHL; ++j)
        {
            for(int k=0; k<dim; ++k)
            {
                for(int l=0; l<dim; ++l)
                    out_hl << mat[j](k,l).real() << " " << mat[j](k,l).imag() << " ";
            }
            out_hl << endl;
        }


        delete [] mat_bk;
    }

    return ar/iter;
}


double Geom24::dual_averaging(const int& Nt, double& dt, const int& iter, gsl_rng* engine, ostream& out_s, ostream& out_hl)
{
    // dual averaging variables
    double Stat = 0;
    double mu = log(10*dt);
    double shr = 0.05;
    int i0 = 10;
    double kappa = 0.75;
    double log_dt_avg = log(dt);


    double Si, Ki, Hi;
    double Sf, Kf, Hf;

    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        sample_mom(engine);

        // store previous configuration
        cx_mat* mat_bk = new cx_mat [nHL];
        for(int j=0; j<nHL; j++)
            mat_bk[j] = mat[j];

        // set potential to previous final value,
        // unless it's the first iteration
        if(i)
            Si = Sf;
        else
            Si = calculate_S();

        Ki = calculate_K();
        Hi = Si+Ki;

        // leapfrog
        leapfrog(Nt, dt);

        // final hamiltonian
        Sf = calculate_S();
        Kf = calculate_K();
        Hf = Sf+Kf;

        // metropolis test
        if(Hf > Hi)
        {
            double r = gsl_rng_uniform(engine);
            double e = exp(Hi-Hf);

            Stat += 0.65-e;

            if(r > e)
            {
                // restore old configuration
                for(int j=0; j<nHL; ++j)
                {
                    mat[j] = mat_bk[j];
                    Sf = Si;
                    Kf = Ki;
                    Hf = Hi;
                }
            }
        }
        else Stat += 0.65-1;

        
        // dual averaging
        double log_dt = mu - Stat*sqrt(i+1)/(shr*(i+1+i0));
        dt = exp(log_dt);
        double eta = pow(i+1, -kappa);
        log_dt_avg = eta*log_dt + (1-eta)*log_dt_avg;


    
        // print S, K, H
        out_s << Sf << " " << Kf << " " << Hf << endl; 

        // print mat
        for(int j=0; j<nHL; ++j)
        {
            for(int k=0; k<dim; ++k)
            {
                for(int l=0; l<dim; ++l)
                    out_hl << mat[j](k,l).real() << " " << mat[j](k,l).imag() << " ";
            }
            out_hl << endl;
        }


        delete [] mat_bk;
    }


    return exp(log_dt_avg);
}


    





cx_mat Geom24::compute_B4(const int& k, const int& i2, const int& i3, const int& i4, const double& cliff, const bool& neg) const
{
    if(neg)
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        double tr234 = trace(M2M3M4).imag();
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cx_double iu(0,1);
        cx_mat res(dim ,dim, fill::eye);
        res *= -2*eps[k]*tr234;
        res += double(dim)*iu*(M2M3M4 - M2M3M4.t());
        res += eps[i2]*tr2*iu*(M3M4 - M3M4.t());
        res += eps[i3]*tr3*iu*(M2M4 - M2M4.t());
        res += eps[i4]*tr4*iu*(M2M3 - M2M3.t());

        return cliff*res.st();
    }
    else
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        double tr234 = trace(M2M3M4).real();
        double tr23 = trace(M2M3).real();
        double tr24 = trace(M2M4).real();
        double tr34 = trace(M3M4).real();
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cx_mat res(dim ,dim, fill::eye);
        res *= 2*eps[k]*tr234;
        res += dim*(M2M3M4 + M2M3M4.t());
        res += eps[i2]*tr2*(M3M4 + M3M4.t());
        res += eps[i3]*tr3*(M2M4 + M2M4.t());
        res += eps[i4]*tr4*(M2M3 + M2M3.t());
        res += 2*eps[k]*eps[i2]*tr34*mat[i2];
        res += 2*eps[k]*eps[i3]*tr24*mat[i3];
        res += 2*eps[k]*eps[i4]*tr23*mat[i4];

        return cliff*res.st();
    }
}


cx_mat Geom24::compute_B2(const int& k, const int& i) const
{
    // clifford product
    double cliff = omega_table_4[i + nHL*(k + nHL*(i + nHL*k))].real();

    // base matrix products
    cx_mat MiMk = mat[i]*mat[k];
    cx_mat MiMi = mat[i]*mat[i];
    cx_mat MiMiMk = mat[i]*MiMk;
    cx_mat MiMkMi = MiMk*mat[i];

    // traces
    double triki = trace(MiMkMi).real();
    double trik = trace(MiMk).real();
    double trii = trace(MiMi).real();
    double tri = trace(mat[i]).real();
    double trk = trace(mat[k]).real();
    
    
    if(cliff < 0)
    {
        // compute sum
        cx_mat res(dim, dim, fill::eye);
        res *= eps[k]*triki;
        res += dim*(MiMiMk + MiMiMk.t() - MiMkMi);
        res += eps[i]*tri*(MiMk + MiMk.t());
        res += 2*eps[k]*eps[i]*trik*mat[i];
        res += eps[k]*trk*MiMi;
        res += trii*mat[k];

        return 2*dim_omega*res.st();
    }
    else
    {
        // compute sum
        cx_mat res(dim, dim, fill::eye);
        res *= 3*eps[k]*triki;
        res += dim*(MiMiMk + MiMiMk.t() + MiMkMi);
        res += 3*eps[i]*tri*(MiMk + MiMk.t());
        res += 6*eps[k]*eps[i]*trik*mat[i];
        res += 3*eps[k]*trk*MiMi;
        res += 3*trii*mat[k];

        return 2*dim_omega*res.st();
    }
}

cx_mat Geom24::compute_B(const int& k) const
{
    // base matrix products
    cx_mat M2 = mat[k]*mat[k];
    cx_mat M3 = mat[k]*M2;

    // traces
    double tr3 = trace(M3).real();
    double tr2 = trace(M2).real();
    double tr1 = trace(mat[k]).real();

    cx_mat res(dim, dim, fill::eye);
    res *= eps[k]*tr3;
    res += dim*M3;
    res += 3*tr2*mat[k];
    res += 3*eps[k]*tr1*M2;

    return 2*dim_omega*res.st();
}

/*
cx_mat Geom24::der_dirac4_new(const int& k, const bool& herm) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        for(int i2=i1+1; i2<nHL; ++i2)
        {
            for(int i3=i2+1; i3<nHL; ++i3)
            {
                for(int i4=i3+1; i4<nHL; ++i4)
                {
                    if(i1==k)
                    {
                        // epsilon factor
                        int e = eps[i1]*eps[i2]*eps[i3]*eps[i4];

                        if(e<0)
                        {
                            // clifford product
                            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].imag(); 

                            if(cliff != 0.)
                            {
                                cx_mat temp = compute_B4(i1,i2,i3,i4, cliff, true) + compute_B4(i1,i2,i4,i3, cliff, true) + compute_B4(i1,i3,i2,i4, cliff, true);
                                res += temp + temp.t();
                            }
                        }
                        else
                        {
                            // clifford product
                            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].real(); 

                            if(cliff != 0.)
                            {
                                cx_mat temp = compute_B4(i1,i2,i3,i4, cliff, false) + compute_B4(i1,i2,i4,i3, cliff, false) + compute_B4(i1,i3,i2,i4, cliff, false);
                                res += temp + temp.t();
                            }
                        }
                    }
                    else if(i2==k)
                    {
                        // epsilon factor
                        int e = eps[i1]*eps[i2]*eps[i3]*eps[i4];

                        if(e<0)
                        {
                            // clifford product
                            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].imag(); 

                            if(cliff != 0.)
                            {
                                cx_mat temp = compute_B4(i2,i3,i4,i1, cliff, true) + compute_B4(i2,i4,i3,i1, cliff, true) + compute_B4(i2,i4,i1,i3, cliff, true);
                                res += temp + temp.t();
                            }
                        }
                        else
                        {
                            // clifford product
                            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].real(); 

                            if(cliff != 0.)
                            {
                                cx_mat temp = compute_B4(i2,i3,i4,i1, cliff, false) + compute_B4(i2,i4,i3,i1, cliff, false) + compute_B4(i2,i4,i1,i3, cliff, false);
                                res += temp + temp.t();
                            }
                        }
                    }
                    else if(i3==k)
                    {
                        // epsilon factor
                        int e = eps[i1]*eps[i2]*eps[i3]*eps[i4];

                        if(e<0)
                        {
                            // clifford product
                            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].imag(); 

                            if(cliff != 0.)
                            {
                                cx_mat temp = compute_B4(i2,i3,i4,i1, cliff, true) + compute_B4(i2,i4,i3,i1, cliff, true) + compute_B4(i2,i4,i1,i3, cliff, true);
                                res += temp + temp.t();
                            }
                        }
                        else
                        {
                            // clifford product
                            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].real(); 

                            if(cliff != 0.)
                            {
                                cx_mat temp = compute_B4(i2,i3,i4,i1, cliff, false) + compute_B4(i2,i4,i3,i1, cliff, false) + compute_B4(i2,i4,i1,i3, cliff, false);
                                res += temp + temp.t();
                            }
                        }
                    }
                    else if(i4==k)
                    {}
                        }
                    }
                }
            }
        }
    }

    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        if(i != k)
            res += compute_B2(k,i);
    }

    // all indices equal
    res += compute_B(k);


    if(herm)
        return 2*(res+res.t());
    else
        return 4*res;

}
*/


cx_mat Geom24::der_dirac4(const int& k, const bool& herm) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        if(i1 != k)
        {
            for(int i2=i1+1; i2<nHL; ++i2)
            {
                if(i2 != k)
                {
                    for(int i3=i2+1; i3<nHL; ++i3)
                    {
                        if(i3 != k)
                        {
                            // epsilon factor
                            int e = eps[k]*eps[i1]*eps[i2]*eps[i3];

                            if(e<0)
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].imag(); 

                                if(cliff1 != 0.)
                                {
                                    cx_mat temp = compute_B4(k,i1,i2,i3, cliff1, true) + compute_B4(k,i1,i3,i2, cliff2, true) + compute_B4(k,i2,i1,i3, cliff3, true);
                                    res += temp + temp.t();
                                }
                            }
                            else
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].real(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].real(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].real(); 

                                if(cliff1 != 0.)
                                {
                                    cx_mat temp = compute_B4(k,i1,i2,i3, cliff1, false) + compute_B4(k,i1,i3,i2, cliff2, false) + compute_B4(k,i2,i1,i3, cliff3, false);
                                    res += temp + temp.t();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        if(i != k)
            res += compute_B2(k,i);
    }

    // all indices equal
    res += compute_B(k);


    if(herm)
        return 2*(res+res.t());
    else
        return 4*res;

}



cx_mat Geom24::der_dirac2(const int& k) const
{
    cx_mat res(dim, dim, fill::eye);

    res *= eps[k]*trace(mat[k]).real();
    res += dim*mat[k].st();

    return 4*dim_omega*res;
}

cx_mat Geom24::der_dirac24(const int& k, const bool& herm) const
{
    return g2*der_dirac2(k) + der_dirac4(k, herm);
}


// old function for trD4

double Geom24::dirac4_old() const
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
                        
                        res += dim_omega*temp;
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
                        
                        res += dim_omega*temp;
                    }
                    

                    else
                    {
                        cx_double cliff = omega_table_4[i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                        
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
                            cx_double T1 = trM0M1M2M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM0M1M2M3);
                            res += 2.*dim*(cliff*T1).real();

                            // tr3tr1 terms

                            cx_double T3 = (double)(eps[i[3]])*trM0M1M2 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]])*conj(trM0M1M2);
                            T3 = trM3*T3;

                            cx_double T4 = (double)(eps[i[2]])*trM0M1M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[3]])*conj(trM0M1M3);
                            T3 += trM2*T4;

                            cx_double T5 = (double)(eps[i[1]])*trM0M2M3 + (double)(eps[i[0]]*eps[i[2]]*eps[i[3]])*conj(trM0M2M3);
                            T3 += trM1*T5;

                            cx_double T6 = (double)(eps[i[0]])*trM1M2M3 + (double)(eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM1M2M3);
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
                    cx_double cliff = omega_table_4[i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                    
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
                        cx_double T1 = trM0M1M2M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM0M1M2M3);
                        res += dim*2.*(cliff*T1).real();

                        // tr3tr1 terms
                        cx_double T3 = (double)(eps[i[3]])*trM0M1M2 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]])*conj(trM0M1M2);
                        T3 = trM3*T3;

                        cx_double T4 = (double)(eps[i[2]])*trM0M1M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[3]])*conj(trM0M1M3);
                        T3 += trM2*T4;

                        cx_double T5 = (double)(eps[i[1]])*trM0M2M3 + (double)(eps[i[0]]*eps[i[2]]*eps[i[3]])*conj(trM0M2M3);
                        T3 += trM1*T5;

                        cx_double T6 = (double)(eps[i[0]])*trM1M2M3 + (double)(eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM1M2M3);
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




double Geom24::calculate_S_old() const
{
    return g2*dirac2() + dirac4_old();
}
