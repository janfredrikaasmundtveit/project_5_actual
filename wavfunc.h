#ifndef _wavefunc
#define _wavefunc


#include <cmath>
#include <iostream>
using namespace  std;

class wavefuncgassian;
class wavefuncexp;

//twoparticle wavefunctions of the form Ce^{-\alpha\omega(r_1^2+r_2^2)/2}
class wavfuncgaussian
{
 private:
 	double alpha;
 	double omega;
public:
	wavfuncgaussian();
	wavfuncgaussian(double alpha,double omega);
	~wavfuncgaussian();
	double Hamiltonian(const Vector r1, const Vector r2); //change to vector.
	double Hamiltonianint(const Vector r1, const Vector r2);
};

//twoparticle wave fvunction of the form  Ce^{-\alpha\omega(r_1^2+r_2^2)/2}e^{\frac{r_{12}}{2(1+\beta)}}
class wavfuncexp
{
 private:
 	double alpha;
 	double omega;
 	double beta;
public:
	wavfuncexp();
	wavfuncexp(double alpha,double beta,double omega);
	~wavfuncexp();
	double Hamiltonian(const Vector r1, const Vector r2);
	double Hamiltonianint(const Vector r1, const Vector r2);
};


#endif