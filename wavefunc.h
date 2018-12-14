#ifndef _wavefunc
#define _wavefunc


#include <cmath>
#include <iostream>
using namespace  std;

class wavefuncgassian;
class wavefuncexp;
class Vector;

//twoparticle wavefunctions of the form Ce^{-\alpha\omega(r_1^2+r_2^2)/2}
class wavefuncgaussian
{
 
public:
	double alpha;
 	double omega;
	wavefuncgaussian();
	wavefuncgaussian(double a,double o);
	~wavefuncgaussian();
	double evaluate(Vector& r1,Vector& r2);
	double Hamiltonian(Vector& r1, Vector& r2); 
	double Hamiltonianint(Vector& r1, Vector& r2);

};

//twoparticle wave fvunction of the form  Ce^{-\alpha\omega(r_1^2+r_2^2)/2}e^{\frac{r_{12}}{2(1+\beta)}}
class wavefuncexp
{

 
public:
	double alpha;
 	double omega;
 	double beta;
	wavefuncexp();
	wavefuncexp(double a,double b,double o);
	~wavefuncexp();
	double evaluate(Vector& r1,Vector& r2);
	double Hamiltonianint(Vector& r1, Vector& r2);
		double Kinetic(Vector& r1, Vector& r2);
	double potential(Vector& r1, Vector& r2);
};



#endif