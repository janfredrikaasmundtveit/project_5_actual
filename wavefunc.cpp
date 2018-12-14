#include "wavefunc.h"
#include "vectormatrixclass.h"

wavefuncgaussian::wavefuncgaussian(){
alpha=0;omega=0;	
}

wavefuncgaussian::wavefuncgaussian(double a,double o){
alpha=a; omega=o; 
}

double wavefuncgaussian::Hamiltonian(Vector &R1, Vector &R2){
	double r1=R1.Norm_l2(); 
	double r2=R2.Norm_l2(); //normsquared
	double E=0.5*omega*omega*(r1+r2)*(1-alpha*alpha)+3*alpha*omega;
	return E;
}
double wavefuncgaussian::Hamiltonianint(Vector &R1, Vector &R2){
	Vector R12=R1-R2;
	double r12=R12.Norm_l2();
	double r1=R1.Norm_l2(); 
	double r2=R2.Norm_l2(); //normsquared
	double E=0.5*omega*omega*(r1+r2)*(1-alpha*alpha)+3*alpha*omega+1/sqrt(r12);
	return E;
}

wavefuncgaussian::~wavefuncgaussian(){};

double wavefuncgaussian::evaluate(Vector& r1,Vector& r2){
	double R1=r1.Norm_l2(); 
	double R2=r2.Norm_l2();
	double a=exp(-alpha*omega*(R1+R2)); // unnormalized psi^2
	return a;
}

wavefuncexp::wavefuncexp(){
alpha=0;omega=0;	beta=0;
}

wavefuncexp::wavefuncexp(double a,double o,double b){
alpha=a; omega=o; beta=b; // andvanced stuff
}


double wavefuncexp::Hamiltonianint(Vector &R1, Vector &R2){
	Vector R12=R1-R2;
	double r12=sqrt(R12.Norm_l2());
	double r1=R1.Norm_l2(); 
	double r2=R2.Norm_l2(); //normsquared
	double E=0.5*omega*omega*(r1+r2)*(1-alpha*alpha)+3*alpha*omega+1/r12+0.5/((1+beta*r12)*(1+beta*r12))*(alpha*omega*r12-0.5/((1+beta*r12)*(1+beta*r12))-2/r12+2*beta/(1+beta*r12));
	return E;
}


double wavefuncexp::Kinetic(Vector &R1, Vector &R2){
	Vector R12=R1-R2;
	double r12=sqrt(R12.Norm_l2());
	double r1=R1.Norm_l2(); 
	double r2=R2.Norm_l2(); //normsquared
	double E=0.5*omega*omega*(r1+r2)*(1-alpha*alpha);
	return E;
}
double wavefuncexp::potential(Vector &R1, Vector &R2){
	Vector R12=R1-R2;
	double r12=sqrt(R12.Norm_l2());
	double r1=R1.Norm_l2(); 
	double r2=R2.Norm_l2(); //normsquared
	double E=3*alpha*omega+1/r12+0.5/((1+beta*r12)*(1+beta*r12))*(alpha*omega*r12-0.5/((1+beta*r12)*(1+beta*r12))-2/r12+2*beta/(1+beta*r12));
	return E;
}

	
	wavefuncexp::~wavefuncexp(){};

double wavefuncexp::evaluate(Vector& r1,Vector& r2){
	Vector R12=r1-r2;
	double R1=r1.Norm_l2(); 
	double R2=r2.Norm_l2();
	double r12=sqrt(R12.Norm_l2());
	double a=exp(-alpha*omega*(R1+R2))*exp(2*r12/(1+beta*r12)); // unnormalized
	return a;
}