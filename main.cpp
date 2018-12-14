#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
//#include  "mpi.h"
#include <ctime>
#include <fstream>
#include <string>
#include <time.h>
#include <random>
#include<armadillo>
#include"wavefunc.h"
#include "vectormatrixclass.h"

using namespace arma;
using namespace std;
ofstream ofile;
void output(double,int,double,double,int);

int main (int argc, char* argv[])
{ string filename;
int MCC,MC;
if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file,cycles" << endl;
    exit(1);
  }
  else{
    filename=argv[1];
    MC=atoi(argv[2]);
    MCC=pow(10,MC);
   }
   ofile.open(filename);
//Vector r1,r2;
//wavefuncgaussian psi;


double omega=0.7; 
double beta=0.3;
double alpha=0.99;

Vector r1(3); Vector r2(3); Vector nr1(3); Vector nr2(3);
wavefuncexp psi(alpha,omega,beta); 
	double E=psi.Hamiltonianint(r1,r2); 
double T=psi.Kinetic(r1,r2);
 double V=psi.potential(r1,r2);
double psi2=psi.evaluate(r1,r2);  //cout << E << endl;
 double d=1.4*exp(-0.7*alpha); //tweek for different wavefuncs or omegas.
double totE=0.0; double totE2=0.0; int acc=0; double totT=0.0; double totV=0.0;
double xr=0.0;
for(int i=0; i < MCC; i++){
	std::random_device rd;   std::mt19937_64 gen(rd());
   // Set up the uniformdistribution for x \in [[0, 1]  
    std::uniform_real_distribution<double>
	RandomNumberGenerator(0.0,1.0); //mapping[0,1] to [-1,1] 
  	double x1=2*RandomNumberGenerator(gen)-1; double x2=2*RandomNumberGenerator(gen)-1;
	double y1=2*RandomNumberGenerator(gen)-1;  double y2=2*RandomNumberGenerator(gen)-1; 
	double z1=2*RandomNumberGenerator(gen)-1;
	double z2=2*RandomNumberGenerator(gen)-1; 

	 double n1[3]={d*x1,d*y1,d*z1};
	Vector d1(3,n1);  
	
	 double n2[3]={d*x2,d*y2,d*z2}; Vector d2(3,n2);
	nr1=r1+d1;  nr2=r2+d2;
 double Tn=psi.Kinetic(nr1,nr2);
 double Vn=psi.potential(nr1,nr2);
  double En=psi.Hamiltonianint(nr1,nr2);
  double psi2n=psi.evaluate(nr1,nr2);   

     if( RandomNumberGenerator(gen) <=psi2n/psi2){   
           r1=nr1; r2=nr2;
           //r12=r1-r2;
            //xr=sqrt(r12.Norm_l2())*psi2n;
            E=En; psi2=psi2n;
			totE+=E;  totT+=Tn; totV+=Vn;
		    totE2+=E*E; acc++;    
 }

}

output(omega,MCC,totT,totV,acc);


   ofile.close();
   return 0;
}


void output(double alpha,int MCC,double totE,double totE2,int acc){
double E=totE/((double)(acc));
double E2=totE2/((double)(acc));
double s=E2-E*E;
double a=(double(acc))/(double(MCC));
double R=E/E2;
ofile << setw(15) << setprecision(8) << alpha;
ofile << setw(15) << setprecision(8) << R;
ofile << setw(15) << setprecision(8) << E2; // Potential in viral theorm part.
//ofile << setw(15) << setprecision(8) << s;
ofile << setw(15) << setprecision(8) << a << "\n";

}