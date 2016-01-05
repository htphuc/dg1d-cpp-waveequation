#include <iostream>
#include <vector>
#include <cmath>
#include "ReferenceElement.h"
using namespace std;

int main()
{
  //  double lowlim =-20.;
  //  double uplim = 45.5451774444796;
  double lowlim = 0.;
  double uplim = 20.;
  int numelem = 20;
  int order = 16;
  int pdenum = 3;
  double sineamp = 1.0;
  double sinewavelength = 20.0;
  double sinephase = 0.;
  double omega = 2.*M_PI/sinewavelength;
  //  double gaussamp = 1.0;
  //  double gausssigma = 1.0;
  //  double gaussmu = 10.1;
  double speed = 1.0;
  
  
  vector<double> thegrid(numelem*(order+1),0.0);
  ReferenceElement refelem(order);
  for(int elemnum = 0; elemnum < numelem; elemnum ++){
    double jacobian = (uplim-lowlim)/double(numelem)/2.;
    double lower = lowlim+jacobian*(elemnum+0.5)*2.;
    for(int nodenum = 0; nodenum <= order; nodenum ++){
      thegrid[elemnum*(order+1)+nodenum] =
	lower+jacobian*refelem.getNode(nodenum);
    }
  }
  int elemcount = numelem*(order+1);
  vector<double> uh(pdenum*numelem*(order+1),0.0);
  vector<double> rhs(pdenum*numelem*(order+1),0.0);
  for(int elemnum = 0; elemnum< numelem; elemnum++){
    for(int nodenum = 0; nodenum <= order; nodenum++){
      int nodecount = elemnum*(order+1)+nodenum;
      //      double gaussian = gaussamp * exp(-pow((thegrid[nodecount]-gaussmu),2.0)
      //				       /2.0/pow(gausssigma,2.0));
      //      double dgauss = -(thegrid[nodecount]-gaussmu)/pow(gausssigma,2.0)
      //	*gaussian;
      //      uh[nodecount] = gaussian;
      //      uh[2*elemcount+nodecount] = dgauss;
      double psivar = sineamp * sin(omega* thegrid[nodecount]+sinephase);
      double rhovar = omega*sineamp*cos(omega*thegrid[nodecount]+sinephase);
      double pivar = - speed * rhovar;
      uh[nodecount] = psivar;
      uh[elemcount + nodecount] = pivar;
      uh[2*elemcount + nodecount] = rhovar;
      rhs[nodecount] = 0.0;
      rhs[elemcount + nodecount] = 0.0;
      rhs[2*elemcount + nodecount] = 0.0;
      
    }
  }

  rightHandSide(thegrid,uh,rhs, order, numelem, pdenum, speed);

  //  int i=0;
  //  while()
  //   {
  //    i = i+1;
  //    rk4(thegrid,uh,rhs);
      
  
  //output initial data (redirect to init.dat)
  /*  cout << setprecision(15);
  for(int elemnum = 0; elemnum< numelem; elemnum++){
    for(int nodenum = 0; nodenum <= order; nodenum++){
      int nodecount = elemnum*(order+1)+nodenum;
      cout <<  thegrid[nodecount];
      for(int eqnnum =0; eqnnum < pdenum; eqnnum ++){
	
	cout <<  " " <<uh[eqnnum*elemcount+nodecount]; 
      }
      cout << endl;
      }
      }*/

  
  
  
}  

void rightHandSide(vector<double>& thegrid, vector<double>& uh,
		   vector<double>& rhs, int order, int numelem, int pdenum,
		   double speed){

  //debug node and elemcount section NEXT TIME. WORKING on Char Flux. Should
  // be zero in first sub time step.and zeroth time step. HERE.
  vector<double> lambda(4,0.0);
  lambda[0]=-speed;
  lambda[3]=speed;

  vector<double> smatrix(4);
  smatrix[0]=speed;
  smatrix[1]=-speed;
  smatrix[2]=1.0;
  smatrix[3]=1.0;

  vector<double> sinv(4);
  sinv[0]=0.5/speed;
  sinv[1]=0.5;
  sinv[2]=-0.5/speed;
  sinv[3]=0.5;

  int pdebase = numelem*(order+1);
  int indL = order;
  int indR = 0;
  for(int elemnum=0; elemnum<numelem; elemnum++){

    int elemcount = elemnum*(order+1);
    
    vector<double> uint((pdenum-1)*2);
    vector<double> uext((pdenum-1)*2); // times 2 for left and right hand side

    uint[0]=uh[pdebase+elemcount+indL];
    uint[1]=uh[pdebase+elemcount+indR];
    uint[2]=uh[2*pdebase+elemcount+indL];
    uint[3]=uh[2*pdebase+elemcount+indR];
    
    
    if(elemnum<=0){
      elemcount = (numelem-1)*(order+1);
    }
    if(elemnum>=numelem-1){
      elemcount= 0;
    }
  
    uext[0]=uh[pdebase+elemcount+indL];
    uext[1]=uh[pdebase+elemcount+indR];
    uext[2]=uh[2*pdebase+elemcount+indL];
    uext[3]=uh[2*pdebase+elemcount+indR];

    vector<double> lambdaplus(4,0.0);
