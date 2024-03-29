#include <iostream>
#include <vector>
#include <cmath>
#include "ReferenceElement.h"
using namespace std;

vector<double> vecsum(vector <double> &A, vector<double> &B, int dimAB1, int dimAB2);

vector<double> vecdiff(vector <double> &A, vector<double> &B, int dimAB1, int dimAB2);

vector<double> scalarmult(double s, vector<double> &A, int dimA1, int dimA2);

vector<double> matmul(vector<double>& A, vector<double>& B, int dimA1,
		       int dimA2B1, int dimB2);

void rightHandSide(vector<double>& thegrid, vector<double>& uh,
		   vector<double>& rhs, int order, int numelem, int pdenum,
		   int diffnum, double speed);

int main()
{
  //  double lowlim =-20.;
  //  double uplim = 45.5451774444796;
  double lowlim = 0.;
  double uplim = 20.;
  int numelem = 20;
  int order = 16;
  int pdenum = 3;
  int diffnum = 2;
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

  rightHandSide(thegrid,uh,rhs, order, numelem, pdenum, diffnum, speed);

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
		   int diffnum, double speed){

  vector<double> Atrimmed(diffnum*diffnum,0.0);
  Atrimmed[1]=-speed*speed;
  Atrimmed[2]=-1.0;
  
  vector<double> lamb(diffnum*diffnum);
  lamb[0]=-speed;
  lamb[1]=speed;
  lamb[2]=-speed;
  lamb[3]=speed;

  vector<double> smatrix(diffnum*diffnum);
  smatrix[0]=-sqrt(1./2.)*speed;
  smatrix[1]=-sqrt(1./2.)*speed;
  smatrix[2]=-sqrt(1./2.);
  smatrix[3]=sqrt(1./2.);

  //CHECK THIS WITH MATHEMATICA WHEN I START USING SPEED NOT EQUAL TO ONE
  vector<double> sinv(diffnum*diffnum);
  sinv[0]=-sqrt(1./2.);
  sinv[1]=-sqrt(1./2.)/speed;
  sinv[2]=-sqrt(1./2.);
  sinv[3]=sqrt(1./2.)/speed;

  int pdebase = numelem*(order+1);
  int indL = 0;
  int indR = order;
  vector<double> nx(2);
  nx[0]=-1.0;
  nx[1]=1.0;


  for(int elemnum=0; elemnum<numelem; elemnum++){
    vector<double> nfluxL(diffnum);
    vector<double> nfluxR(diffnum);
    

    int elemcount = elemnum*(order+1);
    
    vector<double> uintL(diffnum);
    vector<double> uintR(diffnum);
    vector<double> uextL(diffnum);
    vector<double> uextR(diffnum); 

    uintL[0]=uh[pdebase+elemcount+indL];
    uintR[0]=uh[pdebase+elemcount+indR];
    uintL[1]=uh[2*pdebase+elemcount+indL];
    uintR[1]=uh[2*pdebase+elemcount+indR];


    int elemcount2;
    if(elemnum>0){
      elemcount2 =  (elemnum-1)*(order+1);
    } else {
      elemcount2=numelem*(order);
    }
    uextL[0]=uh[pdebase+elemcount2+indR];
    uextL[1]=uh[2*pdebase+elemcount2+indR];
    
    int elemcount3;
    if(elemnum<numelem-1){
      elemcount3 = (elemnum+1)*(order+1);
    }else{
      elemcount3 =0;
    }
    uextR[0]=uh[pdebase+elemcount3+indL];
    uextR[1]=uh[2*pdebase+elemcount3+indL];

    
    for(int j=0; j<2; j++){
      vector<double> lambplus(2*diffnum);
      vector<double> lambminus(2*diffnum);
      for(int k=0; k<diffnum; k++){
	if(nx[j] * lamb[2 * j + k] < 0.0) {
	  lambminus[2 * k + k] = nx[j] * lamb[2 * j + k];
	}else{
	  lambplus[2 * k + k] = nx[j] * lamb[2 * j + k];
	}
      }
      if(j==0){
	vector<double> tempL =  matmul(sinv, uintL,2,2,1);
	vector<double> nfluxL1 = matmul(lambplus, tempL,2,2,1);
	tempL = matmul(sinv, uextL,2,2,1);
	vector<double> nfluxL2 = matmul(lambminus, tempL,2,2,1);
	
	vector<double> nfluxL3 = vecsum(nfluxL1, nfluxL2,2,1);
	nfluxL = matmul(smatrix,nfluxL3,2,2,1); 


      } else {
	vector<double> tempR = matmul(sinv, uintR,2,2,1);
	vector<double> nfluxR1 = matmul(lambplus, tempR,2,2,1);
	tempR = matmul(sinv,uextR,2,2,1);
	vector<double> nfluxR2 = matmul(lambminus, tempR,2,2,1);
	vector<double> nfluxR3 = vecsum(nfluxR1,nfluxR2, 2,1);
	nfluxR = matmul(smatrix,nfluxR3,2,2,1); 
      }
    }

    //    cout << nfluxL[0] << " "<< nfluxL[1] <<" " << nfluxR[0] << " " <<nfluxR[1] << endl;
    vector<double> duL(diffnum);
    vector<double> duR(diffnum);

    vector<double> temp = matmul(Atrimmed,uintL,2,2,1);
    vector<double> temp2 = scalarmult(nx[0], temp, 2, 1);
    duL=vecdiff(temp2,nfluxL,2,1);
    temp = matmul(Atrimmed,uintR,2,2,1);
    temp2 = scalarmult(nx[1], temp, 2, 1);
    duR = vecdiff(temp2,nfluxR,2,1); 

    cout << duL[0] << " " << duL[1] << " " << duR[0] << " " << duR[1] << endl;
    //du's should be zero after first call to rhs-- checks out to within
    // roundoff
  }
  /*  vector<double> ident = matmul(smatrix,sinv,2,2,2);
  for(int i=0; i<4; i++){
    cout << ident[i] << " ";
  }
  cout << endl;
  vector<double> temp10 = matmul(lamb,sinv,2,2,2);
  vector<double> AT = matmul(smatrix,temp10,2,2,2);
  for (int i=0; i<4; i++){
      cout << AT[i] << " " ;
    }
  cout << endl;
  */
}

vector<double> vecsum(vector <double> &A, vector<double> &B, int dimAB1, int dimAB2){
  vector<double> C(dimAB1*dimAB2);
  
  for(int i=0; i<dimAB1; i++){
    for(int j=0; j<dimAB2; j++){
      C[j*dimAB1+i]=A[j*dimAB1+i]+B[j*dimAB1+i];
    }
  }
  return C;
}

vector<double> vecdiff(vector <double> &A, vector<double> &B, int dimAB1, int dimAB2){
  vector<double> C(dimAB1*dimAB2);
  
  for(int i=0; i<dimAB1; i++){
    for(int j=0; j<dimAB2; j++){
      C[j*dimAB1+i]=A[j*dimAB1+i]-B[j*dimAB1+i];
    }
  }
  return C;
}

vector<double> scalarmult(double s, vector<double> &A, int dimA1, int dimA2){
  vector<double> C(dimA1*dimA2);
  
  for(int i=0; i<dimA1; i++){
    for(int j=0; j<dimA2; j++){
      C[j*dimA1+i]=A[j*dimA1+i]*s;
    }
  }
  return C;
}

vector<double> matmul(vector<double> &A, vector<double> &B, int dimA1,
		       int dimA2B1, int dimB2){
  vector<double> C(dimA1*dimB2);
  
  for(int i=0; i<dimA1; i++){
    for(int j=0; j<dimB2; j++){
      double sum = 0.0;
      for(int k=0; k<dimA2B1;k++){
	sum+=A[i*dimA2B1+k]*B[k*dimB2+j];
      }
      C[i*dimB2+j]=sum;
    }
  }
  return C;

}


  
