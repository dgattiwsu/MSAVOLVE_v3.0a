#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include "matrix.c"
#define TINY 1e-20

using namespace std;

static const double log_10=log(10.);

double alpha_exponent=0.;
double K=5;//precision
int A=20;//amino acid alphabet size


// function to convert amino acids to integers
inline int convert(char a){
  switch (a){
  case 'A':  return 0;  break;
  case 'a':  return 0;  break;
  case 'C':  return 1;  break;
  case 'c':  return 1;  break;
  case 'D':  return 2;  break;
  case 'd':  return 2;  break;
  case 'E':  return 3;  break;
  case 'e':  return 3;  break;
  case 'F':  return 4;  break;
  case 'f':  return 4;  break;
  case 'G':  return 5;  break;
  case 'g':  return 5;  break;
  case 'H':  return 6;  break;
  case 'h':  return 6;  break;
  case 'I':  return 7;  break;
  case 'i':  return 7;  break;
  case 'K':  return 8;  break;
  case 'k':  return 8;  break;
  case 'L':  return 9;  break;
  case 'l':  return 9;  break;
  case 'M':  return 10; break;
  case 'm':  return 10; break;
  case 'N':  return 11; break;
  case 'n':  return 11; break;
  case 'P':  return 12; break;
  case 'p':  return 12; break;
  case 'Q':  return 13; break;
  case 'q':  return 13; break;
  case 'R':  return 14; break;
  case 'r':  return 14; break;
  case 'S':  return 15; break;
  case 's':  return 15; break;
  case 'T':  return 16; break;
  case 't':  return 16; break;
  case 'V':  return 17; break;
  case 'v':  return 17; break;
  case 'W':  return 18; break;
  case 'w':  return 18; break;
  case 'Y':  return 19; break;
  case 'y':  return 19; break;
  case '-':  return 20; break; 
  case '.':  return 20; break; 
    //very rare cases
  case 'Z':  return 0; break; // IUPAC degenerate code
  case 'z':  return 0; break;
  case 'X':  return 0; break; // any amino acid, IUPAC
  case 'x':  return 0; break;
  case 'B':  return 0; break; // IUPAC degenerate code
  case 'b':  return 0; break; 
  
  default: cout<< "error, amino acid "<<a<<" not recognized" <<endl;
    exit(1);
    return -1;
  }
}





inline double transf(double log_x, double min_logW, double max_logW){
  //values must be in log-base 10!

  double a=1.;
  if(max_logW-min_logW>K){
    a=K/(max_logW-min_logW);
  }
  double x_transf=pow(10., -K)*pow(10., a*(log_x-min_logW));
  return x_transf;//not in log_scale!
}


inline double determinant(Mat<double>& a, int n){

  //calculates the determinant in logscale, without the sign
  //sign will be printed to STDERR
  //by setting n, one can specify whether one wants to calculate
  //the determinant of the whole matrix or a submatrix a[0....n][0....n]
  if(a.nrows()!=a.ncols()){
    cerr<<"error in determinant: matrix a must be square\n";
    exit(1);
  }
  //  Mat<double> tmp_mat=a;
  vector<int> indx(n);
  int i,imax,j,k;
  double big,dum,sum,temp;
  vector<double> vv(n);
  double d;
  double det=0;
  d=1.0;
  for(i=0;i<n;i++){
    big=0.0;
    for(j=0;j<n;j++)
      if((temp=fabs(a[i][j]))>big) big=temp;
    if(big==0.0){a.print_mat(); cerr<<"Singular matrix in routine ludcmp"<<endl; exit(1);}
    vv[i]=1.0/big;
  }
  for(j=0;j<n;j++){
    for(i=0;i<j;i++){
      sum=a[i][j];
      for(k=0;k<i;k++) sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<n;i++){
      sum=a[i][j];
      for(k=0;k<j;k++)
	sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum))>=big){
	big=dum;
	imax=i;
      }
    }
    if(j!=imax){
      for(k=0;k<n;k++){
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      d=-d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if(a[j][j]==0.0){
      a[j][j]=TINY;
      cerr<<"watch out, matrix is singular\n";
    }
    if(j!=(n-1)){
      dum=1.0/(a[j][j]);
      for(i=j+1;i<n;i++) a[i][j]*=dum;
    }
  }
  
  double signum=d;//the sign of the determinant

  for(int j=0;j<n;j++){
    if(a[j][j]<0){
      signum*=-1;
    }
    det+=log(fabs(a[j][j]));
  }
  return det;
}


void rescale_logR(Mat<double>& M){
  //first get min and max
  double max=M[0][1];
  double min=M[0][1];
  for(int i=0;i<M.nrows();i++){
    for(int j=0;j<i;j++){
      if(M[i][j]>max){max=M[i][j];}
      if(M[i][j]<min){min=M[i][j];}
    }
  }
  max=max/log_10;
  min=min/log_10;

  if(max-min>K){
    alpha_exponent=K/(max-min);
  }
  else{
    alpha_exponent=1.;
  }


  for(int i=0;i<M.nrows();i++){
    for(int j=0;j<i;j++){
      M[i][j]=transf(M[i][j]/log_10, min, max);
      M[j][i]=M[i][j];
    }
  }
}

inline double score(Mat<double>& M){

  double t_score=0.;
  Mat<double> current_Q=M;

  for(int j=0;j<current_Q.ncols();j++){
    double sum=0.;
    for(int i=0;i<current_Q.nrows();i++){
      current_Q[i][j]=-current_Q[i][j];
      if(i!=j){
	sum+=fabs(current_Q[i][j]);
      }
    }
    current_Q[j][j]=sum;
  }
  t_score=determinant(current_Q, current_Q.nrows()-1);
  return t_score;
}





Mat<double> calculate_logR(Mat<char>& alignment){

 
  Mat<double> M(0., alignment.ncols(), alignment.ncols());
  vector<double> correction(alignment.ncols(), 0.);  
  double min;
  for (int l=0;l<M.ncols();l++){//Hdim2
    for (int k=0; k<l;k++){//Rdim2
      
      double logP_dep=0;
      double logP_indep=0;
      vector<double> n_i(A, 0.);
      vector<double> n_j(A, 0.);
      double n=0;   // sum over i and j
      const double alpha=0.5;
      Mat<double> pair_freq(0., A, A);
      //initialising the variables
 
      for (int i=0;i<alignment.nrows();i++){
	if((alignment[i][l] != '-') && (alignment[i][k] != '-')){
	  pair_freq[convert(alignment[i][l])][convert(alignment[i][k])]++; 
	}
      }
      
      for (int i=0; i<A; i++){
	for (int j=0; j<A;j++){
      	  n_i[i]+=pair_freq[i][j];
	  n_j[j]+=pair_freq[i][j];
	  n+=pair_freq[i][j];}
      }
      
      logP_dep=lgamma(alpha*pow(A,2.))-lgamma(n+alpha*pow(A,2.));
      logP_indep=2*lgamma(alpha*pow(A,2.))-2*lgamma(n+alpha*pow(A,2.));
      for (int i=0;i<A;i++){
	logP_indep+=lgamma(n_i[i]+A*alpha)-2*lgamma(A*alpha)+lgamma(n_j[i]+A*alpha);
	for (int j=0;j<A;j++){
	  logP_dep+=lgamma(pair_freq[i][j]+alpha)-lgamma(alpha);
	}
      }
      M[l][k]=logP_dep-logP_indep;
      if((l==1) && (k==0)){
	min=M[l][k];
      }
      if(M[l][k]<min){
	min=M[l][k];
      }
      M[k][l]=M[l][k];
    }
  }
  for (int l=0;l<M.ncols();l++){
    for (int k=0; k<l;k++){
      M[k][l]-=min;
      M[l][k]-=min;
    }
  }

 //apply product correction
  double num2=0;
  double tot_average=0;
  for (int l=0;l<M.ncols();l++){
    double num=0;
    for (int k=0; k<M.ncols();k++){
      if(k!=l){
	correction[l]+=M[l][k];
	num++;
      }
      if(k<l){
	tot_average+=M[l][k];
	num2++;
      }
    }
    correction[l]/=num;
  }
  tot_average/=num2;
  for (int l=0;l<M.ncols();l++){
    for (int k=0; k<l;k++){
      double corr=(correction[l]*correction[k])/tot_average;
      M[k][l]-=corr;
      M[l][k]=M[k][l];
     }
  }
  return M;
}

Mat<char> read_alignment(string file){
  int dim[2];
  ifstream in(file.c_str());
  assert(in);
  string line;
  in>>dim[0]>>dim[1];
  Mat<char> M(dim[0], dim[1]);
  for(int i=0;i<M.nrows();i++){
    in>>line;
    for(int j=0;j<M.ncols();j++){
      M[i][j]=line[j];
    }
    line="";//not necessary
  }
  return M;
}

Mat<double> contract_edge(Mat<double>& logR, int node1, int node2){

  Mat<double> M(0., logR.nrows()-1, logR.ncols()-1);
  for(int i=0;i<logR.nrows();i++){
    if((i!=node1) && (i!=node2)){
      for(int j=0;j<i;j++){
	if((j!=node1) && (j!=node2)){
	  int index1=i;
	  int index2=j;
	  //take out rows/columns given by node1 and node2
	  if(index1>=node1){index1--;}
	  if(index1>=node2){index1--;}
	  if(index2>=node1){index2--;}
	  if(index2>=node2){index2--;}
	  if((index1>=0) && (index2>=0)){
	    M[index1][index2]=logR[i][j];
	    M[index2][index1]=logR[i][j];
	  }
	}
      }
    }
  }
  //add contracted edge in the last column
  int last_index=logR.nrows()-2;
  for(int i=0;i<logR.nrows();i++){
    if((i!=node1) && (i!=node2)){
      int current_index=i;
      if(current_index>=node1){current_index--;}
      if(current_index>=node2){current_index--;}
      if(current_index>=0){
	M[last_index][current_index]=logR[node1][i]+logR[node2][i];
	M[current_index][last_index]=M[last_index][current_index];
      }
    }
  }
  return M;
}


//************** main ***************************//
int main(int argc, char* argv[]){
  if (argc!=2){
    cerr<<"usage: alignment"<<endl;
    exit(1);
  } 
  
  Mat<char> alignment=read_alignment(argv[1]);

   Mat<double> M=calculate_logR(alignment);
   rescale_logR(M);
   double logp_D=score(M);

  for(int i=0;i<M.nrows();i++){
    for(int j=0;j<i;j++){
      Mat<double> Q=contract_edge(M, i, j);
      double logp_edge=score(Q)+log(M[i][j]);
      cout<<i<<" "<<j<<" "<<logp_edge-logp_D<<endl;
    }
  }
}
