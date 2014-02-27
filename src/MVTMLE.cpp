#include "MVTMLE.h"
#include <Rcpp.h>
//#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat SqrtS(mat S){
  unsigned int p = S.n_cols;
  vec eigvalS;
  mat eigvecS;
  eig_sym(eigvalS,eigvecS,S);
  vec SqrtEigvalS = sqrt(eigvalS);
  mat SqrtEigvalSM = repmat(SqrtEigvalS,1,p);
  //mat Sqrt = eigvecS % SqrtEigvalSM.t();
  mat Sqrt = eigvecS;
  Sqrt.each_row() %= SqrtEigvalS.t();
  return(Sqrt);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


SEXP MVTMLE0(SEXP x, SEXP nu, SEXP prewhitened, SEXP delta, SEXP maxiter)
  {  
    mat X = as<arma::mat>(x); 
    double NU = as<double>(nu);
    double DELTA = as<double>(delta);
    bool PREWHITENED = as<bool>(prewhitened);
    int MAXITER = as<int>(maxiter);
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    mat B = eye(p,p);
    
    if(!PREWHITENED){
                 mat S0 = X.t()*X/n;
                 B = SqrtS(S0);
                 X = trans(arma::solve(B, X.t()));
                 }
                 
    vec denom = NU + arma::sum(square(X),1);
    //mat Z = X / repmat(sqrt(denom),1,p);
    mat Z = X;
    Z.each_col() /= sqrt(denom);
    mat Psi = (NU + p) * Z.t()*Z / n;
    vec eigvalPsi;
    mat eigvecPsi;
    eig_sym(eigvalPsi,eigvecPsi,Psi);
    double nG = sqrt(sum(pow(1 - eigvalPsi,2)));
	  int iter = 0;
    
    while (nG > DELTA && iter < MAXITER)
          {
          iter = iter+1;
          B = B*eigvecPsi;
          X = X * eigvecPsi;
  	      //Z = square(X)/repmat(denom,1,p);
          Z= square(X);
          Z.each_col() /= denom;
		      mat Ht = diagmat(eigvalPsi) - (NU + p) * Z.t()*Z/n; //+ (NU == 0)/p;
		      vec a = solve(Ht,eigvalPsi - 1);
          vec ahalf = a/2;
		      //mat Xnew = trans(X.t() % repmat(exp(-ahalf),1,n));
          mat Xnew = X;
          Xnew.each_row() %= trans(exp(-ahalf));
		      vec denomnew = NU + sum(square(Xnew),1);
		      double DL = (NU + p)*mean(log(denomnew/denom)) + sum(a);
          double DL0 = sum(a.t()*(1-eigvalPsi))/4;
          if (DL < DL0)
  	          {
			        //B = trans(B.t()  % repmat(exp(ahalf),1,p));
			        B.each_row() %= trans(exp(ahalf));
              X = Xnew;
			        denom = denomnew;
		          } 
          else
              {
              vec sqrteigvalPsi = sqrt(eigvalPsi);
  		        //B = trans(B.t() % repmat(sqrteigvalPsi,1,p));
			        //X = trans(X.t() / repmat(sqrteigvalPsi,1,n));
              B.each_row() %= trans(sqrteigvalPsi);
              X.each_row() /= trans(sqrteigvalPsi);
			        denom = NU + sum(square(X),1);  
              }
          //Z = X / repmat(sqrt(denom),1,p);
          Z=X;
          Z.each_col() /= sqrt(denom);
  	      Psi = (NU + p) * Z.t()*Z / n;
		      eig_sym(eigvalPsi,eigvecPsi,Psi);
		      nG = sqrt(sum(square(1 - eigvalPsi)));
          }  
    mat S = B*B.t();
    
    return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("nu") = NU,
                            Rcpp::Named("delta") = DELTA,
                            Rcpp::Named("prewhitened") = PREWHITENED,
                            Rcpp::Named("B") = B,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("nG") = nG
                            );
    
  }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


mat MVTMLE0cpp(mat X, double NU, bool PREWHITENED, double DELTA, int MAXITER)
  {     
    int p = X.n_cols;
    int n = X.n_rows;
    
    mat B = eye(p,p);
    
    if(!PREWHITENED){
                 mat S0 = X.t()*X/n;
                 B = SqrtS(S0);
                 X = trans(arma::solve(B, X.t()));
                 }
                 
    vec denom = NU + arma::sum(square(X),1);
    //mat Z = X / repmat(sqrt(denom),1,p);
    mat Z=X;
    Z.each_col() /= sqrt(denom);
    mat Psi = (NU + p) * Z.t()*Z / n;
    vec eigvalPsi;
    mat eigvecPsi;
    eig_sym(eigvalPsi,eigvecPsi,Psi);
    double nG = sqrt(sum(pow(1 - eigvalPsi,2)));
    int iter = 0;
    
    while (nG > DELTA && iter < MAXITER)
          {
          iter = iter+1;
          B = B*eigvecPsi;
          X = X * eigvecPsi;
  	      //Z = square(X)/repmat(denom,1,p);
		      Z=square(X);
          Z.each_col() /= denom;
          mat Ht = diagmat(eigvalPsi) - (NU + p) * Z.t()*Z/n; //+ (NU == 0)/p;
		      vec a = solve(Ht,eigvalPsi - 1);
          vec ahalf = a/2;
		      //mat Xnew = trans(X.t() % repmat(exp(-ahalf),1,n));
		      mat Xnew = X;
          Xnew.each_row() %= trans(exp(-ahalf));
          vec denomnew = NU + sum(square(Xnew),1);
		      double DL = (NU + p)*mean(log(denomnew/denom)) + sum(a);
          double DL0 = sum(a.t()*(1-eigvalPsi))/4;
          if (DL < DL0)
  	          {
			        //B = trans(B.t()  % repmat(exp(ahalf),1,p));
			        B.each_row() %= trans(exp(ahalf));
              X = Xnew;
			        denom = denomnew;
		          } 
          else
              {
              vec sqrteigvalPsi = sqrt(eigvalPsi);
  		        //B = trans(B.t() % repmat(sqrteigvalPsi,1,p));
			        //X = trans(X.t() / repmat(sqrteigvalPsi,1,n));
			        B.each_row() %= sqrteigvalPsi.t();
              X.each_row() /= sqrteigvalPsi.t();
              
              denom = NU + sum(square(X),1);  
              }
          //Z = X / repmat(sqrt(denom),1,p);
  	      Z = X;
          Z.each_col() /= sqrt(denom);
          Psi = (NU + p) * Z.t()*Z / n;
		      eig_sym(eigvalPsi,eigvecPsi,Psi);
		      nG = sqrt(sum(square(1 - eigvalPsi)));
          }  
    
    
    return B;
    
  }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

mat LocalPsi(mat X, double nu, int n, int p, int N)
  {
  mat Yi = X.row(n-1)-X.row(n-2);  
  mat Z = Yi/sqrt(nu + accu(square(Yi)));
  mat Psi =  Z.t() * Z ;
  
  for (int i=0; i<n-2; i++)
    {
    mat rowXi = X.row(i);
    //mat Y = X.rows(i+1,n-1)-repmat(rowXi.t(),1,n-1-i).t();
    mat Y = X.rows(i+1,n-1);
    Y.each_row()-=rowXi;
    vec denom = nu + sum(square(Y),1); 
    //mat Z = Y/ repmat(sqrt(denom),1,p);
    mat Z=Y;
    Z.each_col() /= sqrt(denom);
    Psi = Psi + Z.t()*Z;
    }
  Psi = Psi * ((nu+p)/N);
  return Psi;
  }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

mat LocalA(mat X, double nu, int n, int p, int N, vec evs)
{
  mat Yi = X.row(n-1)-X.row(n-2);  
  mat YiSq = square(Yi);
  mat Z = YiSq/(nu + accu(YiSq));
  mat Ht =  Z.t() * Z;
  
  for (int i=0; i<n-2; i++)
  {
    mat rowXi = X.row(i);
    //mat Y = X.rows(i+1,n-1)-repmat(trans( X.row(i)),1,n-1-i).t();
    mat Y = X.rows(i+1,n-1);
    Y.each_row() -= rowXi;
    //mat YSq = square(Y);
    //vec denom = nu + sum(YSq,1); 
    //mat Z = YSq/ repmat(denom,1,p);
    mat Z = square(Y);
    vec denom = nu + sum(Z,1); 
    Z.each_col() /= denom;
    Ht = Ht + Z.t()*Z;
  }
  
  Ht = diagmat(evs) - Ht*((nu+p)/N);
  mat a = solve(Ht,evs-1);
  
  
  return a;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double LocalDL(mat X, double nu, int n, int p, int N, mat Xnew, mat a)
{
  mat Yi = X.row(n-1)-X.row(n-2); 
  mat Ynewi = Xnew.row(n-1)-Xnew.row(n-2);
  double DL = log((nu + accu(square(Ynewi)))/(nu + accu(square(Yi))));
  
  for (int i=0; i<n-2; i++)
  {
    //mat Y = X.rows(i+1,n-1)-repmat(trans(X.row(i)),1,n-1-i).t();
    //mat Ynew = Xnew.rows(i+1,n-1)-repmat(trans(Xnew.row(i)),1,n-1-i).t();
    mat Y = X.rows(i+1,n-1);
    Y.each_row()-=X.row(i);
    mat Ynew = Xnew.rows(i+1,n-1);
    Ynew.each_row() -= Xnew.row(i);
    vec denom = nu + sum(square(Y),1); 
    vec denomNew = nu + sum(square(Ynew),1);
    DL = DL + accu(log(denomNew/denom));
  }
  
  
  DL = DL*((nu+p)/N) + accu(a);
  
  return DL;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

mat pairdiff(mat X)
  {
    int n=X.n_rows; 
    int p=X.n_cols; 
    int i;
    int j;
    int rowRes=0;
    mat DIFF=zeros(n*(n-1)/2,p);
    for (i=0; i<(n-1); i++)
      {
      for (j=(i+1); j<n; j++)
        {
        DIFF.row(rowRes) = X.row(i)-X.row(j);
        rowRes = rowRes + 1;
        }
    }
    return(DIFF);
  }

/*
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP MVTMLEsymm(SEXP x, SEXP nu, SEXP delta, SEXP maxiter, SEXP nmax)
  {  
    mat X = as<arma::mat>(x);
    double NU = as<double>(nu);
    double DELTA = as<double>(delta);
    int MAXITER = as<int>(maxiter);
    int NMAX = as<int>(nmax);
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    mat S0 = 2*cov(X);
    mat B = SqrtS(S0);
    mat Xs = trans(arma::solve(B, X.t()));
    
    mat Xrev = flipud(Xs);
    mat X0 = Xs-Xrev;
    
    mat C = MVTMLE0cpp(X0, NU, true, DELTA, MAXITER);
    B = B*C;
    Xs = trans(arma::solve(C, X.t()));
    
    if (n < NMAX)
        {
          mat XXs = pairdiff(Xs);
          C = MVTMLE0cpp(XXs, NU, true, DELTA, MAXITER);
          B = B*C;
        }
    else
        {
          int N = n*(n-1)/2;
          mat Psi = LocalPsi(Xs,NU,n,p,N);
          vec TMPeigvalPsi;
          mat TMPeigvecPsi;
          eig_sym(TMPeigvalPsi,TMPeigvecPsi,Psi);
          int nG = sqrt(accu(square(1-TMPeigvalPsi)));
          
          int iter = 0;
          
          while (nG > DELTA && iter < MAXITER)
            {
              iter = iter+1;
              B = B*TMPeigvecPsi;
              Xs = Xs * TMPeigvecPsi;
              mat a = LocalA(Xs,NU,n,p,N,TMPeigvalPsi);
              mat ahalf = a/2;
              mat Xsnew = trans(Xs.t() % repmat(exp(-ahalf),1,n));

		          double DL = LocalDL(Xs,NU,n,p,N,Xsnew,a);
              double DL0 = sum(a.t()*(1-TMPeigvalPsi))/4;
              
              if (DL < DL0)
                {
			          B = trans(B.t()  % repmat(exp(ahalf),1,p));
			          Xs = Xsnew;
		            } 
              else
                {
                vec sqrtTMPeigvalPsi = sqrt(TMPeigvalPsi);
  		          B = trans(B.t() % repmat(sqrtTMPeigvalPsi,1,p));
			          Xs = trans(Xs.t() / repmat(TMPeigvalPsi,1,n));
                }
              Psi = LocalPsi(Xs,NU,n,p,N);
        	    eig_sym(TMPeigvalPsi,TMPeigvecPsi,Psi);
		          nG = sqrt(accu(square(1-TMPeigvalPsi))); 
            }
        }
          mat S = B*B.t();
          return Rcpp::List::create(Rcpp::Named("B") = B,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("nG") = nG
                            );
        
  }
 */
  

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP MVTMLEsymm1(SEXP x, SEXP nu, SEXP delta, SEXP maxiter)
  {  
    mat X = as<arma::mat>(x);
    double NU = as<double>(nu);
    double DELTA = as<double>(delta);
    int MAXITER = as<int>(maxiter);
    
    int p = X.n_cols;
    int n = X.n_rows;
  
    mat S0 = 2*cov(X);
    mat B = SqrtS(S0);
   
    mat Xs = trans(arma::solve(B, X.t()));
   
    mat Xperm = Xs;
    Xperm.shed_row(0);
    Xperm.insert_rows(n-1,Xs.row(0));
    mat X0 = Xs-Xperm;
  
    
    mat C = MVTMLE0cpp(X0, NU, true, DELTA, MAXITER);
 
    B = B*C;
  
    Xs = trans(arma::solve(C, Xs.t()));
    
    mat XXs = pairdiff(Xs);
    int N = XXs.n_rows;
    
  
    mat C2 = eye(p,p);
    vec denom = NU + arma::sum(square(XXs),1);
    //mat Z = XXs / repmat(sqrt(denom),1,p);
    mat Z = XXs;
    Z.each_col() /= sqrt(denom);
    
    mat Psi = (NU + p) * Z.t()*Z / N;
    vec eigvalPsi;
    mat eigvecPsi;
    eig_sym(eigvalPsi,eigvecPsi,Psi);
    double nG = sqrt(sum(pow(1 - eigvalPsi,2)));
    int iter = 0;
    
    while (nG > DELTA && iter < MAXITER)
          {
          iter = iter+1;
          C2 = C2*eigvecPsi;
          XXs = XXs * eigvecPsi;
          //Z = square(XXs)/repmat(denom,1,p);
          Z=square(XXs);
          Z.each_col() /= denom;
          mat Ht = diagmat(eigvalPsi) - (NU + p) * Z.t()*Z/N; //+ (NU == 0)/p;
		      vec a = solve(Ht,eigvalPsi - 1);
          vec ahalf = a/2;
		      //mat XXsnew = trans(XXs.t() % repmat(exp(-ahalf),1,N));
          mat XXsnew = XXs;
          XXsnew.each_row() %= trans(exp(-ahalf));
		      vec denomnew = NU + sum(square(XXsnew),1);
		      double DL = (NU + p)*mean(log(denomnew/denom)) + sum(a);
          double DL0 = sum(a.t()*(1-eigvalPsi))/4;
          if (DL < DL0)
  	          {
			        //C2 = trans(C2.t()  % repmat(exp(ahalf),1,p));
              C2.each_row() %= trans(exp(ahalf));
			        XXs = XXsnew;
			        denom = denomnew;
		          } 
          else
              {
              vec sqrteigvalPsi = sqrt(eigvalPsi);
  		        //C2 = trans(C2.t() % repmat(sqrteigvalPsi,1,p));
			        //XXs = trans(XXs.t() / repmat(sqrteigvalPsi,1,N));
			        C2.each_row() %= sqrteigvalPsi.t();
              XXs.each_row() /= sqrteigvalPsi.t();
              denom = NU + sum(square(XXs),1);  
              }
          //Z = XXs / repmat(sqrt(denom),1,p);
          Z=XXs;
          Z.each_col() /= sqrt(denom);
  	      Psi = (NU + p) * Z.t()*Z / N;
		      eig_sym(eigvalPsi,eigvecPsi,Psi);
		      nG = sqrt(sum(square(1 - eigvalPsi)));
          } 
          
    B = B*C2;
        
    
    mat S = B*B.t();
    return Rcpp::List::create(Rcpp::Named("B") = B,
           Rcpp::Named("S") = S,
           Rcpp::Named("iter") = iter,
           Rcpp::Named("nG") = nG
          );
  }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP MVTMLEsymm2(SEXP x, SEXP nu, SEXP delta, SEXP maxiter)
  {  
    mat X = as<arma::mat>(x);
    double NU = as<double>(nu);
    double DELTA = as<double>(delta);
    int MAXITER = as<int>(maxiter);
    //int NMAX = as<int>(nmax);
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    mat S0 = 2*cov(X);
    mat B = SqrtS(S0);
    mat Xs = trans(arma::solve(B, X.t()));
    
    mat Xperm = Xs;
    Xperm.shed_row(0);
    Xperm.insert_rows(n-1,Xs.row(0));
    mat X0 = Xs-Xperm;
    
    mat C = MVTMLE0cpp(X0, NU, true, DELTA, MAXITER);
    B = B*C;
    Xs = trans(arma::solve(C, Xs.t()));
    

    int N = n*(n-1)/2.0;
    mat Psi = LocalPsi(Xs,NU,n,p,N);
    vec TMPeigvalPsi;
    mat TMPeigvecPsi;
    eig_sym(TMPeigvalPsi,TMPeigvecPsi,Psi);
    double nG = sqrt(accu(square(1-TMPeigvalPsi)));
          
    int iter = 0;
          
    while (nG > DELTA && iter < MAXITER)
        {
        iter = iter+1;
        B = B*TMPeigvecPsi;
        Xs = Xs * TMPeigvecPsi;
        mat a = LocalA(Xs,NU,n,p,N,TMPeigvalPsi);
        mat ahalf = a/2;
        //mat Xsnew = trans(Xs.t() % repmat(exp(-ahalf),1,n));
        mat Xsnew = Xs;
        Xsnew.each_row() %= trans(exp(-ahalf));
        
  	    double DL = LocalDL(Xs,NU,n,p,N,Xsnew,a);
        double DL0 = sum(a.t()*(1-TMPeigvalPsi))/4;
              
              if (DL < DL0)
                {
			          //B = trans(B.t()  % repmat(exp(ahalf),1,p));
                B.each_row() %= trans(exp(ahalf));
			          Xs = Xsnew;
		            } 
              else
                {
                vec sqrtTMPeigvalPsi = sqrt(TMPeigvalPsi);
  		          //B = trans(B.t() % repmat(sqrtTMPeigvalPsi,1,p));
			          //Xs = trans(Xs.t() / repmat(sqrtTMPeigvalPsi,1,n));
                B.each_row() %= trans(sqrtTMPeigvalPsi);
                Xs.each_row() /= trans(sqrtTMPeigvalPsi);
                }
              Psi = LocalPsi(Xs,NU,n,p,N);
        	    eig_sym(TMPeigvalPsi,TMPeigvecPsi,Psi);
		          nG = sqrt(accu(square(1-TMPeigvalPsi))); 
          }
      mat S = B*B.t();
      return Rcpp::List::create(Rcpp::Named("B") = B,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("nG") = nG
                            );
      }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


SEXP Tyler0(SEXP x, SEXP prewhitened, SEXP delta, SEXP maxiter)
  {  
    mat X = as<arma::mat>(x); 
    double DELTA = as<double>(delta);
    bool PREWHITENED = as<bool>(prewhitened);
    int MAXITER = as<int>(maxiter);
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    mat B = eye(p,p);
    
    if(!PREWHITENED){
                 mat S0 = X.t()*X/n;
                 B = SqrtS(S0);
                 X = trans(arma::solve(B, X.t()));
                 }
                 
    vec denom = arma::sum(square(X),1);
    //mat Z = X / repmat(sqrt(denom),1,p);
    mat Z = X;
    Z.each_col() /= sqrt(denom);
    mat Psi = p * Z.t()*Z / n;
    vec eigvalPsi;
    mat eigvecPsi;
    eig_sym(eigvalPsi,eigvecPsi,Psi);
    double nG = sqrt(sum(pow(1 - eigvalPsi,2)));
    int iter = 0;
    
    while (nG > DELTA && iter < MAXITER)
          {
          iter = iter+1;
          B = B*eigvecPsi;
          X = X * eigvecPsi;
          //Z = square(X)/repmat(denom,1,p);
          Z= square(X);
          Z.each_col() /= denom;
		      mat Ht = diagmat(eigvalPsi) - p * Z.t()*Z/n + 1/(double)p;
		      vec a = solve(Ht,eigvalPsi - 1);
          vec ahalf = a/2;
		      //mat Xnew = trans(X.t() % repmat(exp(-ahalf),1,n));
          mat Xnew = X;
          Xnew.each_row() %= trans(exp(-ahalf));
		      vec denomnew =  sum(square(Xnew),1);
		      double DL = p * mean(log(denomnew/denom)) + sum(a);
          double DL0 = sum(a.t()*(1-eigvalPsi))/4;
          if (DL < DL0)
  	          {
			        //B = trans(B.t()  % repmat(exp(ahalf),1,p));
			        B.each_row() %= trans(exp(ahalf));
              X = Xnew;
			        denom = denomnew;
		          } 
          else
              {
              vec sqrteigvalPsi = sqrt(eigvalPsi);
  		        //B = trans(B.t() % repmat(sqrteigvalPsi,1,p));
			        //X = trans(X.t() / repmat(sqrteigvalPsi,1,n));
              B.each_row() %= trans(sqrteigvalPsi);
              X.each_row() /= trans(sqrteigvalPsi);
			        denom = sum(square(X),1);  
              }
          //Z = X / repmat(sqrt(denom),1,p);
          Z=X;
          Z.each_col() /= sqrt(denom);
  	      Psi = p * Z.t()*Z / n;
		      eig_sym(eigvalPsi,eigvecPsi,Psi);
		      nG = sqrt(sum(square(1 - eigvalPsi)));
          }  
    mat S1 = B*B.t();
    mat S = S1 / pow(det(S1),(1/(double)p));
    
    return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("delta") = DELTA,
                            Rcpp::Named("prewhitened") = PREWHITENED,
                            Rcpp::Named("B") = B,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("nG") = nG
                            );
    
  }
  



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

mat Tyler0cpp(mat X, bool PREWHITENED, double DELTA, int MAXITER)
  {     
    int p = X.n_cols;
    int n = X.n_rows;
    
    mat B = eye(p,p);
    
    if(!PREWHITENED){
                 mat S0 = X.t()*X/n;
                 B = SqrtS(S0);
                 X = trans(arma::solve(B, X.t()));
                 }
                 
    vec denom = arma::sum(square(X),1);
    //mat Z = X / repmat(sqrt(denom),1,p);
    mat Z=X;
    Z.each_col() /= sqrt(denom);
    mat Psi = p * Z.t()*Z / n;
    vec eigvalPsi;
    mat eigvecPsi;
    eig_sym(eigvalPsi,eigvecPsi,Psi);
    double nG = sqrt(sum(pow(1 - eigvalPsi,2)));
    int iter = 0;
    
    while (nG > DELTA && iter < MAXITER)
          {
          iter = iter+1;
          B = B*eigvecPsi;
          X = X * eigvecPsi;
          //Z = square(X)/repmat(denom,1,p);
		      Z=square(X);
          Z.each_col() /= denom;
          mat Ht = diagmat(eigvalPsi) -  p * Z.t()*Z/n + 1/(double)p;
		      vec a = solve(Ht,eigvalPsi - 1);
          vec ahalf = a/2;
		      //mat Xnew = trans(X.t() % repmat(exp(-ahalf),1,n));
		      mat Xnew = X;
          Xnew.each_row() %= trans(exp(-ahalf));
          vec denomnew = sum(square(Xnew),1);
		      double DL = p*mean(log(denomnew/denom)) + sum(a);
          double DL0 = sum(a.t()*(1-eigvalPsi))/4;
          if (DL < DL0)
  	          {
			        //B = trans(B.t()  % repmat(exp(ahalf),1,p));
			        B.each_row() %= trans(exp(ahalf));
              X = Xnew;
			        denom = denomnew;
		          } 
          else
              {
              vec sqrteigvalPsi = sqrt(eigvalPsi);
  		        //B = trans(B.t() % repmat(sqrteigvalPsi,1,p));
			        //X = trans(X.t() / repmat(sqrteigvalPsi,1,n));
			        B.each_row() %= sqrteigvalPsi.t();
              X.each_row() /= sqrteigvalPsi.t();
              
              denom = sum(square(X),1);  
              }
          //Z = X / repmat(sqrt(denom),1,p);
  	      Z = X;
          Z.each_col() /= sqrt(denom);
          Psi =  p * Z.t()*Z / n;
		      eig_sym(eigvalPsi,eigvecPsi,Psi);
		      nG = sqrt(sum(square(1 - eigvalPsi)));
          }  
    
    
    return B;
    
  }
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

mat TylerLocalPsi(mat X, int n, int p, int N)
  {

  mat Yi = X.row(n-1)-X.row(n-2);  
  mat Z = Yi/sqrt(accu(square(Yi)));
  mat Psi =  Z.t() * Z ;
  
  for (int i=0; i<n-2; i++)
    {
    mat rowXi = X.row(i);
    //mat Y = X.rows(i+1,n-1)-repmat(rowXi.t(),1,n-1-i).t();
    mat Y = X.rows(i+1,n-1);
    Y.each_row()-=rowXi;
    vec denom = sum(square(Y),1); 
    //mat Z = Y/ repmat(sqrt(denom),1,p);
    mat Z=Y;
    Z.each_col() /= sqrt(denom);
    Psi = Psi + Z.t()*Z;
    }

  Psi = Psi * ((double)p/N);
  return Psi;
  }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat TylerLocalA(mat X, int n, int p, int N, vec evs)
{
  mat Yi = X.row(n-1)-X.row(n-2);  
  mat YiSq = square(Yi);
  mat Z = YiSq/(accu(YiSq));
  mat Ht =  Z.t() * Z;
  
  for (int i=0; i<n-2; i++)
  {
    mat rowXi = X.row(i);
    //mat Y = X.rows(i+1,n-1)-repmat(trans( X.row(i)),1,n-1-i).t();
    mat Y = X.rows(i+1,n-1);
    Y.each_row() -= rowXi;
    //mat YSq = square(Y);
    //vec denom =  sum(YSq,1); 
    //mat Z = YSq/ repmat(denom,1,p);
    mat Z = square(Y);
    vec denom = sum(Z,1); 
    Z.each_col() /= denom;
    Ht = Ht + Z.t()*Z;
  }
  Ht = diagmat(evs) - Ht*((double)p/N)+1/(double)p;
  mat a = solve(Ht,evs-1);
  
  
  return a;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double TylerLocalDL(mat X, int n, int p, int N, mat Xnew, mat a)
{
  mat Yi = X.row(n-1)-X.row(n-2); 
  mat Ynewi = Xnew.row(n-1)-Xnew.row(n-2);
  double DL = log((accu(square(Ynewi)))/(accu(square(Yi))));
  
  for (int i=0; i<n-2; i++)
  {
    //mat Y = X.rows(i+1,n-1)-repmat(trans(X.row(i)),1,n-1-i).t();
    //mat Ynew = Xnew.rows(i+1,n-1)-repmat(trans(Xnew.row(i)),1,n-1-i).t();
    mat Y = X.rows(i+1,n-1);
    Y.each_row()-=X.row(i);
    mat Ynew = Xnew.rows(i+1,n-1);
    Ynew.each_row() -= Xnew.row(i);
    vec denom = sum(square(Y),1); 
    vec denomNew = sum(square(Ynew),1);
    DL = DL + accu(log(denomNew/denom));
  }
  
  
  DL = DL*((double)p/N) + accu(a);
  
  return DL;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP Tylersymm1(SEXP x, SEXP delta, SEXP maxiter)
  {  
    mat X = as<arma::mat>(x);
    double DELTA = as<double>(delta);
    int MAXITER = as<int>(maxiter);
    
    int p = X.n_cols;
    int n = X.n_rows;
  
    mat S0 = 2*cov(X);
    mat B = SqrtS(S0);
    
    mat Xs = trans(arma::solve(B, X.t()));
    
    mat Xperm = Xs;
    Xperm.shed_row(0);
    Xperm.insert_rows(n-1,Xs.row(0));
    mat X0 = Xs-Xperm;
    
    mat C = Tyler0cpp(X0, true, DELTA, MAXITER);
   
    B = B*C;
  
    Xs = trans(arma::solve(C, Xs.t()));
    
    mat XXs = pairdiff(Xs);
    int N = XXs.n_rows;
    
  
    mat C2 = eye(p,p);
    vec denom =  arma::sum(square(XXs),1);
    //mat Z = XXs / repmat(sqrt(denom),1,p);
    mat Z = XXs;
    Z.each_col() /= sqrt(denom);
    mat Psi =  p * Z.t()*Z / N;
    vec eigvalPsi;
    mat eigvecPsi;
    eig_sym(eigvalPsi,eigvecPsi,Psi);
    double nG = sqrt(sum(pow(1 - eigvalPsi,2)));
    int iter = 0;
    
    while (nG > DELTA && iter < MAXITER)
          {
          iter = iter+1;
          C2 = C2*eigvecPsi;
          XXs = XXs * eigvecPsi;
          //Z = square(XXs)/repmat(denom,1,p);
          Z=square(XXs);
          Z.each_col() /= denom;
          mat Ht = diagmat(eigvalPsi) - p * Z.t()*Z/N + 1/double(p);
  	      vec a = solve(Ht,eigvalPsi - 1);
          vec ahalf = a/2;
		      //mat XXsnew = trans(XXs.t() % repmat(exp(-ahalf),1,N));
		      mat XXsnew = XXs;
          XXsnew.each_row() %= trans(exp(-ahalf));
          vec denomnew = sum(square(XXsnew),1);
		      double DL = p*mean(log(denomnew/denom)) + sum(a);
          double DL0 = sum(a.t()*(1-eigvalPsi))/4;
          if (DL < DL0)
  	          {
			        //C2 = trans(C2.t()  % repmat(exp(ahalf),1,p));
              C2.each_row() %= trans(exp(ahalf));
			        XXs = XXsnew;
			        denom = denomnew;
		          } 
          else
              {
              vec sqrteigvalPsi = sqrt(eigvalPsi);
  		        //C2 = trans(C2.t() % repmat(sqrteigvalPsi,1,p));
			        //XXs = trans(XXs.t() / repmat(sqrteigvalPsi,1,N));
			        C2.each_row() %= sqrteigvalPsi.t();
              XXs.each_row() /= sqrteigvalPsi.t();
              denom = sum(square(XXs),1);  
              }
          //Z = XXs / repmat(sqrt(denom),1,p);
          Z=XXs;
          Z.each_col() /= sqrt(denom);
  	      Psi = p * Z.t()*Z / N;
		      eig_sym(eigvalPsi,eigvecPsi,Psi);
		      nG = sqrt(sum(square(1 - eigvalPsi)));
          } 
          
    B = B*C2;
        
    
    mat S1 = B*B.t();
    mat S = S1 / pow(det(S1),(1/(double)p));
    return Rcpp::List::create(Rcpp::Named("B") = B,
           Rcpp::Named("S") = S,
           Rcpp::Named("iter") = iter,
           Rcpp::Named("nG") = nG
          );
  }

  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP Tylersymm2(SEXP x, SEXP delta, SEXP maxiter)
{  
  mat X = as<arma::mat>(x);
  double DELTA = as<double>(delta);
  int MAXITER = as<int>(maxiter);
  //int NMAX = as<int>(nmax);
  
  int p = X.n_cols;
  int n = X.n_rows;
  
  mat S0 = 2*cov(X);
  mat B = SqrtS(S0);
  mat Xs = trans(arma::solve(B, X.t()));
  mat Xperm = Xs;
  Xperm.shed_row(0);
  Xperm.insert_rows(n-1,Xs.row(0));
  mat X0 = Xs-Xperm;
  
  mat C = Tyler0cpp(X0, true, DELTA, MAXITER);
  B = B*C;
  Xs = trans(arma::solve(C, Xs.t()));
 
  int N = n*(n-1)/2.0;
  mat Psi = TylerLocalPsi(Xs,n,p,N);
  vec TMPeigvalPsi;
  mat TMPeigvecPsi;
  eig_sym(TMPeigvalPsi,TMPeigvecPsi,Psi);
  double nG = sqrt(accu(square(1-TMPeigvalPsi)));
  
  int iter = 0;
  
  while (nG > DELTA && iter < MAXITER)
  {
    iter = iter+1;
    B = B*TMPeigvecPsi;
    Xs = Xs * TMPeigvecPsi;
      
    mat a = TylerLocalA(Xs,n,p,N,TMPeigvalPsi);
    mat ahalf = a/2;
    //mat Xsnew = trans(Xs.t() % repmat(exp(-ahalf),1,n));
    mat Xsnew = Xs;
    Xsnew.each_row() %= trans(exp(-ahalf));
    
    double DL = TylerLocalDL(Xs,n,p,N,Xsnew,a);
    double DL0 = sum(a.t()*(1-TMPeigvalPsi))/4;
    
    if (DL < DL0)
    {
      //B = trans(B.t()  % repmat(exp(ahalf),1,p));
      B.each_row() %= trans(exp(ahalf));
      Xs = Xsnew;
    } 
    else
    {
      vec sqrtTMPeigvalPsi = sqrt(TMPeigvalPsi);
      //B = trans(B.t() % repmat(sqrtTMPeigvalPsi,1,p));
      //Xs = trans(Xs.t() / repmat(sqrtTMPeigvalPsi,1,n));
      B.each_row() %= trans(sqrtTMPeigvalPsi);
      Xs.each_row() /= trans(sqrtTMPeigvalPsi);
    }
    Psi = TylerLocalPsi(Xs,n,p,N);
    eig_sym(TMPeigvalPsi,TMPeigvecPsi,Psi);
    nG = sqrt(accu(square(1-TMPeigvalPsi))); 
  }
  mat S1 = B*B.t();
  mat S = S1 / pow(det(S1),(1/(double)p));
  
  return Rcpp::List::create(Rcpp::Named("B") = B,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("nG") = nG
  );
}
