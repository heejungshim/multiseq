#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix cxxParentTItable(SEXP sig){
  IntegerMatrix signal=sig; 
  int nsig=(int) signal.nrow();
  int n=(int) signal.ncol();
  int J=(int) log2((double)n);
  
  IntegerMatrix parent(nsig,J*2*n);
  IntegerMatrix TItable(J+1,n);
  //List listTItable;
  for (int s=0; s<nsig; s++){
    TItable(0,_) = signal(s,_);
    for (int D=0; D<J; D++){
      int nD=(int) pow(2., (double) (J-D)), pD=(int) pow(2.,(double) D);
      for (int l=0; l<pD; l++){
	int a=l*nD+1, b=2*l*nD+1, c=2*D*n+b, d;
	for (int i=0; i<nD-1; i++){//todo: check what happens for nD-1<=0 
	  d=TItable(D,a+i-1);
	  parent(s,c+i-1)=d;
	  parent(s,c+i+nD)=d;
	}
	//i=nD-1
	d=TItable(D,a+nD-2);
	parent(s,c+nD-2)=d;
	parent(s,c+nD-1)=d;
	
	for (int i=0; i<nD; i++)
	  TItable(D+1,a+i-1)=parent(s,c+2*i-1)+parent(s,c+2*i);
      }
    }
    //listTItable.push_back(TItable);
  }
  //if you need the list of TI tables
  //return(List::create(Named("listTItable")=listTItable, Named("parent")=parent));
  return(parent);
}

// [[Rcpp::export]]
NumericVector cxxSParentTItable(SEXP sig){
  NumericVector signal=sig; 
  int n=(int) signal.size();
  int J=(int) log2((double)n);

  NumericVector parent(2*J*n);
  NumericMatrix TItable(J+1,n);
  TItable(0,_) = signal;
  for (int D=0; D<J; D++){
     int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
     for (int l=0; l<pD; l++){
        int a=l*nD+1, b=2*l*nD+1, c=2*D*n+b, d;
        for (int i=0; i<nD-1; i++){
           d=(int) TItable(D,a+i-1);
           parent(c+i-1)=d;
           parent(c+i+nD)=d;
        }
      //i=nD-1
      d=(int) TItable(D,a+nD-2);
      parent(c+nD-2)=d;
      parent(c+nD-1)=d;


      for (int i=0; i<nD; i++)
        TItable(D+1,a+i-1)=parent(c+2*i-1)+parent(c+2*i);
    }
  }
  return(parent);
}


// [[Rcpp::export]]
NumericVector cxxreverse_pwave(SEXP estimate, SEXP pmat, SEXP qmat){
  NumericMatrix pp=pmat;
  NumericMatrix qq=qmat;
  NumericVector est1=estimate;
  int np=(int) pp.ncol();
  int J=(int) pp.nrow();
  NumericVector est(np,est1(0));
  //for(int D=J; D-->0;){
  for(int D=0; D<J; D++){
    //int nD=pow(2., (int) (J-D)), pD=pow(2., (int) D);
    int nD=(int) pow(2., (int) (D+1)), pD=(int) pow(2., (int) (J-1-D));
    int nDo2=nD/2;
    NumericVector tempvecl(nD), tempvecr(nD);
    for(int l=0; l<pD; l++){
      int a=l*nD+1; 
      double dep, deq, dp, dq;
      for (int i=0; i<nDo2; i++){
	dep=est(a+i-1);
	//dp=pp(D,a+i-1);
	dp=pp(J-1-D,a+i-1);
	//dq=qq(D,a+i-1);
	dq=qq(J-1-D,a+i-1);
	tempvecl(2*i)=dep+dp;
	tempvecl(2*i+1)=dep+dq;
      }
      for (int i=nDo2; i<nD-1; i++){
	dep=est(a+i);
	//dp=pp(D,a+i);
	dp=pp(J-1-D,a+i);
	deq=est(a+i-1);
	//dq=qq(D,a+i-1);
	dq=qq(J-1-D,a+i-1);
	tempvecr(2*(i-nDo2))=deq+dq;
	tempvecr(2*(i-nDo2)+1)=dep+dp;
      }         
      //i=nD-1
      dep=est(a+nDo2-1);
      //dp=pp(D,a+nDo2-1);
      dp=pp(J-1-D,a+nDo2-1);
      deq=est(a+nD-2); 
      //dq=qq(D,a+nD-2);
      dq=qq(J-1-D,a+nD-2);
      tempvecr(nD-2)=deq+dq;
      tempvecr(nD-1)=dep+dp;
      for(int i=0; i<nD; i++)
               est(a+i-1)=0.5*(tempvecl(i)+tempvecr(i));
    }
  }
  return(est);
}
