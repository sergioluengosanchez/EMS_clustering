#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double poly_volume(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol= x.ncol();

  NumericVector p1(3);
  NumericVector p2(3);
  NumericVector p3(3);

  double total=0;
	for(int i = 0; i < nrow; i+=3) {
		for(int j = 0; j < ncol; ++j) {
  		p1[j]= x(i,j);
  		p2[j]= x((i+1),j);
  		p3[j]= x((i+2),j);
	  }

	  total+=(-(p3[0]*p2[1]*p1[2])+(p2[0]*p3[1]*p1[2])+(p3[0]*p1[1]*p2[2])-(p1[0]*p3[1]*p2[2])-(p2[0]*p1[1]*p3[2])+(p1[0]*p2[1]*p3[2]))/6;
  }
  return abs(total);
}

