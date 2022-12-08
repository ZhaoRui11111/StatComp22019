#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

//' @title A function helps to get terms using Rcpp
//' @description A function helps to get terms using Rcpp
//' @param x the word in one document
//' @param vocab the all word list
//' @return a index
//' 
//' @export
// [[Rcpp::export]]
IntegerVector get_termsC(CharacterVector x,CharacterVector vocab) {
  IntegerVector index = match(x, vocab);
  index = index[!is_na(index)];
  index = index-1;
  return(index);
}


