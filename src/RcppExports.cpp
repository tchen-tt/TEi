// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// processbam
void processbam(std::string bamfile, std::string outfq, double quantile, int length, int tsd);
RcppExport SEXP _TEi_processbam(SEXP bamfileSEXP, SEXP outfqSEXP, SEXP quantileSEXP, SEXP lengthSEXP, SEXP tsdSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bamfile(bamfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type outfq(outfqSEXP);
    Rcpp::traits::input_parameter< double >::type quantile(quantileSEXP);
    Rcpp::traits::input_parameter< int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< int >::type tsd(tsdSEXP);
    processbam(bamfile, outfq, quantile, length, tsd);
    return R_NilValue;
END_RCPP
}
// processSam2bed
void processSam2bed(std::string alignmentfile, std::string outbedfile, float ratio);
RcppExport SEXP _TEi_processSam2bed(SEXP alignmentfileSEXP, SEXP outbedfileSEXP, SEXP ratioSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type alignmentfile(alignmentfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type outbedfile(outbedfileSEXP);
    Rcpp::traits::input_parameter< float >::type ratio(ratioSEXP);
    processSam2bed(alignmentfile, outbedfile, ratio);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TEi_processbam", (DL_FUNC) &_TEi_processbam, 5},
    {"_TEi_processSam2bed", (DL_FUNC) &_TEi_processSam2bed, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_TEi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
