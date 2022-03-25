// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// extract
int extract(const std::string input_file_prefix, const std::string class_label, int max_number_of_frames, int maximum_boundary_length, int maximum_cell_area, int cooccurrence_levels, int number_of_wavelet_levels);
RcppExport SEXP _CellPhe_extract(SEXP input_file_prefixSEXP, SEXP class_labelSEXP, SEXP max_number_of_framesSEXP, SEXP maximum_boundary_lengthSEXP, SEXP maximum_cell_areaSEXP, SEXP cooccurrence_levelsSEXP, SEXP number_of_wavelet_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type input_file_prefix(input_file_prefixSEXP);
    Rcpp::traits::input_parameter< const std::string >::type class_label(class_labelSEXP);
    Rcpp::traits::input_parameter< int >::type max_number_of_frames(max_number_of_framesSEXP);
    Rcpp::traits::input_parameter< int >::type maximum_boundary_length(maximum_boundary_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type maximum_cell_area(maximum_cell_areaSEXP);
    Rcpp::traits::input_parameter< int >::type cooccurrence_levels(cooccurrence_levelsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_wavelet_levels(number_of_wavelet_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(extract(input_file_prefix, class_label, max_number_of_frames, maximum_boundary_length, maximum_cell_area, cooccurrence_levels, number_of_wavelet_levels));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CellPhe_extract", (DL_FUNC) &_CellPhe_extract, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_CellPhe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}