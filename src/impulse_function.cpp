
#include "hydrorecipes.h"

//==============================================================================
//' @title
//' impulse_function
//'
//' @description
//' Calculation of the impulse function from a well function.
//'
//' @param u well function
//' @param flow_time_interval time between flow rate measurements in samples
//'
//' @return impulse function for convolution
//'
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<double> impulse_function(std::vector<double> u,
                                     int flow_time_interval)
{
    
    int n = u.size();

    if (flow_time_interval >= n) {
       flow_time_interval = n - 1;
    }

    // calculate the pulse
    for (int i = (n - flow_time_interval - 1); i > -1; --i)
    {
        u[i + flow_time_interval] -= u[i];
    }

    return u;
}