#include "../all.h"

using namespace std;

int sdz(int n, double err,double* v)
{
    double a = 0.;
    int info_sdz = 0;//0 : ok , 1 : error occured
    a = 1. / cblas_dnrm2(n, v, 1);
    cblas_dscal(n, a, v, 1);

    // Check norm
    double norm = 0.;
    norm = cblas_dnrm2(n, v, 1);
    double eps = abs(1.0 - norm);
    if (eps > 1.0e-15)
    {
        info_sdz = 1;
        printf("\x1b[31m NORM ERROR\x1b[39m\n");
    }
    return info_sdz;
}