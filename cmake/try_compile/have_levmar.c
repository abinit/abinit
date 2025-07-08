#include <levmar.h>

void dfit_function(double *p, double *y, int m, int n, void *adata)
{
 p = 0;
}

int main(int argc, char* argv[])
{
 int ret;
 int c_npoles = 1;
 int c_nvals  = 1;
 int nparam   = 1;

 double adata[1];
 double p[1];
 double yvals[1];
 double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

 ret=dlevmar_dif(dfit_function, p, yvals, nparam, c_nvals, 5000, \
        opts, info, NULL, NULL, (void *)&adata);
}
