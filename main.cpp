#include <stdio.h>
#include <math.h>

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/builtin-defun-decls.h>
#include <octave/CMatrix.h>
#include <octave/Array.h>


//vb = a, mu, sigma

template <int dim_b, int dim_data>
class GaussNewton
{
    double g(double x, RowVector b)
    {
        double A = b(0);
        double mu = b(1);
        double s = b(2);
        double ret = -(x - mu)*(x - mu)/(2*s*s);
        return ret;
    }

public:
    double f(double x, RowVector b)
    {
        return b(0) * exp(g(x, b));
    }



};



int main(int argc, const char * argv[])
{
    GaussNewton<-1,-1> gn;

    for (double x = -5; x < 5; x += 0.05)
    {
        RowVector b(3);
        b(0) = 1;
        b(1) = 0;
        b(2) = 1;

        double vf = gn.f(x, b);
        printf("%.3lf\n", vf);

    }

    return 0;
}
