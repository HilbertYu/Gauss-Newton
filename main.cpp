#include <stdio.h>
#include <math.h>

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/builtin-defun-decls.h>
#include <octave/CMatrix.h>
#include <octave/Array.h>


//vb = A, mu, sigma

template <int dim_b, int dim_data>
class GaussNewton
{
    double g(double x, RowVector b)
    {
     //   double A = b(0);
        double mu = b(1);
        double s = b(2);
        double ret = -(x - mu)*(x - mu)/(2*s*s);
        return ret;
    }

    double g_mu(double x, RowVector b)
    {
        double mu = b(1);
        double s = b(2);

        return (1.0/s)*(x - mu);
    }

    double g_s(double x, RowVector b)
    {
        double mu = b(1);
        double s = b(2);
        return (1.0/(s*s*s))*(x - mu)*(x - mu);
    }
public:
    double f(double x, RowVector b)
    {
        return b(0) * exp(g(x, b));
    }

    RowVector grad_f(double x, RowVector b)
    {
        double A = b(0);
        double mu = b(1);
        double s = b(2);

        double f_A = f(x, b)/A;
        double f_mu = f(x, b) * g_mu(x, b);
        double f_s = f(x, b) * g_s(x, b);

        RowVector ret_v(3);
        ret_v(0) = f_A;
        ret_v(1) = f_mu;
        ret_v(2) = f_s;

        return ret_v;
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
