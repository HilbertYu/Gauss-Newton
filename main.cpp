#include <stdio.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/builtin-defun-decls.h>
#include <octave/CMatrix.h>
#include <octave/Array.h>


//vb = A, mu, sigma

//template <int dim_b, int dim_data>
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
    class Point
    {
    public:
        double x;
        double y;
        Point(double _x, double _y):
            x(_x),
            y(_y)
        {}

    };

    int dim_data;


    std::vector<GaussNewton::Point> pts;

    double f(double x, RowVector b)
    {
        return b(0) * exp(g(x, b));
    }

    ColumnVector res(RowVector b)
    {
        const int N = pts.size();
        ColumnVector ret(N);

        for (int i = 0; i < N; ++i)
        {
            double x = pts[i].x;
            double y = pts[i].y;
            ret(i) = y - f(x, b);
        }
        return ret;

    }

    RowVector grad_f(double x, RowVector b)
    {
        double A = b(0);
        //double mu = b(1);
       // double s = b(2);

        double f_A = f(x, b)/A;
        double f_mu = f(x, b) * g_mu(x, b);
        double f_s = f(x, b) * g_s(x, b);

        RowVector ret_v(3);
        ret_v(0) = f_A;
        ret_v(1) = f_mu;
        ret_v(2) = f_s;

        return ret_v;
    }

    Matrix Jf(RowVector b)
    {
        const int N = pts.size();
        Matrix ret_f(N, 3);

        for (int ir = 0; ir < N; ++ ir)
        {
            for (int j = 0; j < 3; ++j)
            {
                ret_f(ir, j) = grad_f(pts[ir].x, b)(j);
            }

        }
        return ret_f;

    }

    RowVector run(RowVector b)
    {
        const int max_iter = 20;

        for (int i = 0; i < max_iter; ++i)
        {
            Matrix J = Jf(b);
            Matrix Jt = J.transpose();

            Matrix T = Jt*J;
            RowVector bn = T.solve(Jt*res(b));

            double err = sqrt(bn(0)*bn(0) + bn(1)*bn(1) + bn(2)*bn(2));
        //    printf("err = %.5lf\n", err);

            b = b + bn;
        }

        return b;


    }


};





int main(int argc, const char * argv[])
{
    using namespace std;

    GaussNewton gn;


    int val;
    int ix = 0;
    while (cin >> val)
    {
       // cout << val << endl;
        gn.pts.push_back(GaussNewton::Point(ix, val));
        ++ix;

    };

    gn.dim_data = gn.pts.size();

    RowVector b(3);
    b(0) = 1500;
    b(1) = 90;
    b(2) = 1;

    RowVector fb = gn.run(b);

    for (auto pt: gn.pts)
    {
        double vf = gn.f(pt.x, fb);

        printf("%.3lf, %.3lf\n",pt.x, vf);
    }

    return 0;
}
