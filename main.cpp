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
        const int max_iter = 10;

        for (int i = 0; i < max_iter; ++i)
        {
            Matrix J = Jf(b);
            Matrix Jt = J.transpose();

            Matrix T = Jt*J;
            RowVector bn = T.solve(Jt*res(b));

            b = b + bn;
        }

        return b;


    }


};


    // int val;
    // while (cin >> val)
    // {
    //     cout << val << endl;
    //
    // };
    //
    // return 0;



int main(int argc, const char * argv[])
{
    using namespace std;

    GaussNewton gn;
    gn.dim_data = 3;


    gn.pts.push_back(GaussNewton::Point(-1, 5));
    gn.pts.push_back(GaussNewton::Point(0, 10));
    gn.pts.push_back(GaussNewton::Point(1, 5));

    RowVector b(3);
    b(0) = 1;
    b(1) = 0;
    b(2) = 1;

    RowVector fb = gn.run(b);

    for (double x = -5; x < 5; x += 0.05)
    {
        double vf = gn.f(x, fb);

        printf("%.3lf, %.3lf\n",x, vf);
    }

    return 0;
}
