#include "lambert_gooding.h"

void tlamb(double m, double q, double qsqfm1, double x, double n, double *t, double *dt, double *d2t, double *d3t)
{
    // Gooding support routine
    double sw = 0.4;
    double tv = 0;
    double dtv = 0;
    double d2tv = 0;
    double d3tv = 0.0;

    bool lm1 = (n == -1);
    bool l1 = (n >= 1);
    bool l2 = (n >= 2);
    bool l3 = (n == 3);
    double qsq = q * q;
    double xsq = x * x;
    double u = (1.0 - x) * (1.0 + x);

    double y = 0, z = 0, qx = 0, a = 0, b = 0, aa = 0, bb = 0, g = 0, f = 0, fg1 = 0, term = 0, fg1sq = 0, twoi1 = 0, told = 0, qz = 0, qz2 = 0, u0i = 0, u1i = 0, u2i = 0, u3i = 0, tq = 0, i = 0, tqsum = 0, ttmold = 0, p = 0, tterm = 0, tqterm = 0;

    if (!lm1)
    {
        // (needed if series, and otherwise useful when z = 0)
        dtv = 0.0;
        d2tv = 0.0;
        d3tv = 0.0;
    }
    if (lm1 || m > 0 || x < 0.0 || fabs(u) > sw)
    {
        // direct computation (not series)
        y = sqrt(fabs(u));
        z = sqrt(qsqfm1 + qsq * xsq);
        qx = q * x;
        if (qx <= 0.0)
        {
            a = z - qx;
            b = q * z - x;
        }
        if (qx < 0.0 && lm1)
        {
            aa = qsqfm1 / a;
            bb = qsqfm1 * (qsq * u - xsq) / b;
        }
        if (qx == 0.0 && lm1 || qx > 0.0)
        {
            aa = z + qx;
            bb = q * z + x;
        }
        if (qx > 0.0)
        {
            a = qsqfm1 / aa;
            b = qsqfm1 * (qsq * u - xsq) / bb;
        }
        if (!lm1)
        {
            if (qx * u >= 0.0)
            {
                g = x * z + q * u;
            }
            else
            {
                g = (xsq - qsq * u) / (x * z - q * u);
            }
            f = a * y;
            if (x <= 1.0)
            {
                tv = m * M_PI + atan2(f, g);
            }
            else
            {
                if (f > sw)
                {
                    tv = log(f + g);
                }
                else
                {
                    fg1 = f / (g + 1.0);
                    term = 2.0 * fg1;
                    fg1sq = fg1 * fg1;
                    tv = term;
                    twoi1 = 1.0;

                    while (1)
                    {
                        twoi1 = twoi1 + 2.0;
                        term = term * fg1sq;
                        told = tv;
                        tv = tv + term / twoi1;
                        if (tv != told)
                        {
                            // cycle
                            break;
                        }
                    } // (continue looping for inverse tanh)
                }
            }
            tv = 2.0 * (tv / y + b) / u;
            if (l1 && z != 0.0)
            {
                qz = q / z;
                qz2 = qz * qz;
                qz = qz * qz2;
                dtv = (3.0 * x * tv - 4.0 * (a + qx * qsqfm1) / z) / u;
                if (l2)
                {
                    d2tv = (3.0 * tv + 5.0 * x * dtv + 4.0 * qz * qsqfm1) / u;
                }
                if (l3)
                {
                    d3tv = (8.0 * dtv + 7.0 * x * d2tv - 12.0 * qz * qz2 * x * qsqfm1) / u;
                }
            }

        }
        else{
                dtv = b;
                d2tv = bb;
                d3tv = aa;
        }



    }
    else
        {
            // compute by series
            u0i = 1.0;
            if (l1)
            {
                u1i = 1.0;
            }
            if (l2)
            {
                u2i = 1.0;
            }
            if (l3)
            {
                u3i = 1.0;
            }
            term = 4.0;
            tq = q * qsqfm1;
            i = 0;
            if (q < 0.5)
            {
                tqsum = 1.0 - q * qsq;
            }
            if (q >= 0.5)
            {
                tqsum = (1.0 / (1.0 + q) + q) * qsqfm1;
            }
            ttmold = term / 3.0;
            tv = ttmold * tqsum;
            while (1)
            {
                i = i + 1;
                p = i;
                u0i = u0i * u;
                if (l1 && i > 1)
                {
                    u1i = u1i * u;
                }
                if (l2 && i > 2)
                {
                    u2i = u2i * u;
                }
                if (l3 && i > 3)
                {
                    u3i = u3i * u;
                }
                term = term * (p - 0.5) / p;
                tq = tq * qsq;
                tqsum = tqsum + tq;
                told = tv;
                tterm = term / (2.0 * p + 3.0);
                tqterm = tterm * tqsum;
                tv = tv - u0i * ((1.5 * p + 0.25) * tqterm / (p * p - 0.25) - ttmold * tq);
                ttmold = tterm;
                tqterm = tqterm * p;
                if (l1)
                {
                    dtv = dtv + tqterm * u1i;
                }
                if (l2)
                {
                    d2tv = d2tv + tqterm * u2i * (p - 1.0);
                }
                if (l3)
                {
                    d3tv = d3tv + tqterm * u3i * (p - 1.0) * (p - 2.0);
                }
                if (i < n || tv != told)
                {
                    // cycle
                    break;
                }
            }
            if (l3)
            {
                d3tv = 8.0 * x * (1.5 * d2tv - xsq * d3tv);
            }
            if (l2)
            {
                d2tv = 2.0 * (2.0 * xsq * d2tv - dtv);
            }
            if (l1)
            {
                dtv = -2.0 * x * dtv;
            }
            tv = tv / xsq;
        }

        *t = tv;
        *dt = dtv;
        *d2t = d2tv;
        *d3t = d3tv;
}