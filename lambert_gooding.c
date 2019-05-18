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
                            break;;
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
                    break;;
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

double d8rt(double x){
    return sqrt(sqrt(sqrt(x)));
}

void xlamb(double m,double q,double qsqfm1,double tin, double * n,double * x,double * xpl){
    double tol = 3e-7;
    double c0  = 1.7;
    double c1  = 0.5;
    double c2  = 0.03;
    double c3  = 0.15;
    double c41 = 1.0;
    double c42 = 0.24;

    double t0=0,dt=0,d2t=0,d3t=0,d2t2=0;

    double thr2 = atan2(qsqfm1, 2.0*q)/M_PI;
    double tdiff=0,w=0,xm=0,xtest=0,xmold=0,tdiffm=0,tmin=0,t=0,tdiff0=0,ij=0;

    (*xpl) = 0;
    (*x) = 0;

    if (m==0){
        // single-rev starter from t (at (*x) = 0) & bilinear (usually)

        (*n) = 1;
        tlamb(m,q,qsqfm1,0,0,&t0,&dt,&d2t,&d3t);
        tdiff = tin - t0;
        if (tdiff<=0.0){
            (*x) = t0*tdiff/(-4.0*tin);
            // (-4 is the value of dt, for (*x) = 0)
        }
        else{
            (*x) = -tdiff/(tdiff + 4.0);
            w = (*x) + c0*sqrt(2.0*(1.0 - thr2));
            if (w<0.0){
                (*x) = (*x) - sqrt(d8rt(-w))*((*x) + sqrt(tdiff/(tdiff + 1.5*t0)));
            }
            w = 4.0/(4.0 + tdiff);
            (*x) = (*x)*(1.0 + (*x)*(c1*w - c2*(*x)*sqrt(w)));
        }
    }else{
        // with multirevs, first get t(min) as basis for starter
        xm = 1.0/(1.5*(m + 0.5)*M_PI);
        if (thr2<0.5){
            xm = d8rt(2.0*thr2)*xm;
        }
        if (thr2>0.5){
            xm = (2.0 - d8rt(2.0 - 2.0*thr2))*xm;
        }
        // (starter for tmin)

        int i=1;
        for (i;i<=12;i++){
            tlamb(m,q,qsqfm1,xm,3,&tmin,&dt,&d2t,&d3t);
            
            if (d2t==0.0){
                break;
            }
            xmold = xm;
            xm = xm - dt*d2t/(d2t*d2t - dt*d3t/2.0);
            xtest = abs(xmold/xm - 1.0);
            if (xtest<=tol){
                break;
            }
        }
        
        if (i>12){
            // (break; off & exit if tmin not located - should never happen)
            // now proceed from t(min) to full starter
            (*n) = -1;
            return;
        }
        tdiffm = tin - tmin;
        if (tdiffm<0.0){
            (*n) = 0;
            return;
            // (exit if no solution with this m)
        }else if (tdiffm==0.0){
            (*x) = xm;
            (*n) = 1;
            return;
            // (exit if unique solution already from (*x)(tmin))
        }else{
            (*n) = 3;
            if (d2t==0.0){
                d2t = 6.0*m*M_PI;
            }
            (*x) = sqrt(tdiffm/(d2t/2.0 + tdiffm/((1.0 - xm)*(1.0 - xm))));
            w = xm + (*x);
            w = w*4.0/(4.0 + tdiffm) + (1.0 - w)*(1.0 - w);
            (*x) = (*x)*(1.0 - (1.0 + m + c41*(thr2 - 0.5))/(1.0 + c3*m)*(*x)*(c1*w + c2*(*x)*sqrt(w))) + xm;
            d2t2 = d2t/2.0;
            if ((*x)>=1.0){
                (*n) = 1;
                // goto 3
            tlamb(m,q,qsqfm1,0.0,0,&t0,&dt,&d2t,&d3t);
            tdiff0 = t0 - tmin;
            tdiff = tin - t0;
            if (tdiff<=0){
                (*x) = xm - sqrt(tdiffm/(d2t2 - tdiffm*(d2t2/tdiff0 - 1.0/(xm*xm))));
            }else{
                (*x) = -tdiff/(tdiff + 4.0);
                ij = 200;
                w = (*x) + c0*sqrt(2.0*(1.0 - thr2));
                if (w<0.0){
                    (*x) = (*x) - sqrt(d8rt(-w))*((*x) + sqrt(tdiff/(tdiff+1.5*t0)));
                }
                w = 4.0/(4.0 + tdiff);
                (*x) = (*x)*(1.0 + (1.0 + m + c42*(thr2 - 0.5))/(1.0 + c3*m)*(*x)*(c1*w - c2*(*x)*sqrt(w)));
                if ((*x)<=-1.0){
                    (*n) = (*n) - 1;
                    // (no finite solution with (*x) < xm)
                    if ((*n)==1){
                        (*x) = (*xpl);
                    }
                }
            } // 3
            }
            // (no finite solution with (*x) > xm)
        }
    }
        // (now have a starter, so proceed by halley)
    // while(1)
        for (int i=1;i<=3;i++){
            
            tlamb(m,q,qsqfm1,(*x),2,&t,&dt,&d2t,&d3t);
            t = tin - t;
            if (dt!=0.0){
                (*x) = (*x) + t*dt/(dt*dt + t*d2t/2.0);
            }
        }
        if ((*n)!=3){
            return;
        }
        // (exit if only one solution, normally when m = 0)
        (*n) = 2;
        (*xpl) = (*x);
        // (second multi-rev starter)
        
        tlamb(m,q,qsqfm1,0.0,0,&t0,&dt,&d2t,&d3t); //3
        tdiff0 = t0 - tmin;
        tdiff = tin - t0;
        if (tdiff<=0){
            (*x) = xm - sqrt(tdiffm/(d2t2 - tdiffm*(d2t2/tdiff0 - 1.0/(xm*xm))));
        }else{
            (*x) = -tdiff/(tdiff + 4.0);
            ij = 200;
            w = (*x) + c0*sqrt(2.0*(1.0 - thr2));
            if (w<0.0){
                (*x) = (*x) - sqrt(d8rt(-w))*((*x) + sqrt(tdiff/(tdiff+1.5*t0)));
            }
            w = 4.0/(4.0 + tdiff);
            (*x) = (*x)*(1.0 + (1.0 + m + c42*(thr2 - 0.5))/(1.0 + c3*m)*(*x)*(c1*w - c2*(*x)*sqrt(w)));
            if ((*x)<=-1.0){
                (*n) = (*n) - 1;
                // (no finite solution with (*x) < xm)
                if ((*n)==1){
                    (*x) = (*xpl);
                }
            }
        }
}