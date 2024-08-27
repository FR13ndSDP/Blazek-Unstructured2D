#ifndef KFVS
#define KFVS

#include <valarray>

using label   = long long int;
using scalar  = double;
using tensor1 = std::valarray<double>;
using tensor2 = std::valarray<std::valarray<double>>;

void cal_Moment_xi(const double &lambda, tensor1 &Mxi)
{
    Mxi[0] = 1.0;                                     // <\xi^0>
    Mxi[1] = 3/(2.0*lambda);                          // <\xi^2>
    Mxi[2] = (3.0*3 + 3*(3-1.0))/(4*lambda*lambda);   // <\xi^4>
}

void cal_Moment_uPositive(const double *prim, tensor1 &MuL, tensor1 &Mxi)// x<0, left
{
    
    const double U = prim[1], lambda = prim[3];

    MuL[0] = 0.5 * erfc( -sqrt(lambda)* U ) ;
    MuL[1] = U * MuL[0] + 0.5*exp( -lambda*U*U ) / sqrt( M_PI*lambda );

    for (std::size_t i = 2; i < (MuL).size(); ++i)
    {
        MuL[i] = U * MuL[i - 1] + 0.5 * (i - 1) * MuL[i - 2] / lambda;
    }

    cal_Moment_xi(lambda, Mxi);
}

void cal_Moment_uNegative(const double *prim, tensor1 &MuR, tensor1 &Mxi)// x>0, right
{

    const double U = prim[1], lambda = prim[3];

    MuR[0] = 0.5 * erfc( sqrt(lambda)* U ) ;
    MuR[1] = U * MuR[0] - 0.5*exp( -lambda*U*U ) / sqrt( M_PI*lambda );

    for (std::size_t i = 2; i < (MuR).size(); ++i)
    {
        MuR[i] = U * MuR[i - 1] + 0.5 * (i - 1) * MuR[i - 2] / lambda;
    }
   
    cal_Moment_xi(lambda, Mxi);
}

void cal_Moment_v(const double *prim, tensor1 &Mu, tensor1 &Mxi)
{
    
    const double U = prim[2], lambda = prim[3];
    
    Mu[0] = 1.0;                                     // <u^0>
    Mu[1] = U;                                       // <u^1>

    for (std::size_t i = 2; i < (Mu).size(); ++i)
    {
        Mu[i] = U * Mu[i - 1] + 0.5 * (i - 1) * Mu[i - 2] / lambda;
    }

    cal_Moment_xi(lambda, Mxi);
}

tensor1 Moment_half(const tensor1 &slope, label m, label n, label k, const tensor1 &Mu, const tensor1 &Mv, const tensor1 &Mxi)
{
    tensor1 moment(4);

    moment[0] = slope[0] * Mu[m] * Mv[k] * Mxi[n] \
                + slope[1] * Mu[m+1] * Mv[k] * Mxi[n] \
                + slope[2] * Mu[m] * Mv[k+1] * Mxi[n] \
                + slope[3] * 0.5 * (Mu[m+2] * Mv[k] * Mxi[n] + Mu[m] * Mv[k+2] * Mxi[n] + Mu[m] * Mv[k] * Mxi[n+1]);

    moment[1] = slope[0] * Mu[m+1] * Mv[k] * Mxi[n] \
                + slope[1] * Mu[m+2] * Mv[k] * Mxi[n] \
                + slope[2] * Mu[m+1] * Mv[k+1] * Mxi[n] \
                + slope[3] * 0.5 * (Mu[m+3] * Mv[k] * Mxi[n] + Mu[m+1] * Mv[k+2] * Mxi[n] + Mu[m+1] * Mv[k] * Mxi[n+1]);

    moment[2] = slope[0] * Mu[m] * Mv[k+1] * Mxi[n] \
                + slope[1] * Mu[m+1] * Mv[k+1] * Mxi[n] \
                + slope[2] * Mu[m] * Mv[k+2] * Mxi[n] \
                + slope[3] * 0.5 * (Mu[m+2] * Mv[k+1] * Mxi[n] + Mu[m] * Mv[k+3] * Mxi[n] + Mu[m] * Mv[k+1] * Mxi[n+1]);

    moment[3] = 0.5 * (slope[0] * (Mu[m+2] * Mv[k] * Mxi[n] + Mu[m] * Mv[k+2] * Mxi[n] + Mu[m] * Mv[k] * Mxi[n+1]) \
                + slope[1] * (Mu[m+3] * Mv[k] * Mxi[n] + Mu[m+1] * Mv[k+2] * Mxi[n] + Mu[m+1] * Mv[k] * Mxi[n+1]) \
                + slope[2] * (Mu[m+2] * Mv[k+1] * Mxi[n] + Mu[m] * Mv[k+3] * Mxi[n] + Mu[m] * Mv[k+1] * Mxi[n+1]) \
                + slope[3] * 0.5 * (Mu[m+4] * Mv[k] * Mxi[n] + Mu[m] * Mv[k+4] * Mxi[n] + Mu[m] * Mv[k] * Mxi[n+2] \
                + 2.0 * (Mu[m+2] * Mv[k+2] * Mxi[n] + Mu[m+2] * Mv[k] * Mxi[n+1] + Mu[m] * Mv[k+2] * Mxi[n+1])));

    return moment;
}

void
riemann (const double gamma, const double smallp, const double smallr,
         const double rl, const double ul, const double pl,
         const double ut1l,
         const double rr, const double ur, const double pr,
         const double ut1r, 
         double& flxrho, double& flxu, double& flxut,
         double& flxe) noexcept
{
    constexpr double weakwv = double(1.e-3);
    constexpr double small = double(1.e-6);

    double clsql = gamma*pl*rl;
    double clsqr = gamma*pr*rr;
    double wl = std::sqrt(clsql);
    double wr = std::sqrt(clsqr);
    double cleft = wl/rl;
    double cright = wr/rr;
    double ccsmall = small*(cleft+cright);

    double pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = std::max(pstar,smallp);
    double pstnm1 = pstar;

    double wlsq = (double(0.5)*(gamma-double(1.))*(pstar+pl)+pstar)*rl;
    double wrsq = (double(0.5)*(gamma-double(1.))*(pstar+pr)+pstar)*rr;

    wl = std::sqrt(wlsq);
    wr = std::sqrt(wrsq);
    double ustarp = ul - (pstar-pl)/wl;
    double ustarm = ur + (pstar-pr)/wr;

    pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
    pstar = std::max(pstar,smallp);

    double ustar;
    for (int iter = 0; iter < 3; ++iter)
    {
        wlsq = (double(0.5)*(gamma-double(1.))*(pstar+pl)+pstar)*rl;
        wrsq = (double(0.5)*(gamma-double(1.))*(pstar+pr)+pstar)*rr;

        wl = double(1.)/std::sqrt(wlsq);
        wr = double(1.)/std::sqrt(wrsq);

        double ustnm1 = ustarm;
        double ustnp1 = ustarp;

        ustarm = ur - (pr - pstar)*wr;
        ustarp = ul + (pl - pstar)*wl;

        double dpditer = abs(pstnm1-pstar);
        double zp = abs(ustarp-ustnp1);
        if (zp-weakwv*cleft < double(0.0) ) {
            zp = dpditer*wl;
        }
        double zm = abs(ustarm-ustnm1);
        if (zm-weakwv*cright < double(0.0) ) {
            zm = dpditer*wr;
        }

        double zz = zp+zm;
        double denom = dpditer/ std::max(zz,ccsmall);
        pstnm1 = pstar;
        pstar = pstar - denom*(ustarm-ustarp);
        pstar = std::max(pstar,smallp);
        ustar = double(0.5)*(ustarm+ustarp);
    }

    double ro, uo, po, sgnm, utrans1;
    if (ustar > double(0.)) {
        ro = rl;
        uo = ul;
        po = pl;
        sgnm = double(1.);
        utrans1 = ut1l;
    } else if (ustar < double(0.)) {
        ro = rr;
        uo = ur;
        po = pr;
        sgnm = double(-1.);
        utrans1 = ut1r;
    } else {
        uo = double(0.5)*(ur+ul);
        po = double(0.5)*(pr+pl);
        ro = double(2.)*(rl*rr)/(rl+rr);
        sgnm = double(1.);
        utrans1 = double(0.5)*(ut1l+ut1r);
    }
    double wosq = (double(0.5)*(gamma-double(1.))*(pstar+po)+pstar)*ro;
    double co = std::sqrt(gamma * po / ro);
    double wo = std::sqrt(wosq);
    double dpjmp = pstar-po;
    double rstar = ro/(double(1.)-ro*dpjmp/wosq);
    double cstar = std::sqrt(gamma * pstar / rstar);
    double spout = co-sgnm*uo;
    double spin = cstar - sgnm*uo;
    if(pstar >= po) {
        spin = wo/ro-sgnm*uo;
        spout = spin;
    }
    double ss = std::max(spout-spin, spout+spin);
    double frac = double(0.5)*(double(1.)+(spin+spout)/std::max(ss,ccsmall));

    double rgdnv, ugdnv, pgdnv;
    if (spout < double(0.)) {
        rgdnv = ro;
        ugdnv = uo;
        pgdnv = po;
    } else if(spin >= double(0.)) {
        rgdnv = rstar;
        ugdnv = ustar;
        pgdnv = pstar;
    } else {
        rgdnv = frac*rstar + (double(1.) - frac)* ro;
        ugdnv = frac*ustar + (double(1.) - frac)* uo;
        pgdnv = frac*pstar + (double(1.) - frac)* po;
    }

    flxrho = rgdnv*ugdnv;
    flxu = rgdnv*ugdnv*ugdnv+pgdnv;
    flxut = rgdnv*ugdnv*utrans1;
    flxe = ugdnv*(double(0.5)*rgdnv*(ugdnv*ugdnv+utrans1*utrans1) + pgdnv/(gamma -double(1.)) + pgdnv);
}

#endif