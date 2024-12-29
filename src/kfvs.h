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

#endif