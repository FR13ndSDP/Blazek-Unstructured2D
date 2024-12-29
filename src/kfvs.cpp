/// @file kfvs.cpp
///
/// Computation of 1st-order upwind dissipation.
///
//*****************************************************************************
//
//  (c) J. Blazek, CFD Consulting & Analysis, www.cfd-ca.de
//  Created February 15, 2014
//  Last modification: July 7, 2014
// 
//=============================================================================
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//*****************************************************************************

#include "spaceDiscr.h"
#include "kfvs.h"
#include <iostream>

/// Computes KFVS flux 1st order
///
/// @param geometry    geometrical data
/// @param fluidProps  fluid properties
///
void SpaceDiscr::FluxKFVS1st( const Geometry &geometry, const FluidProps &fluidProps,
                              REAL beta )
{
  int  i, j, ie;
  REAL ds, nx, ny, rrho, rl, ul, vl, pl, rr, ur, vr, pr, beta5;
  REAL fd[4];

  CONSVARS       *cv   = fluidProps.cv;
  DEPVARS        *dv   = fluidProps.dv;
  NODE           *sij  = geometry.sij;
  Geometry::EDGE *edge = geometry.edge;

  double *primL = new double [4];
  double *primR = new double [4];

  beta5 = 0.5*beta;

  for (ie=0; ie<geometry.nEdges; ie++)
  {
    i  = edge[ie].i;
    j  = edge[ie].j;
    ds = SQRT(sij[ie].x*sij[ie].x+sij[ie].y*sij[ie].y);
    nx = sij[ie].x/ds;
    ny = sij[ie].y/ds;

    // left & right state

    rrho = 1.0/cv[i].dens;
    rl   = cv[i].dens;
    ul   = cv[i].xmom*rrho;
    vl   = cv[i].ymom*rrho;
    pl   = dv[i].press;

    rrho = 1.0/cv[j].dens;
    rr   = cv[j].dens;
    ur   = cv[j].xmom*rrho;
    vr   = cv[j].ymom*rrho;
    pr   = dv[j].press;

    tensor1 slope = {1, 0, 0, 0};

    tensor1 MuL (7), MuR (7);
    tensor1 MvL (7), MvR (7);
    tensor1 MxiL(3), MxiR(3);

    // first order part
    primL[0] = rl;
    primL[1] = ul*nx + vl*ny;
    primL[2] = -ul*ny + vl*nx;
    primL[3] = rl/(2*pl);

    primR[0] = rr;
    primR[1] = ur*nx + vr*ny;
    primR[2] = -ur*ny + vr*nx;
    primR[3] = rr/(2*pr);

    cal_Moment_uPositive(primL, MuL, MxiL);
    cal_Moment_uNegative(primR, MuR, MxiR);
    cal_Moment_v(primL, MvL, MxiL);
    cal_Moment_v(primR, MvR, MxiR);

    tensor1 flux(4);

    flux = ( rl * Moment_half(slope, 1, 0, 0, MuL, MvL, MxiL) + rr * Moment_half(slope, 1, 0, 0, MuR, MvR, MxiR));

    fd[0] = flux[0];
    fd[1] = flux[1]*nx-flux[2]*ny;
    fd[2] = flux[1]*ny+flux[2]*nx;
    fd[3] = flux[3];


    ds          *= beta5;
    rhs[i].dens += fd[0]*ds;
    rhs[i].xmom += fd[1]*ds;
    rhs[i].ymom += fd[2]*ds;
    rhs[i].ener += fd[3]*ds;

    rhs[j].dens -= fd[0]*ds;
    rhs[j].xmom -= fd[1]*ds;
    rhs[j].ymom -= fd[2]*ds;
    rhs[j].ener -= fd[3]*ds;
  }

  FluxWalls( geometry,fluidProps );

  delete [] primL;
  delete [] primR;
}

/// Computes KFVS flux 2nd order
///
/// @param geometry    geometrical data
/// @param fluidProps  fluid properties
///
void SpaceDiscr::FluxKFVS2nd( const Geometry &geometry, const FluidProps &fluidProps,
                              REAL beta )
{
  int  i, j, ie;
  REAL ds, nx, ny, rx, ry, rrho, rl, ul, vl, pl, rr, ur, vr, pr, beta5;
  REAL fd[4];

  CONSVARS       *cv   = fluidProps.cv;
  DEPVARS        *dv   = fluidProps.dv;
  NODE           *sij  = geometry.sij;
  NODE           *coords = geometry.coords;
  Geometry::EDGE *edge = geometry.edge;

  double *primL = new double [4];
  double *primR = new double [4];

  beta5 = 0.5*beta;

  for (ie=0; ie<geometry.nEdges; ie++)
  {
    i  = edge[ie].i;
    j  = edge[ie].j;
    ds = SQRT(sij[ie].x*sij[ie].x+sij[ie].y*sij[ie].y);
    nx = sij[ie].x/ds;
    ny = sij[ie].y/ds;
    rx = 0.5*(coords[j].x-coords[i].x);
    ry = 0.5*(coords[j].y-coords[i].y);

    // left & right state

    rrho = 1.0/cv[i].dens;
    rl   = cv[i].dens      + lim[i].dens *(gradx[i].dens *rx+grady[i].dens *ry);
    ul   = cv[i].xmom*rrho + lim[i].uvel *(gradx[i].uvel *rx+grady[i].uvel *ry);
    vl   = cv[i].ymom*rrho + lim[i].vvel *(gradx[i].vvel *rx+grady[i].vvel *ry);
    pl   = dv[i].press     + lim[i].press*(gradx[i].press*rx+grady[i].press*ry);

    rrho = 1.0/cv[j].dens;
    rr   = cv[j].dens      - lim[j].dens *(gradx[j].dens *rx+grady[j].dens *ry);
    ur   = cv[j].xmom*rrho - lim[j].uvel *(gradx[j].uvel *rx+grady[j].uvel *ry);
    vr   = cv[j].ymom*rrho - lim[j].vvel *(gradx[j].vvel *rx+grady[j].vvel *ry);
    pr   = dv[j].press     - lim[j].press*(gradx[j].press*rx+grady[j].press*ry);

    tensor1 slope = {1, 0, 0, 0};

    tensor1 MuL (7), MuR (7);
    tensor1 MvL (7), MvR (7);
    tensor1 MxiL(3), MxiR(3);

    // first order part
    primL[0] = rl;
    primL[1] = ul*nx + vl*ny;
    primL[2] = -ul*ny + vl*nx;
    primL[3] = rl/(2*pl);

    primR[0] = rr;
    primR[1] = ur*nx + vr*ny;
    primR[2] = -ur*ny + vr*nx;
    primR[3] = rr/(2*pr);

    cal_Moment_uPositive(primL, MuL, MxiL);
    cal_Moment_uNegative(primR, MuR, MxiR);
    cal_Moment_v(primL, MvL, MxiL);
    cal_Moment_v(primR, MvR, MxiR);

    tensor1 flux(4);

    flux = ( rl * Moment_half(slope, 1, 0, 0, MuL, MvL, MxiL) + rr * Moment_half(slope, 1, 0, 0, MuR, MvR, MxiR));

    fd[0] = flux[0];
    fd[1] = flux[1]*nx-flux[2]*ny;
    fd[2] = flux[1]*ny+flux[2]*nx;
    fd[3] = flux[3];

    ds          *= beta5;
    rhs[i].dens += fd[0]*ds;
    rhs[i].xmom += fd[1]*ds;
    rhs[i].ymom += fd[2]*ds;
    rhs[i].ener += fd[3]*ds;

    rhs[j].dens -= fd[0]*ds;
    rhs[j].xmom -= fd[1]*ds;
    rhs[j].ymom -= fd[2]*ds;
    rhs[j].ener -= fd[3]*ds;
  }

  FluxWalls( geometry,fluidProps );

  delete [] primL;
  delete [] primR;
}