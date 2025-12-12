//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB SP code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

#include "header.h"

//---------------------------------------------------------------------
// this function performs the solution of the approximate factorization
// step in the x-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the x-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void solve()
{
  int i, j, k, i1, i2, j1, j2, k1, k2, m;
  double ru1, fac1, fac2;
  double r1, r2, r3, r4, r5, t1, t2, t3;
  double btuz, ac2u, uzik1, xvel, yvel, zvel, ac;

  if (timeron) timer_start(t_xsolve);
  #pragma dvm actual(u,rhs,us,vs,ws,qs,rho_i,speed,square)
  #pragma dvm interval

  #pragma dvm region local(lhs, lhsp, lhsm) inout(rhs)
  {
    #pragma dvm parallel([k][j] on rhs[k][j][][]) \
            private(i,j,k,m,i1,i2,ru1,fac1,fac2)
    for (k = 1; k <= nz2; k++) {
      for (j = 1; j <= ny2; j++) {
        //---------------------------------------------------------------------
        //location part
        //---------------------------------------------------------------------
        double cv_loc[IMAXP+1];
        double rhon_loc[IMAXP+1];

        int ni = nx2 + 1; 
        for (m = 0; m < 5; m++) {
          lhs [k][j][0 ][m] = 0.0;
          lhsp[k][j][0 ][m] = 0.0;
          lhsm[k][j][0 ][m] = 0.0;
          lhs [k][j][ni][m] = 0.0;
          lhsp[k][j][ni][m] = 0.0;
          lhsm[k][j][ni][m] = 0.0;
        }
        lhs [k][j][0 ][2] = 1.0;
        lhsp[k][j][0 ][2] = 1.0;
        lhsm[k][j][0 ][2] = 1.0;
        lhs [k][j][ni][2] = 1.0;
        lhsp[k][j][ni][2] = 1.0;
        lhsm[k][j][ni][2] = 1.0;
        //---------------------------------------------------------------------
        // Computes the left hand side for the three x-factors  
        //---------------------------------------------------------------------

        //---------------------------------------------------------------------
        // first fill the lhs for the u-eigenvalue                   
        //---------------------------------------------------------------------
        for (i = 0; i <= grid_points[0]-1; i++) {
          ru1 = c3c4*rho_i[k][j][i];
          cv_loc[i] = us[k][j][i];
          rhon_loc[i] = max(max(dx2+con43*ru1,dx5+c1c5*ru1), max(dxmax+ru1,dx1));
        }

        for (i = 1; i <= nx2; i++) {
          lhs[k][j][i][0] =  0.0;
          lhs[k][j][i][1] = -dttx2 * cv_loc[i-1] - dttx1 * rhon_loc[i-1];
          lhs[k][j][i][2] =  1.0 + c2dttx1 * rhon_loc[i];
          lhs[k][j][i][3] =  dttx2 * cv_loc[i+1] - dttx1 * rhon_loc[i+1];
          lhs[k][j][i][4] =  0.0;
        }
        //---------------------------------------------------------------------
        // add fourth order dissipation                             
        //---------------------------------------------------------------------
        i = 1;
        lhs[k][j][i][2] += comz5;
        lhs[k][j][i][3] -= comz4;
        lhs[k][j][i][4] += comz1;
        lhs[k][j][i+1][1] -= comz4;
        lhs[k][j][i+1][2] += comz6;
        lhs[k][j][i+1][3] -= comz4;
        lhs[k][j][i+1][4] += comz1;

        for (i = 3; i <= grid_points[0]-4; i++) {
          lhs[k][j][i][0] += comz1;
          lhs[k][j][i][1] -= comz4;
          lhs[k][j][i][2] += comz6;
          lhs[k][j][i][3] -= comz4;
          lhs[k][j][i][4] += comz1;
        }

        i = grid_points[0]-3;
        lhs[k][j][i][0] += comz1;
        lhs[k][j][i][1] -= comz4;
        lhs[k][j][i][2] += comz6;
        lhs[k][j][i][3] -= comz4;
        lhs[k][j][i+1][0] += comz1;
        lhs[k][j][i+1][1] -= comz4;
        lhs[k][j][i+1][2] += comz5;
        //---------------------------------------------------------------------
        // subsequently, fill the other factors (u+c), (u-c) by adding to 
        // the first  
        //---------------------------------------------------------------------
        for (i = 1; i <= nx2; i++) {
          lhsp[k][j][i][0] = lhs[k][j][i][0];
          lhsp[k][j][i][1] = lhs[k][j][i][1] - dttx2 * speed[k][j][i-1];
          lhsp[k][j][i][2] = lhs[k][j][i][2];
          lhsp[k][j][i][3] = lhs[k][j][i][3] + dttx2 * speed[k][j][i+1];
          lhsp[k][j][i][4] = lhs[k][j][i][4];
          lhsm[k][j][i][0] = lhs[k][j][i][0];
          lhsm[k][j][i][1] = lhs[k][j][i][1] + dttx2 * speed[k][j][i-1];
          lhsm[k][j][i][2] = lhs[k][j][i][2];
          lhsm[k][j][i][3] = lhs[k][j][i][3] - dttx2 * speed[k][j][i+1];
          lhsm[k][j][i][4] = lhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // FORWARD ELIMINATION  
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // perform the Thomas algorithm; first, FORWARD ELIMINATION     
        //---------------------------------------------------------------------
        for (i = 0; i <= grid_points[0]-3; i++) {
          i1 = i + 1;
          i2 = i + 2;
          fac1 = 1.0/lhs[k][j][i][2];
          lhs[k][j][i][3] *= fac1;
          lhs[k][j][i][4] *= fac1;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
          }
          lhs[k][j][i1][2] -= lhs[k][j][i1][1]*lhs[k][j][i][3];
          lhs[k][j][i1][3] -= lhs[k][j][i1][1]*lhs[k][j][i][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j][i1][m] -= lhs[k][j][i1][1]*rhs[k][j][i][m];
          }
          lhs[k][j][i2][1] -= lhs[k][j][i2][0]*lhs[k][j][i][3];
          lhs[k][j][i2][2] -= lhs[k][j][i2][0]*lhs[k][j][i][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j][i2][m] -= lhs[k][j][i2][0]*rhs[k][j][i][m];
          }
          //---------------------------------------------------------------------
          // for (the u+c and the u-c factors                 
          //---------------------------------------------------------------------
          //m = 3;
          fac1 = 1.0/lhsp[k][j][i][2];
          lhsp[k][j][i][3]   *= fac1;
          lhsp[k][j][i][4]   *= fac1;
          rhs[k][j][i][3]  *= fac1;
          lhsp[k][j][i1][2]  -= lhsp[k][j][i1][1]*lhsp[k][j][i][3];
          lhsp[k][j][i1][3]  -= lhsp[k][j][i1][1]*lhsp[k][j][i][4];
          rhs[k][j][i1][3] -= lhsp[k][j][i1][1]*rhs[k][j][i][3];
          lhsp[k][j][i2][1]  -= lhsp[k][j][i2][0]*lhsp[k][j][i][3];
          lhsp[k][j][i2][2]  -= lhsp[k][j][i2][0]*lhsp[k][j][i][4];
          rhs[k][j][i2][3] -= lhsp[k][j][i2][0]*rhs[k][j][i][3];

          //m = 4;
          fac1 = 1.0/lhsm[k][j][i][2];
          lhsm[k][j][i][3]   *= fac1;
          lhsm[k][j][i][4]   *= fac1;
          rhs[k][j][i][4]  *= fac1;
          lhsm[k][j][i1][2]  -= lhsm[k][j][i1][1]*lhsm[k][j][i][3];
          lhsm[k][j][i1][3]  -= lhsm[k][j][i1][1]*lhsm[k][j][i][4];
          rhs[k][j][i1][4] -= lhsm[k][j][i1][1]*rhs[k][j][i][4];
          lhsm[k][j][i2][1]  -= lhsm[k][j][i2][0]*lhsm[k][j][i][3];
          lhsm[k][j][i2][2]  -= lhsm[k][j][i2][0]*lhsm[k][j][i][4];
          rhs[k][j][i2][4] -= lhsm[k][j][i2][0]*rhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        i  = grid_points[0]-2;
        i1 = grid_points[0]-1;
        fac1 = 1.0/lhs[k][j][i][2];
        lhs[k][j][i][3] *= fac1;
        lhs[k][j][i][4] *= fac1;
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[k][j][i1][2] -= lhs[k][j][i1][1]*lhs[k][j][i][3];
        lhs[k][j][i1][3] -= lhs[k][j][i1][1]*lhs[k][j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i1][m] -= lhs[k][j][i1][1]*rhs[k][j][i][m];
        }

        //---------------------------------------------------------------------
        // scale the last row immediately 
        //---------------------------------------------------------------------
        fac2 = 1.0/lhs[k][j][i1][2];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i1][m] = fac2*rhs[k][j][i1][m];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        //m = 3;
        fac1 = 1.0/lhsp[k][j][i][2];
        lhsp[k][j][i][3]   *= fac1;
        lhsp[k][j][i][4]   *= fac1;
        rhs[k][j][i][3]  *= fac1;
        lhsp[k][j][i1][2]  -= lhsp[k][j][i1][1]*lhsp[k][j][i][3];
        lhsp[k][j][i1][3]  -= lhsp[k][j][i1][1]*lhsp[k][j][i][4];
        rhs[k][j][i1][3] -= lhsp[k][j][i1][1]*rhs[k][j][i][3];
        rhs[k][j][i1][3] /= lhsp[k][j][i1][2];

        //m = 4;
        fac1 = 1.0/lhsm[k][j][i][2];
        lhsm[k][j][i][3]   *= fac1;
        lhsm[k][j][i][4]   *= fac1;
        rhs[k][j][i][4]  *= fac1;
        lhsm[k][j][i1][2]  -= lhsm[k][j][i1][1]*lhsm[k][j][i][3];
        lhsm[k][j][i1][3]  -= lhsm[k][j][i1][1]*lhsm[k][j][i][4];
        rhs[k][j][i1][4] -= lhsm[k][j][i1][1]*rhs[k][j][i][4];
        rhs[k][j][i1][4] /= lhsm[k][j][i1][2];
        //---------------------------------------------------------------------
        // BACKSUBSTITUTION 
        //---------------------------------------------------------------------
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] -= lhs[k][j][i][3]*rhs[k][j][i1][m];
        }
        rhs[k][j][i][3] -= lhsp[k][j][i][3]*rhs[k][j][i1][3];
        rhs[k][j][i][4] -= lhsm[k][j][i][3]*rhs[k][j][i1][4];
        //---------------------------------------------------------------------
        // The first three factors
        //---------------------------------------------------------------------
        for (i = grid_points[0]-3; i >= 0; i--) {
          i1 = i + 1;
          i2 = i + 2;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] -= lhs[k][j][i][3]*rhs[k][j][i1][m] +
                               lhs[k][j][i][4]*rhs[k][j][i2][m];
          }
          //-------------------------------------------------------------------
          // And the remaining two
          //-------------------------------------------------------------------
          rhs[k][j][i][3] -= lhsp[k][j][i][3]*rhs[k][j][i1][3] +
                             lhsp[k][j][i][4]*rhs[k][j][i2][3];
          rhs[k][j][i][4] -= lhsm[k][j][i][3]*rhs[k][j][i1][4] +
                             lhsm[k][j][i][4]*rhs[k][j][i2][4];
        }
      }
    }
  } 

  //---------------------------------------------------------------------
  // ninvr        
  //---------------------------------------------------------------------
  #pragma dvm region in(u, us, vs, ws, qs, speed, rhs) out(rhs)
  {
    #pragma dvm parallel([k][j][i] on rhs[k][j][i][]) private(k,j,i,r1,r2,r3,r4,r5,t1,t2)
    for (k = 1; k <= nz2; k++) {
      for (j = 1; j <= ny2; j++) {
        for (i = 1; i <= nx2; i++) {
          r1 = rhs[k][j][i][0];
          r2 = rhs[k][j][i][1];
          r3 = rhs[k][j][i][2];
          r4 = rhs[k][j][i][3];
          r5 = rhs[k][j][i][4];

          t1 = bt * r3;
          t2 = 0.5 * ( r4 + r5 );

          rhs[k][j][i][0] = -r2;
          rhs[k][j][i][1] =  r1;
          rhs[k][j][i][2] = bt * ( r4 - r5 );
          rhs[k][j][i][3] = -t1 + t2;
          rhs[k][j][i][4] =  t1 + t2;
        }
      }
    }
  }

  //---------------------------------------------------------------------
  //y_solve
  //---------------------------------------------------------------------
  #pragma dvm region local(lhs, lhsp, lhsm) inout(rhs)
  {
    #pragma dvm parallel([k][i] on rhs[k][][i][]) private(i,k,j, j1, j2, m, ru1, fac1, fac2)
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (i = 1; i <= grid_points[0]-2; i++) {

        //---------------------------------------------------------------------
        //location part
        //---------------------------------------------------------------------
        double cv_loc[JMAXP+1];
        double rhoq_loc[JMAXP+1];

        int nj = ny2+1;
        for (m = 0; m < 5; m++) {
          lhs [k][0][i][m] = 0.0;
          lhsp[k][0][i][m] = 0.0;
          lhsm[k][0][i][m] = 0.0;
          lhs [k][nj][i][m] = 0.0;
          lhsp[k][nj][i][m] = 0.0;
          lhsm[k][nj][i][m] = 0.0;
        }
        lhs [k][0][i][2] = 1.0;
        lhsp[k][0][i][2] = 1.0;
        lhsm[k][0][i][2] = 1.0;
        lhs [k][nj][i][2] = 1.0;
        lhsp[k][nj][i][2] = 1.0;
        lhsm[k][nj][i][2] = 1.0;

        //---------------------------------------------------------------------
        // Computes the left hand side for the three y-factors   
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // first fill the lhs for the u-eigenvalue         
        //---------------------------------------------------------------------
        for (j = 0; j <= grid_points[1]-1; j++) {
          ru1 = c3c4*rho_i[k][j][i];
          cv_loc[j] = vs[k][j][i];
          rhoq_loc[j] = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
        }

        for (j = 1; j <= grid_points[1]-2; j++) {
          lhs[k][j][i][0] =  0.0;
          lhs[k][j][i][1] = -dtty2 * cv_loc[j-1] - dtty1 * rhoq_loc[j-1];
          lhs[k][j][i][2] =  1.0 + c2dtty1 * rhoq_loc[j];
          lhs[k][j][i][3] =  dtty2 * cv_loc[j+1] - dtty1 * rhoq_loc[j+1];
          lhs[k][j][i][4] =  0.0;
        }
        //---------------------------------------------------------------------
        // add fourth order dissipation                             
        //---------------------------------------------------------------------
        j = 1;
        lhs[k][j][i][2] += comz5;
        lhs[k][j][i][3] -= comz4;
        lhs[k][j][i][4] += comz1;
        lhs[k][j+1][i][1] -= comz4;
        lhs[k][j+1][i][2] += comz6;
        lhs[k][j+1][i][3] -= comz4;
        lhs[k][j+1][i][4] += comz1;
        
        for (j = 3; j <= grid_points[1]-4; j++) {
          lhs[k][j][i][0] += comz1;
          lhs[k][j][i][1] -= comz4;
          lhs[k][j][i][2] += comz6;
          lhs[k][j][i][3] -= comz4;
          lhs[k][j][i][4] += comz1;
        }

        j = grid_points[1]-3;
        lhs[k][j][i][0] += comz1;
        lhs[k][j][i][1] -= comz4;
        lhs[k][j][i][2] += comz6;
        lhs[k][j][i][3] -= comz4;
        lhs[k][j+1][i][0] += comz1;
        lhs[k][j+1][i][1] -= comz4;
        lhs[k][j+1][i][2] += comz5;

        for (j = 1; j <= grid_points[1]-2; j++) {
          lhsp[k][j][i][0] = lhs[k][j][i][0];
          lhsp[k][j][i][1] = lhs[k][j][i][1] - dtty2 * speed[k][j-1][i];
          lhsp[k][j][i][2] = lhs[k][j][i][2];
          lhsp[k][j][i][3] = lhs[k][j][i][3] + dtty2 * speed[k][j+1][i];
          lhsp[k][j][i][4] = lhs[k][j][i][4];
          lhsm[k][j][i][0] = lhs[k][j][i][0];
          lhsm[k][j][i][1] = lhs[k][j][i][1] + dtty2 * speed[k][j-1][i];
          lhsm[k][j][i][2] = lhs[k][j][i][2];
          lhsm[k][j][i][3] = lhs[k][j][i][3] - dtty2 * speed[k][j+1][i];
          lhsm[k][j][i][4] = lhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // FORWARD ELIMINATION  
        //---------------------------------------------------------------------
        for (j = 0; j <= grid_points[1]-3; j++) {
          j1 = j + 1;
          j2 = j + 2;
          fac1 = 1.0/lhs[k][j][i][2];
          lhs[k][j][i][3] *= fac1;
          lhs[k][j][i][4] *= fac1;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
          }
          lhs[k][j1][i][2] -= lhs[k][j1][i][1]*lhs[k][j][i][3];
          lhs[k][j1][i][3] -= lhs[k][j1][i][1]*lhs[k][j][i][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[k][j1][i][1]*rhs[k][j][i][m];
          }
          lhs[k][j2][i][1] -= lhs[k][j2][i][0]*lhs[k][j][i][3];
          lhs[k][j2][i][2] -= lhs[k][j2][i][0]*lhs[k][j][i][4];
          for (m = 0; m < 3; m++) {
            rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhs[k][j2][i][0]*rhs[k][j][i][m];
          }
          //---------------------------------------------------------------------
          // for (the u+c and the u-c factors                 
          //---------------------------------------------------------------------
          //m = 3;
          fac1 = 1.0/lhsp[k][j][i][2];
          lhsp[k][j][i][3]   *= fac1;
          lhsp[k][j][i][4]   *= fac1;
          rhs[k][j][i][3]  *= fac1;
          lhsp[k][j1][i][2]  -= lhsp[k][j1][i][1]*lhsp[k][j][i][3];
          lhsp[k][j1][i][3]  -= lhsp[k][j1][i][1]*lhsp[k][j][i][4];
          rhs[k][j1][i][3] -= lhsp[k][j1][i][1]*rhs[k][j][i][3];
          lhsp[k][j2][i][1]  -= lhsp[k][j2][i][0]*lhsp[k][j][i][3];
          lhsp[k][j2][i][2]  -= lhsp[k][j2][i][0]*lhsp[k][j][i][4];
          rhs[k][j2][i][3] -= lhsp[k][j2][i][0]*rhs[k][j][i][3];

          //m = 4;
          fac1 = 1.0/lhsm[k][j][i][2];
          lhsm[k][j][i][3]   *= fac1;
          lhsm[k][j][i][4]   *= fac1;
          rhs[k][j][i][4]  *= fac1;
          lhsm[k][j1][i][2]  -= lhsm[k][j1][i][1]*lhsm[k][j][i][3];
          lhsm[k][j1][i][3]  -= lhsm[k][j1][i][1]*lhsm[k][j][i][4];
          rhs[k][j1][i][4] -= lhsm[k][j1][i][1]*rhs[k][j][i][4];
          lhsm[k][j2][i][1]  -= lhsm[k][j2][i][0]*lhsm[k][j][i][3];
          lhsm[k][j2][i][2]  -= lhsm[k][j2][i][0]*lhsm[k][j][i][4];
          rhs[k][j2][i][4] -= lhsm[k][j2][i][0]*rhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // The last two rows in this grid block are a bit different, 
        // since they for (not have two more rows available for the
        // elimination of off-diagonal entries
        //---------------------------------------------------------------------
        j  = grid_points[1]-2;
        j1 = grid_points[1]-1;
        fac1 = 1.0/lhs[k][j][i][2];
        lhs[k][j][i][3] *= fac1;
        lhs[k][j][i][4] *= fac1;
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[k][j1][i][2] -= lhs[k][j1][i][1]*lhs[k][j][i][3];
        lhs[k][j1][i][3] -= lhs[k][j1][i][1]*lhs[k][j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[k][j1][i][1]*rhs[k][j][i][m];
        }
        //---------------------------------------------------------------------
        // scale the last row immediately 
        //---------------------------------------------------------------------
        fac2 = 1.0/lhs[k][j1][i][2];
        for (m = 0; m < 3; m++) {
          rhs[k][j1][i][m] = fac2*rhs[k][j1][i][m];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        //m = 3;
        fac1 = 1.0/lhsp[k][j][i][2];
        lhsp[k][j][i][3]   *= fac1;
        lhsp[k][j][i][4]   *= fac1;
        rhs[k][j][i][3]  *= fac1;
        lhsp[k][j1][i][2]  -= lhsp[k][j1][i][1]*lhsp[k][j][i][3];
        lhsp[k][j1][i][3]  -= lhsp[k][j1][i][1]*lhsp[k][j][i][4];
        rhs[k][j1][i][3] -= lhsp[k][j1][i][1]*rhs[k][j][i][3];
        rhs[k][j1][i][3] /= lhsp[k][j1][i][2];

        //m = 4;
        fac1 = 1.0/lhsm[k][j][i][2];
        lhsm[k][j][i][3]   *= fac1;
        lhsm[k][j][i][4]   *= fac1;
        rhs[k][j][i][4]  *= fac1;
        lhsm[k][j1][i][2]  -= lhsm[k][j1][i][1]*lhsm[k][j][i][3];
        lhsm[k][j1][i][3]  -= lhsm[k][j1][i][1]*lhsm[k][j][i][4];
        rhs[k][j1][i][4] -= lhsm[k][j1][i][1]*rhs[k][j][i][4];
        rhs[k][j1][i][4] /= lhsm[k][j1][i][2];

        //---------------------------------------------------------------------
        // BACKSUBSTITUTION 
        //---------------------------------------------------------------------

        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] -= lhs[k][j][i][3]*rhs[k][j1][i][m];
        }

        rhs[k][j][i][3] -= lhsp[k][j][i][3]*rhs[k][j1][i][3];
        rhs[k][j][i][4] -= lhsm[k][j][i][3]*rhs[k][j1][i][4];

        //---------------------------------------------------------------------
        // The first three factors
        //---------------------------------------------------------------------
        for (j = grid_points[1]-3; j >= 0; j--) {
          j1 = j + 1;
          j2 = j + 2;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] -= lhs[k][j][i][3]*rhs[k][j1][i][m] +
                              lhs[k][j][i][4]*rhs[k][j2][i][m];
          }

          //-------------------------------------------------------------------
          // And the remaining two
          //-------------------------------------------------------------------
          rhs[k][j][i][3] -= lhsp[k][j][i][3]*rhs[k][j1][i][3] +
                            lhsp[k][j][i][4]*rhs[k][j2][i][3];
          rhs[k][j][i][4] -= lhsm[k][j][i][3]*rhs[k][j1][i][4] +
                            lhsm[k][j][i][4]*rhs[k][j2][i][4];
        }
      }
    }
  }

  //---------------------------------------------------------------------
  // pinvr                    
  //---------------------------------------------------------------------
  #pragma dvm region in(u, us, vs, ws, qs, speed, rhs) out(rhs)
  {
    #pragma dvm parallel([k][j][i] on rhs[k][j][i][]) private(k,j,i,r1,r2,r3,r4,r5,t1,t2)
    for (k = 1; k <= nz2; k++) {
      for (j = 1; j <= ny2; j++) {
        for (i = 1; i <= nx2; i++) {
          r1 = rhs[k][j][i][0];
          r2 = rhs[k][j][i][1];
          r3 = rhs[k][j][i][2];
          r4 = rhs[k][j][i][3];
          r5 = rhs[k][j][i][4];

          t1 = bt * r1;
          t2 = 0.5 * ( r4 + r5 );

          rhs[k][j][i][0] =  bt * ( r4 - r5 );
          rhs[k][j][i][1] = -r3;
          rhs[k][j][i][2] =  r2;
          rhs[k][j][i][3] = -t1 + t2;
          rhs[k][j][i][4] =  t1 + t2;
        }
      }
    }
  }

  //---------------------------------------------------------------------
  //z_sovle
  //---------------------------------------------------------------------
  #pragma dvm region targets(HOST) local(lhs, lhsp, lhsm) inout(rhs)
  {
    #pragma dvm parallel([j][i] on rhs[][j][i][]) private(i,j,k, k1, k2, m, ru1, fac1, fac2)
    for (j = 1; j <= ny2; j++) {
      for (i = 1; i <= nx2; i++) {
        //---------------------------------------------------------------------
        //location part
        //---------------------------------------------------------------------
        double cv_loc[KMAXP+1];
        double rhos_loc[KMAXP+1];

        int nj = nz2+1;
        for (m = 0; m < 5; m++) {
          lhs [0][j][i][m] = 0.0;
          lhsp[0][j][i][m] = 0.0;
          lhsm[0][j][i][m] = 0.0;
          lhs [nj][j][i][m] = 0.0;
          lhsp[nj][j][i][m] = 0.0;
          lhsm[nj][j][i][m] = 0.0;
        }
        lhs [0][j][i][2] = 1.0;
        lhsp[0][j][i][2] = 1.0;
        lhsm[0][j][i][2] = 1.0;
        lhs [nj][j][i][2] = 1.0;
        lhsp[nj][j][i][2] = 1.0;
        lhsm[nj][j][i][2] = 1.0;
        //---------------------------------------------------------------------
        // Computes the left hand side for the three z-factors   
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // first fill the lhs for the u-eigenvalue                          
        //---------------------------------------------------------------------
        for (k = 0; k <= nz2+1; k++) {
          ru1 = c3c4*rho_i[k][j][i];
          cv_loc[k] = ws[k][j][i];
          rhos_loc[k] = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));
        }

        for (k = 1; k <= nz2; k++) {
          lhs[k][j][i][0] =  0.0;
          lhs[k][j][i][1] = -dttz2 * cv_loc[k-1] - dttz1 * rhos_loc[k-1];
          lhs[k][j][i][2] =  1.0 + c2dttz1 * rhos_loc[k];
          lhs[k][j][i][3] =  dttz2 * cv_loc[k+1] - dttz1 * rhos_loc[k+1];
          lhs[k][j][i][4] =  0.0;
        }
        //---------------------------------------------------------------------
        // add fourth order dissipation                                  
        //---------------------------------------------------------------------

        k = 1;
        lhs[k][j][i][2] += comz5;
        lhs[k][j][i][3] -= comz4;
        lhs[k][j][i][4] += comz1;
        lhs[k+1][j][i][1] -= comz4;
        lhs[k+1][j][i][2] += comz6;
        lhs[k+1][j][i][3] -= comz4;
        lhs[k+1][j][i][4] += comz1;

        for (k = 3; k <= grid_points[2]-4; k++) {
          lhs[k][j][i][0] += comz1;
          lhs[k][j][i][1] -= comz4;
          lhs[k][j][i][2] += comz6;
          lhs[k][j][i][3] -= comz4;
          lhs[k][j][i][4] += comz1;
        }

        k = nz2-1;
        lhs[k][j][i][0] += comz1;
        lhs[k][j][i][1] -= comz4;
        lhs[k][j][i][2] += comz6;
        lhs[k][j][i][3] -= comz4;
        lhs[k+1][j][i][0] += comz1;
        lhs[k+1][j][i][1] -= comz4;
        lhs[k+1][j][i][2] += comz5;
        //---------------------------------------------------------------------
        // subsequently, fill the other factors (u+c), (u-c) 
        //---------------------------------------------------------------------
        for (k = 1; k <= grid_points[2]-2; k++) {
          lhsp[k][j][i][0] = lhs[k][j][i][0];
          lhsp[k][j][i][1] = lhs[k][j][i][1] - dttz2 * speed[k-1][j][i];
          lhsp[k][j][i][2] = lhs[k][j][i][2];
          lhsp[k][j][i][3] = lhs[k][j][i][3] + dttz2 * speed[k+1][j][i];
          lhsp[k][j][i][4] = lhs[k][j][i][4];
          lhsm[k][j][i][0] = lhs[k][j][i][0];
          lhsm[k][j][i][1] = lhs[k][j][i][1] + dttz2 * speed[k-1][j][i];
          lhsm[k][j][i][2] = lhs[k][j][i][2];
          lhsm[k][j][i][3] = lhs[k][j][i][3] - dttz2 * speed[k+1][j][i];
          lhsm[k][j][i][4] = lhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // FORWARD ELIMINATION  
        //---------------------------------------------------------------------
        for (k = 0; k <= grid_points[2]-3; k++) {
          k1 = k + 1;
          k2 = k + 2;

          // LHS
          fac1 = 1.0/lhs[k][j][i][2];
          lhs[k][j][i][3] *= fac1;
          lhs[k][j][i][4] *= fac1;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] *= fac1;
          }
          lhs[k1][j][i][2] -= lhs[k1][j][i][1]*lhs[k][j][i][3];
          lhs[k1][j][i][3] -= lhs[k1][j][i][1]*lhs[k][j][i][4];
          for (m = 0; m < 3; m++) {
            rhs[k1][j][i][m] -= lhs[k1][j][i][1]*rhs[k][j][i][m];
          }
          lhs[k2][j][i][1] -= lhs[k2][j][i][0]*lhs[k][j][i][3];
          lhs[k2][j][i][2] -= lhs[k2][j][i][0]*lhs[k][j][i][4];
          for (m = 0; m < 3; m++) {
            rhs[k2][j][i][m] -= lhs[k2][j][i][0]*rhs[k][j][i][m];
          }
          //---------------------------------------------------------------------
          // for (the u+c and the u-c factors               
          //---------------------------------------------------------------------
          //m = 3;
          fac1 = 1.0/lhsp[k][j][i][2];
          lhsp[k][j][i][3]   *= fac1;
          lhsp[k][j][i][4]   *= fac1;
          rhs[k][j][i][3]  *= fac1;
          lhsp[k1][j][i][2]  -= lhsp[k1][j][i][1]*lhsp[k][j][i][3];
          lhsp[k1][j][i][3]  -= lhsp[k1][j][i][1]*lhsp[k][j][i][4];
          rhs[k1][j][i][3] -= lhsp[k1][j][i][1]*rhs[k][j][i][3];
          lhsp[k2][j][i][1]  -= lhsp[k2][j][i][0]*lhsp[k][j][i][3];
          lhsp[k2][j][i][2]  -= lhsp[k2][j][i][0]*lhsp[k][j][i][4];
          rhs[k2][j][i][3] -= lhsp[k2][j][i][0]*rhs[k][j][i][3];

          //m = 4;
          fac1 = 1.0/lhsm[k][j][i][2];
          lhsm[k][j][i][3]   *= fac1;
          lhsm[k][j][i][4]   *= fac1;
          rhs[k][j][i][4]  *= fac1;
          lhsm[k1][j][i][2]  -= lhsm[k1][j][i][1]*lhsm[k][j][i][3];
          lhsm[k1][j][i][3]  -= lhsm[k1][j][i][1]*lhsm[k][j][i][4];
          rhs[k1][j][i][4] -= lhsm[k1][j][i][1]*rhs[k][j][i][4];
          lhsm[k2][j][i][1]  -= lhsm[k2][j][i][0]*lhsm[k][j][i][3];
          lhsm[k2][j][i][2]  -= lhsm[k2][j][i][0]*lhsm[k][j][i][4];
          rhs[k2][j][i][4] -= lhsm[k2][j][i][0]*rhs[k][j][i][4];
        }
        //---------------------------------------------------------------------
        // The last two rows in this grid block are a bit different, 
        // since they for (not have two more rows available for the
        // elimination of off-diagonal entries
        //---------------------------------------------------------------------
        k  = grid_points[2]-2;
        k1 = grid_points[2]-1;
        fac1 = 1.0/lhs[k][j][i][2];
        lhs[k][j][i][3] *= fac1;
        lhs[k][j][i][4] *= fac1;
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[k1][j][i][2] -= lhs[k1][j][i][1]*lhs[k][j][i][3];
        lhs[k1][j][i][3] -= lhs[k1][j][i][1]*lhs[k][j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k1][j][i][m] -= lhs[k1][j][i][1]*rhs[k][j][i][m];
        }

        //---------------------------------------------------------------------
        // scale the last row immediately
        //---------------------------------------------------------------------
        fac2 = 1.0/lhs[k1][j][i][2];
        for (m = 0; m < 3; m++) {
          rhs[k1][j][i][m] = fac2*rhs[k1][j][i][m];
        }
        //---------------------------------------------------------------------
        // And again the last two rows separately
        //---------------------------------------------------------------------
        //m = 3;
        fac1 = 1.0/lhsp[k][j][i][2];
        lhsp[k][j][i][3]   *= fac1;
        lhsp[k][j][i][4]   *= fac1;
        rhs[k][j][i][3]  *= fac1;
        lhsp[k1][j][i][2]  -= lhsp[k1][j][i][1]*lhsp[k][j][i][3];
        lhsp[k1][j][i][3]  -= lhsp[k1][j][i][1]*lhsp[k][j][i][4];
        rhs[k1][j][i][3] -= lhsp[k1][j][i][1]*rhs[k][j][i][3];
        rhs[k1][j][i][3] /= lhsp[k1][j][i][2];

        //m = 4;
        fac1 = 1.0/lhsm[k][j][i][2];
        lhsm[k][j][i][3]   *= fac1;
        lhsm[k][j][i][4]   *= fac1;
        rhs[k][j][i][4]  *= fac1;
        lhsm[k1][j][i][2]  -= lhsm[k1][j][i][1]*lhsm[k][j][i][3];
        lhsm[k1][j][i][3]  -= lhsm[k1][j][i][1]*lhsm[k][j][i][4];
        rhs[k1][j][i][4] -= lhsm[k1][j][i][1]*rhs[k][j][i][4];
        rhs[k1][j][i][4] /= lhsm[k1][j][i][2];
        //---------------------------------------------------------------------
        // BACKSUBSTITUTION 
        //---------------------------------------------------------------------
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] -= lhs[k][j][i][3]*rhs[k1][j][i][m];
        }

        rhs[k][j][i][3] -= lhsp[k][j][i][3]*rhs[k1][j][i][3];
        rhs[k][j][i][4] -= lhsm[k][j][i][3]*rhs[k1][j][i][4];
        //---------------------------------------------------------------------
        // Whether or not this is the last processor, we always have
        // to complete the back-substitution 
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        // The first three factors
        //---------------------------------------------------------------------
        for (k = grid_points[2]-3; k >= 0; k--) {
          k1 = k + 1;
          k2 = k + 2;
          for (m = 0; m < 3; m++) {
            rhs[k][j][i][m] -= lhs[k][j][i][3]*rhs[k1][j][i][m] +
                               lhs[k][j][i][4]*rhs[k2][j][i][m];
          }

          //-------------------------------------------------------------------
          // And the remaining two
          //-------------------------------------------------------------------
          rhs[k][j][i][3] -= lhsp[k][j][i][3]*rhs[k1][j][i][3] +
                             lhsp[k][j][i][4]*rhs[k2][j][i][3];
          rhs[k][j][i][4] -= lhsm[k][j][i][3]*rhs[k1][j][i][4] +
                             lhsm[k][j][i][4]*rhs[k2][j][i][4];
        }
      }
    }
  }

  //-------------------------------------------------------------------
  //tzeter
  //-------------------------------------------------------------------
  #pragma dvm region targets(HOST) in(u, us, vs, ws, qs, speed, rhs) out(rhs)
  {
    #pragma dvm parallel([k][j][i] on rhs[k][j][i][]) private(i,j,k,xvel,yvel,zvel,ac,r1,r2,r3,r4,r5,t1,t2,t3,btuz, ac2u, uzik1)
    for (k = 1; k <= nz2; k++) {
      for (j = 1; j <= ny2; j++) {
        for (i = 1; i <= nx2; i++) {
          xvel = us[k][j][i];
          yvel = vs[k][j][i];
          zvel = ws[k][j][i];
          ac   = speed[k][j][i];

          ac2u = ac*ac;

          r1 = rhs[k][j][i][0];
          r2 = rhs[k][j][i][1];
          r3 = rhs[k][j][i][2];
          r4 = rhs[k][j][i][3];
          r5 = rhs[k][j][i][4];     

          uzik1 = u[k][j][i][0];
          btuz  = bt * uzik1;

          t1 = btuz/ac * (r4 + r5);
          t2 = r3 + t1;
          t3 = btuz * (r4 - r5);

          rhs[k][j][i][0] = t2;
          rhs[k][j][i][1] = -uzik1*r2 + xvel*t2;
          rhs[k][j][i][2] =  uzik1*r1 + yvel*t2;
          rhs[k][j][i][3] =  zvel*t2  + t3;
          rhs[k][j][i][4] =  uzik1*(-xvel*r2 + yvel*r1) + 
                            qs[k][j][i]*t2 + c2iv*ac2u*t1 + zvel*t3;
        }
      }
    }
  }

  #pragma dvm endinterval
  #pragma dvm get_actual(rhs)
  if (timeron) timer_stop(t_xsolve);

  if (timeron) timer_start(t_ninvr);
  if (timeron) timer_stop(t_ninvr);
  if (timeron) timer_start(t_pinvr);
  if (timeron) timer_stop(t_pinvr);
  if (timeron) timer_start(t_tzetar);
  if (timeron) timer_stop(t_tzetar);
  if (timeron) timer_start(t_ysolve);
  if (timeron) timer_stop(t_ysolve);
  if (timeron) timer_start(t_zsolve);
  if (timeron) timer_stop(t_zsolve);
}

