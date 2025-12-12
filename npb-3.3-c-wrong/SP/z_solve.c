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
// step in the z-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the z-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void z_solve()
{
  // int i, j, k, k1, k2, m;
  // double ru1, fac1, fac2;
  // double lhs__[5][3], lhsp__[5][3], lhsm__[5][3];
  // double rhs__[5][3];
  // double t1, t2, t3, ac, xvel, yvel, zvel,btuz, ac2u, uzik1;

  // if (timeron) timer_start(t_zsolve);
  // #pragma dvm region local(lhs)
  // {
  //   #pragma dvm parallel([j][i] on u[][j][i][]) \
  //   private(m,k,ru1,k1,k2,fac1,fac2,lhs__,lhsp__,lhsm__,rhs__,t1,t2,t3,ac,xvel,yvel,zvel,btuz,ac2u,uzik1)
  //     for (j = 1; j <= ny2; j++) {
  //       for (i = 1; i <= nx2; i++) {

  //         lhs__[0][0]  = 0.0;
  //         lhsp__[0][0] = 0.0;
  //         lhsm__[0][0] = 0.0;
  //         lhs__[1][0]  = 0.0;
  //         lhsp__[1][0] = 0.0;
  //         lhsm__[1][0] = 0.0;
  //         lhs__[2][0]  = 1.0;
  //         lhsp__[2][0] = 1.0;
  //         lhsm__[2][0] = 1.0;
  //         lhs__[3][0]  = 0.0;
  //         lhsp__[3][0] = 0.0;
  //         lhsm__[3][0] = 0.0;
  //         lhs__[4][0]  = 0.0;
  //         lhsp__[4][0] = 0.0;
  //         lhsm__[4][0] = 0.0;

  //         lhs__[0][1] = 0.0;
  //         ru1 = c3c4 / u[1-1][j][i][0];
  //         ru1 = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));

  //         lhs__[1][1] = -dttz2 * ws[1-1][j][i] - dttz1 * ru1;
  //         ru1 = c3c4 / u[1][j][i][0];
  //         ru1 = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));
          
  //         lhs__[2][1] = 1.0 + c2dttz1 * ru1;
  //         ru1 = c3c4 / u[1+1][j][i][0];
  //         ru1 = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));

  //         lhs__[3][1] = dttz2 * ws[1+1][j][i] - dttz1 * ru1;
  //         lhs__[4][1] = 0.0;

  //         lhs__[2][1] += comz5;
  //         lhs__[3][1] -= comz4;
  //         lhs__[4][1] += comz1;

  //         lhsp__[0][1] = lhs__[0][1];
  //         lhsp__[1][1] = lhs__[1][1] - dttz2 * speed[0][j][i];
  //         lhsp__[2][1] = lhs__[2][1];
  //         lhsp__[3][1] = lhs__[3][1] + dttz2 * speed[2][j][i];
  //         lhsp__[4][1] = lhs__[4][1];

  //         lhsm__[0][1] = lhs__[0][1];
  //         lhsm__[1][1] = lhs__[1][1] + dttz2 * speed[0][j][i];
  //         lhsm__[2][1] = lhs__[2][1];
  //         lhsm__[3][1] = lhs__[3][1] - dttz2 * speed[2][j][i];
  //         lhsm__[4][1] = lhs__[4][1];

  //         for (k = 0; k <= nz2+1; k++) {
  //           if (k + 2 < nz2+1) {
  //               m = k + 2;

  //               lhs__[0][2] = 0.0;
  //               ru1 = c3c4 / u[m-1][j][i][0];
  //               ru1 = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));

  //               lhs__[1][2] = -dttz2 * ws[m-1][j][i] - dttz1 * ru1;
  //               ru1 = c3c4 / u[m][j][i][0];
  //               ru1 = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));
                
  //               lhs__[2][2] = 1.0 + c2dttz1 * ru1;
  //               ru1 = c3c4 / u[m+1][j][i][0];
  //               ru1 = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));

  //               lhs__[3][2] = dttz2 * ws[m+1][j][i] - dttz1 * ru1;
  //               lhs__[4][2] = 0.0;

  //               /* 根据 m_i 位置加 comz 修正 */
  //               if (m == 1) {
  //                   lhs__[2][2] += comz5;
  //                   lhs__[3][2] -= comz4;
  //                   lhs__[4][2] += comz1;
  //               } 
  //               else if (m == 2) {
  //                   lhs__[1][2] -= comz4;
  //                   lhs__[2][2] += comz6;
  //                   lhs__[3][2] -= comz4;
  //                   lhs__[4][2] += comz1;
  //               } 
  //               else if (m >= 3 && m <= nz2 - 2) {
  //                   lhs__[0][2] += comz1;
  //                   lhs__[1][2] -= comz4;
  //                   lhs__[2][2] += comz6;
  //                   lhs__[3][2] -= comz4;
  //                   lhs__[4][2] += comz1;
  //               } 
  //               else if (m == nz2 - 1) {
  //                   lhs__[0][2] += comz1;
  //                   lhs__[1][2] -= comz4;
  //                   lhs__[2][2] += comz6;
  //                   lhs__[3][2] -= comz4;
  //               } 
  //               else if (m == nz2) {
  //                   lhs__[0][2] += comz1;
  //                   lhs__[1][2] -= comz4;
  //                   lhs__[2][2] += comz5;
  //               }

  //               /* lhsp__/lhsm__ 列 2 */
  //               lhsp__[0][2] = lhs__[0][2];
  //               lhsp__[1][2] = lhs__[1][2] - dttz2 * speed[m-1][j][i];
  //               lhsp__[2][2] = lhs__[2][2];
  //               lhsp__[3][2] = lhs__[3][2] + dttz2 * speed[m+1][j][i];
  //               lhsp__[4][2] = lhs__[4][2];

  //               lhsm__[0][2] = lhs__[0][2];
  //               lhsm__[1][2] = lhs__[1][2] + dttz2 * speed[m-1][j][i];
  //               lhsm__[2][2] = lhs__[2][2];
  //               lhsm__[3][2] = lhs__[3][2] - dttz2 * speed[m+1][j][i];
  //               lhsm__[4][2] = lhs__[4][2];

  //           } 
  //           else if (k + 2 == nz2 + 1) {
  //             lhs__[0][2]  = 0.0;
  //             lhsp__[0][2] = 0.0;
  //             lhsm__[0][2] = 0.0;
  //             lhs__[1][2]  = 0.0;
  //             lhsp__[1][2] = 0.0;
  //             lhsm__[1][2] = 0.0;
  //             lhs__[2][2]  = 1.0;
  //             lhsp__[2][2] = 1.0;
  //             lhsm__[2][2] = 1.0;
  //             lhs__[3][2]  = 0.0;
  //             lhsp__[3][2] = 0.0;
  //             lhsm__[3][2] = 0.0;
  //             lhs__[4][2]  = 0.0;
  //             lhsp__[4][2] = 0.0;
  //             lhsm__[4][2] = 0.0;
  //           }

  //           k1 = k + 1;
  //           k2 = k + 2;

  //           fac1 = 1.0 / lhs__[2][0];     /* lhs__(3,0) */
  //           lhs__[3][0] *= fac1;          /* lhs__(4,0) */
  //           lhs__[4][0] *= fac1;          /* lhs__(5,0) */

  //           for (m = 0; m < 3; m++) {
  //               rhs[k][j][i][m] *= fac1;  /* rhs(m+1,i,j,k) */
  //           }

  //           if (k <= nz2 - 1) {
  //               lhs__[2][1] -= lhs__[1][1] * lhs__[3][0];
  //               lhs__[3][1] -= lhs__[1][1] * lhs__[4][0];

  //               lhs__[1][2] -= lhs__[0][2] * lhs__[3][0];
  //               lhs__[2][2] -= lhs__[0][2] * lhs__[4][0];

  //               for (m = 0; m < 3; m++) {
  //                   rhs[k1][j][i][m] -= lhs__[1][1] * rhs[k][j][i][m];
  //                   rhs[k2][j][i][m] -= lhs__[0][2] * rhs[k][j][i][m];
  //               }
  //           } 
  //           else {
  //               lhs__[2][1] -= lhs__[1][1] * lhs__[3][0];
  //               lhs__[3][1] -= lhs__[1][1] * lhs__[4][0];

  //               fac2 = 1.0 / lhs__[2][1];
  //               for (m = 0; m < 3; m++) {
  //                   rhs[k1][j][i][m] -= lhs__[1][1] * rhs[k][j][i][m];
  //                   rhs[k1][j][i][m] *= fac2;
  //               }
  //           }

  //           m = 3;
  //           fac1 = 1.0 / lhsp__[2][0];    /* lhsp__(3,0) */
  //           lhsp__[3][0] *= fac1;         /* lhsp__(4,0) */
  //           lhsp__[4][0] *= fac1;         /* lhsp__(5,0) */
  //           rhs[k][j][i][m] *= fac1;

  //           lhsp__[2][1] -= lhsp__[1][1] * lhsp__[3][0];
  //           lhsp__[3][1] -= lhsp__[1][1] * lhsp__[4][0];
  //           rhs[k1][j][i][m] -= lhsp__[1][1] * rhs[k][j][i][m];

  //           if (k < nz2) {
  //               lhsp__[1][2] -= lhsp__[0][2] * lhsp__[3][0];
  //               lhsp__[2][2] -= lhsp__[0][2] * lhsp__[4][0];
  //               rhs[k2][j][i][m] -= lhsp__[0][2] * rhs[k][j][i][m];
  //           }

  //           m = 4; 
  //           fac1 = 1.0 / lhsm__[2][0];
  //           lhsm__[3][0] *= fac1;
  //           lhsm__[4][0] *= fac1;
  //           rhs[k][j][i][m] *= fac1;

  //           lhsm__[2][1] -= lhsm__[1][1] * lhsm__[3][0];
  //           lhsm__[3][1] -= lhsm__[1][1] * lhsm__[4][0];
  //           rhs[k1][j][i][m] -= lhsm__[1][1] * rhs[k][j][i][m];

  //           if (k < nz2) {
  //               lhsm__[1][2] -= lhsm__[0][2] * lhsm__[3][0];
  //               lhsm__[2][2] -= lhsm__[0][2] * lhsm__[4][0];
  //               rhs[k2][j][i][m] -= lhsm__[0][2] * rhs[k][j][i][m];
  //           }

  //           if (k == nz2) {
  //               rhs[k1][j][i][3] /= lhsp__[2][1];   /* rhs(4,i1,j,k) / lhsp__(3,1) */
  //               rhs[k1][j][i][4] /= lhsm__[2][1];   /* rhs(5,i1,j,k) / lhsm__(3,1) */

  //               for (m = 0; m < 3; m++) {
  //                   rhs[k][j][i][m] -= lhs__[3][0] * rhs[k1][j][i][m];
  //               }
  //               rhs[k][j][i][3] -= lhsp__[3][0] * rhs[k1][j][i][3];
  //               rhs[k][j][i][4] -= lhsm__[3][0] * rhs[k1][j][i][4];
  //           }


  //           lhs[k][j][i][3][0] = lhs__[3][0];   /* lhs(0,4,i,j,k) */
  //           lhs[k][j][i][3][1] = lhsp__[3][0];  /* lhs(1,4,i,j,k) */
  //           lhs[k][j][i][3][2] = lhsm__[3][0];  /* lhs(2,4,i,j,k) */

  //           lhs[k][j][i][4][0] = lhs__[4][0];   /* lhs(0,5,i,j,k) */
  //           lhs[k][j][i][4][1] = lhsp__[4][0];  /* lhs(1,5,i,j,k) */
  //           lhs[k][j][i][4][2] = lhsm__[4][0];  /* lhs(2,5,i,j,k) */

  //           for (m = 0; m < 5; m++) {
  //               lhs__[m][0]  = lhs__[m][1];
  //               lhsp__[m][0] = lhsp__[m][1];
  //               lhsm__[m][0] = lhsm__[m][1];

  //               lhs__[m][1]  = lhs__[m][2];
  //               lhsp__[m][1] = lhsp__[m][2];
  //               lhsm__[m][1] = lhsm__[m][2];
  //           }
  //         } /* end for k */

  //       k = PROBLEM_SIZE - 3;

  //       for (m = 0; m < 5; m++) {
  //           rhs__[m][2] = rhs[k+2][j][i][m];
  //           rhs__[m][1] = rhs[k+1][j][i][m];
  //           rhs__[m][0] = rhs[k][j][i][m];
  //       }

  //       for (m = 0; m < 3; m++) {
  //           rhs__[m][0] -= lhs[k][j][i][3][0] * rhs__[m][1] +
  //                           lhs[k][j][i][4][0] * rhs__[m][2];
  //       }
  //       rhs__[3][0] -= lhs[k][j][i][3][1] * rhs__[3][1]
  //                     + lhs[k][j][i][4][1] * rhs__[3][2];
  //       rhs__[4][0] -= lhs[k][j][i][3][2] * rhs__[4][1]
  //                     + lhs[k][j][i][4][2] * rhs__[4][2];

  //       for (m = 0; m < 5; m++) {
  //           rhs__[m][2] = rhs__[m][1];
  //           rhs__[m][1] = rhs__[m][0];
  //       }

  //       for (k = PROBLEM_SIZE - 4; k >= 0; k--) {

  //           for (m = 0; m < 5; m++) {
  //               rhs__[m][0] = rhs[k][j][i][m];
  //           }

  //           for (m = 0; m < 3; m++) {
  //               rhs__[m][0] -= lhs[k][j][i][3][0] * rhs__[m][1]
  //                             + lhs[k][j][i][4][0] * rhs__[m][2];
  //           }
  //           rhs__[3][0] -= lhs[k][j][i][3][1] * rhs__[3][1]
  //                         + lhs[k][j][i][4][1] * rhs__[3][2];
  //           rhs__[4][0] -= lhs[k][j][i][3][2] * rhs__[4][1]
  //                         + lhs[k][j][i][4][2] * rhs__[4][2];

  //           xvel = us[k+2][j][i];
  //           yvel = vs[k+2][j][i];
  //           zvel = ws[k+2][j][i];
  //           ac = speed[k+2][j][i];
  //           ac2u = ac*ac;
  //           uzik1 = u[k+2][j][i][0];
  //           btuz = bt * uzik1;
  //           t1 = btuz / ac * (rhs__[3][2] + rhs__[4][2]);
  //           t2 = rhs__[2][2] + t1;
  //           t3 = btuz * (rhs__[3][2] - rhs__[4][2]);

  //           rhs__[2][2] = uzik1 * rhs__[0][2] + yvel * t2;
  //           rhs__[3][2] = zvel * t2 + t3;
  //           rhs__[4][2] = uzik1*(-xvel*rhs__[1][2]+yvel*rhs__[0][2]) +
  //                         qs[k+2][j][i] * t2 + c2iv*ac2u*t1 + zvel*t3;
  //           rhs__[0][2] = t2;
  //           rhs__[1][2] = -uzik1 * rhs__[1][2] + xvel * t2;

  //           u[k+2][j][i][0] += rhs__[0][2];
  //           u[k+2][j][i][1] += rhs__[1][2];
  //           u[k+2][j][i][2] += rhs__[2][2];
  //           u[k+2][j][i][3] += rhs__[3][2];
  //           u[k+2][j][i][4] += rhs__[4][2];

  //           for (m = 0; m < 5; m++) {
  //               rhs__[m][2] = rhs__[m][1];
  //               rhs__[m][1] = rhs__[m][0];
  //           }
  //       }

  //       xvel = us[k+2][j][i];
  //       yvel = vs[k+2][j][i];
  //       zvel = ws[k+2][j][i];
  //       ac = speed[k+2][j][i];
  //       ac2u = ac*ac;
  //       uzik1 = u[k+2][j][i][0];
  //       btuz = bt * uzik1;
  //       t1 = btuz / ac * (rhs__[3][2] + rhs__[4][2]);
  //       t2 = rhs__[2][2] + t1;
  //       t3 = btuz * (rhs__[3][2] - rhs__[4][2]);

  //       rhs__[2][2] = uzik1 * rhs__[0][2] + yvel * t2;
  //       rhs__[3][2] = zvel * t2 + t3;
  //       rhs__[4][2] = uzik1*(-xvel*rhs__[1][2]+yvel*rhs__[0][2]) +
  //                     qs[k+2][j][i] * t2 + c2iv*ac2u*t1 + zvel*t3;
  //       rhs__[0][2] = t2;
  //       rhs__[1][2] = -uzik1 * rhs__[1][2] + xvel * t2;

  //       u[k+2][j][i][0] += rhs__[0][2];
  //       u[k+2][j][i][1] += rhs__[1][2];
  //       u[k+2][j][i][2] += rhs__[2][2];
  //       u[k+2][j][i][3] += rhs__[3][2];
  //       u[k+2][j][i][4] += rhs__[4][2];
  //     }
  //   }
  // }
  // if (timeron) timer_stop(t_zsolve);
}
