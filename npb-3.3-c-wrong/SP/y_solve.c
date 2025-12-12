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
// step in the y-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the y-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void y_solve()
{
//   int i, j, k, j1, j2, m;
//   double ru1, fac1, fac2,t1,t2;

//   double lhs__[5][3], lhsp__[5][3], lhsm__[5][3];
//   double rhs__[5][3];

//   if (timeron) timer_start(t_ysolve);

//   #pragma dvm region local(lhs)
//   {
//     #pragma dvm parallel([k][i] on u[k][][i][]) \
//     private(m,j,ru1,j1,j2,fac1,fac2,lhs__,lhsp__,lhsm__,rhs__,t1,t2)
//       for (k = 1; k <= nz2; k++) {
//         for (i = 1; i <= nx2; i++) {

//           lhs__[0][0]  = 0.0;
//           lhsp__[0][0] = 0.0;
//           lhsm__[0][0] = 0.0;
//           lhs__[1][0]  = 0.0;
//           lhsp__[1][0] = 0.0;
//           lhsm__[1][0] = 0.0;
//           lhs__[2][0]  = 1.0;
//           lhsp__[2][0] = 1.0;
//           lhsm__[2][0] = 1.0;
//           lhs__[3][0]  = 0.0;
//           lhsp__[3][0] = 0.0;
//           lhsm__[3][0] = 0.0;
//           lhs__[4][0]  = 0.0;
//           lhsp__[4][0] = 0.0;
//           lhsm__[4][0] = 0.0;

//           lhs__[0][1] = 0.0;
//           ru1 = c3c4 / u[k][1-1][i][0];
//           ru1 = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));

//           lhs__[1][1] = -dtty2 * vs[k][1-1][i] - dtty1 * ru1;
//           ru1 = c3c4 / u[k][1][i][0];
//           ru1 = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
          
//           lhs__[2][1] = 1.0 + c2dtty1 * ru1;
//           ru1 = c3c4 / u[k][1+1][i][0];
//           ru1 = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));

//           lhs__[3][1] = dtty2 * vs[k][1+1][i] - dtty1 * ru1;
//           lhs__[4][1] = 0.0;

//           lhs__[2][1] += comz5;
//           lhs__[3][1] -= comz4;
//           lhs__[4][1] += comz1;

//           lhsp__[0][1] = lhs__[0][1];
//           lhsp__[1][1] = lhs__[1][1] - dtty2 * speed[k][0][i];
//           lhsp__[2][1] = lhs__[2][1];
//           lhsp__[3][1] = lhs__[3][1] + dtty2 * speed[k][2][i];
//           lhsp__[4][1] = lhs__[4][1];

//           lhsm__[0][1] = lhs__[0][1];
//           lhsm__[1][1] = lhs__[1][1] + dtty2 * speed[k][0][i];
//           lhsm__[2][1] = lhs__[2][1];
//           lhsm__[3][1] = lhs__[3][1] - dtty2 * speed[k][2][i];
//           lhsm__[4][1] = lhs__[4][1];

//           for (j = 0; j <= ny2+1; j++) {
//             if (j + 2 < ny2+1) {
//                 m = j + 2;

//                 lhs__[0][2] = 0.0;
//                 ru1 = c3c4 / u[k][m-1][i][0];
//                 ru1 = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));

//                 lhs__[1][2] = -dtty2 * vs[k][m-1][i] - dtty1 * ru1;
//                 ru1 = c3c4 / u[k][m][i][0];
//                 ru1 = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
                
//                 lhs__[2][2] = 1.0 + c2dtty1 * ru1;
//                 ru1 = c3c4 / u[k][m+1][i][0];
//                 ru1 = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));

//                 lhs__[3][2] = dtty2 * vs[k][m+1][i] - dtty1 * ru1;
//                 lhs__[4][2] = 0.0;

//                 /* 根据 m_i 位置加 comz 修正 */
//                 if (m == 1) {
//                     lhs__[2][2] += comz5;
//                     lhs__[3][2] -= comz4;
//                     lhs__[4][2] += comz1;
//                 } 
//                 else if (m == 2) {
//                     lhs__[1][2] -= comz4;
//                     lhs__[2][2] += comz6;
//                     lhs__[3][2] -= comz4;
//                     lhs__[4][2] += comz1;
//                 } 
//                 else if (m >= 3 && m <= ny2 - 2) {
//                     lhs__[0][2] += comz1;
//                     lhs__[1][2] -= comz4;
//                     lhs__[2][2] += comz6;
//                     lhs__[3][2] -= comz4;
//                     lhs__[4][2] += comz1;
//                 } 
//                 else if (m == ny2 - 1) {
//                     lhs__[0][2] += comz1;
//                     lhs__[1][2] -= comz4;
//                     lhs__[2][2] += comz6;
//                     lhs__[3][2] -= comz4;
//                 } 
//                 else if (m == ny2) {
//                     lhs__[0][2] += comz1;
//                     lhs__[1][2] -= comz4;
//                     lhs__[2][2] += comz5;
//                 }

//                 /* lhsp__/lhsm__ 列 2 */
//                 lhsp__[0][2] = lhs__[0][2];
//                 lhsp__[1][2] = lhs__[1][2] - dtty2 * speed[k][m-1][i];
//                 lhsp__[2][2] = lhs__[2][2];
//                 lhsp__[3][2] = lhs__[3][2] + dtty2 * speed[k][m+1][i];
//                 lhsp__[4][2] = lhs__[4][2];

//                 lhsm__[0][2] = lhs__[0][2];
//                 lhsm__[1][2] = lhs__[1][2] + dtty2 * speed[k][m-1][i];
//                 lhsm__[2][2] = lhs__[2][2];
//                 lhsm__[3][2] = lhs__[3][2] - dtty2 * speed[k][m+1][i];
//                 lhsm__[4][2] = lhs__[4][2];

//             } 
//             else if (j + 2 == ny2 + 1) {
//               lhs__[0][2]  = 0.0;
//               lhsp__[0][2] = 0.0;
//               lhsm__[0][2] = 0.0;
//               lhs__[1][2]  = 0.0;
//               lhsp__[1][2] = 0.0;
//               lhsm__[1][2] = 0.0;
//               lhs__[2][2]  = 1.0;
//               lhsp__[2][2] = 1.0;
//               lhsm__[2][2] = 1.0;
//               lhs__[3][2]  = 0.0;
//               lhsp__[3][2] = 0.0;
//               lhsm__[3][2] = 0.0;
//               lhs__[4][2]  = 0.0;
//               lhsp__[4][2] = 0.0;
//               lhsm__[4][2] = 0.0;
//             }

//             j1 = j + 1;
//             j2 = j + 2;

//             fac1 = 1.0 / lhs__[2][0];     /* lhs__(3,0) */
//             lhs__[3][0] *= fac1;          /* lhs__(4,0) */
//             lhs__[4][0] *= fac1;          /* lhs__(5,0) */

//             for (m = 0; m < 3; m++) {
//                 rhs[k][j][i][m] *= fac1;  /* rhs(m+1,i,j,k) */
//             }

//             if (j <= ny2 - 1) {
//                 lhs__[2][1] -= lhs__[1][1] * lhs__[3][0];
//                 lhs__[3][1] -= lhs__[1][1] * lhs__[4][0];

//                 lhs__[1][2] -= lhs__[0][2] * lhs__[3][0];
//                 lhs__[2][2] -= lhs__[0][2] * lhs__[4][0];

//                 for (m = 0; m < 3; m++) {
//                     rhs[k][j1][i][m] -= lhs__[1][1] * rhs[k][j][i][m];
//                     rhs[k][j2][i][m] -= lhs__[0][2] * rhs[k][j][i][m];
//                 }
//             } 
//             else {
//                 lhs__[2][1] -= lhs__[1][1] * lhs__[3][0];
//                 lhs__[3][1] -= lhs__[1][1] * lhs__[4][0];

//                 fac2 = 1.0 / lhs__[2][1];
//                 for (m = 0; m < 3; m++) {
//                     rhs[k][j1][i][m] -= lhs__[1][1] * rhs[k][j][i][m];
//                     rhs[k][j1][i][m] *= fac2;
//                 }
//             }

//             m = 3;
//             fac1 = 1.0 / lhsp__[2][0];    /* lhsp__(3,0) */
//             lhsp__[3][0] *= fac1;         /* lhsp__(4,0) */
//             lhsp__[4][0] *= fac1;         /* lhsp__(5,0) */
//             rhs[k][j][i][m] *= fac1;

//             lhsp__[2][1] -= lhsp__[1][1] * lhsp__[3][0];
//             lhsp__[3][1] -= lhsp__[1][1] * lhsp__[4][0];
//             rhs[k][j1][i][m] -= lhsp__[1][1] * rhs[k][j][i][m];

//             if (j < ny2) {
//                 lhsp__[1][2] -= lhsp__[0][2] * lhsp__[3][0];
//                 lhsp__[2][2] -= lhsp__[0][2] * lhsp__[4][0];
//                 rhs[k][j2][i][m] -= lhsp__[0][2] * rhs[k][j][i][m];
//             }

//             m = 4; 
//             fac1 = 1.0 / lhsm__[2][0];
//             lhsm__[3][0] *= fac1;
//             lhsm__[4][0] *= fac1;
//             rhs[k][j][i][m] *= fac1;

//             lhsm__[2][1] -= lhsm__[1][1] * lhsm__[3][0];
//             lhsm__[3][1] -= lhsm__[1][1] * lhsm__[4][0];
//             rhs[k][j1][i][m] -= lhsm__[1][1] * rhs[k][j][i][m];

//             if (j < ny2) {
//                 lhsm__[1][2] -= lhsm__[0][2] * lhsm__[3][0];
//                 lhsm__[2][2] -= lhsm__[0][2] * lhsm__[4][0];
//                 rhs[k][j2][i][m] -= lhsm__[0][2] * rhs[k][j][i][m];
//             }

//             if (j == ny2) {
//                 rhs[k][j1][i][3] /= lhsp__[2][1];   /* rhs(4,i1,j,k) / lhsp__(3,1) */
//                 rhs[k][j1][i][4] /= lhsm__[2][1];   /* rhs(5,i1,j,k) / lhsm__(3,1) */

//                 for (m = 0; m < 3; m++) {
//                     rhs[k][j][i][m] -= lhs__[3][0] * rhs[k][j1][i][m];
//                 }
//                 rhs[k][j][i][3] -= lhsp__[3][0] * rhs[k][j1][i][3];
//                 rhs[k][j][i][4] -= lhsm__[3][0] * rhs[k][j1][i][4];
//             }


//             lhs[k][j][i][3][0] = lhs__[3][0];   /* lhs(0,4,i,j,k) */
//             lhs[k][j][i][3][1] = lhsp__[3][0];  /* lhs(1,4,i,j,k) */
//             lhs[k][j][i][3][2] = lhsm__[3][0];  /* lhs(2,4,i,j,k) */

//             lhs[k][j][i][4][0] = lhs__[4][0];   /* lhs(0,5,i,j,k) */
//             lhs[k][j][i][4][1] = lhsp__[4][0];  /* lhs(1,5,i,j,k) */
//             lhs[k][j][i][4][2] = lhsm__[4][0];  /* lhs(2,5,i,j,k) */

//             for (m = 0; m < 5; m++) {
//                 lhs__[m][0]  = lhs__[m][1];
//                 lhsp__[m][0] = lhsp__[m][1];
//                 lhsm__[m][0] = lhsm__[m][1];

//                 lhs__[m][1]  = lhs__[m][2];
//                 lhsp__[m][1] = lhsp__[m][2];
//                 lhsm__[m][1] = lhsm__[m][2];
//             }
//           } /* end for i */

//         j = PROBLEM_SIZE - 3;

//         for (m = 0; m < 5; m++) {
//             rhs__[m][2] = rhs[k][j+2][i][m];
//             rhs__[m][1] = rhs[k][j+1][i][m];
//             rhs__[m][0] = rhs[k][j][i][m];
//         }

//         for (m = 0; m < 3; m++) {
//             rhs__[m][0] -= lhs[k][j][i][3][0] * rhs__[m][1] +
//                             lhs[k][j][i][4][0] * rhs__[m][2];
//         }
//         rhs__[3][0] -= lhs[k][j][i][3][1] * rhs__[3][1]
//                       + lhs[k][j][i][4][1] * rhs__[3][2];
//         rhs__[4][0] -= lhs[k][j][i][3][2] * rhs__[4][1]
//                       + lhs[k][j][i][4][2] * rhs__[4][2];

//         for (m = 0; m < 5; m++) {
//             rhs__[m][2] = rhs__[m][1];
//             rhs__[m][1] = rhs__[m][0];
//         }

//         for (j = PROBLEM_SIZE - 4; j >= 0; j--) {

//             for (m = 0; m < 5; m++) {
//                 rhs__[m][0] = rhs[k][j][i][m];
//             }

//             for (m = 0; m < 3; m++) {
//                 rhs__[m][0] -= lhs[k][j][i][3][0] * rhs__[m][1]
//                               + lhs[k][j][i][4][0] * rhs__[m][2];
//             }
//             rhs__[3][0] -= lhs[k][j][i][3][1] * rhs__[3][1]
//                           + lhs[k][j][i][4][1] * rhs__[3][2];
//             rhs__[4][0] -= lhs[k][j][i][3][2] * rhs__[4][1]
//                           + lhs[k][j][i][4][2] * rhs__[4][2];

//             t1 = bt * rhs__[2][2];
//             t2 = 0.5 * (rhs__[3][2] + rhs__[4][2]);

//             rhs[k][j+2][i][0] = -rhs__[1][2];
//             rhs[k][j+2][i][1] =  rhs__[0][2];
//             rhs[k][j+2][i][2] =  bt * (rhs__[3][2] - rhs__[4][2]);
//             rhs[k][j+2][i][3] = -t1 + t2;
//             rhs[k][j+2][i][4] =  t1 + t2;

//             for (m = 0; m < 5; m++) {
//                 rhs__[m][2] = rhs__[m][1];
//                 rhs__[m][1] = rhs__[m][0];
//             }
//         }

//         t1 = bt * rhs__[2][2];
//         t2 = 0.5 * (rhs__[3][2] + rhs__[4][2]);

//         rhs[k][j+2][i][0] = -rhs__[1][2];
//         rhs[k][j+2][i][1] =  rhs__[0][2];
//         rhs[k][j+2][i][2] =  bt * (rhs__[3][2] - rhs__[4][2]);
//         rhs[k][j+2][i][3] = -t1 + t2;
//         rhs[k][j+2][i][4] =  t1 + t2;
//       }
//     }
//   }
//   if (timeron) timer_stop(t_ysolve);
}

