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
#define Z(z) ((z) + 2) 

//---------------------------------------------------------------------
// compute the right hand side based on exact solution
//---------------------------------------------------------------------
void exact_rhs()
{
  double dtemp[5], xi, eta, zeta, dtpp;
  int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1, z;
  double ue_loc[5][5], buf_loc[5][5], cuf_loc[5], q_loc[5];

  //---------------------------------------------------------------------
  // initialize                                  
  //---------------------------------------------------------------------
  #pragma dvm region targets(HOST) 
  {
    #pragma dvm parallel([k][j][i] on tmp[k][j][i]) private(m)
    for (k = 0; k <= grid_points[2]-1; k++) {
      for (j = 0; j <= grid_points[1]-1; j++) {
        for (i = 0; i <= grid_points[0]-1; i++) {
          for (m = 0; m < 5; m++) {
            forcing[k][j][i][m] = 0.0;
          }
        }
      }
    }
  

  //---------------------------------------------------------------------
  // xi-direction flux differences                      
  //---------------------------------------------------------------------
    #pragma dvm parallel([k][j][i] on forcing[k][j][i][]) private(z,m,zeta,eta,xi,dtemp,dtpp,buf_loc,cuf_loc,q_loc,ue_loc)
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        for (i = 1; i <= grid_points[0]-2; i++) {
          zeta = (double)k * dnzm1;
          eta = (double)j * dnym1;

          for(z = -2; z < 3 ; z++) {
            xi = (double)(i+z) * dnxm1;

            for (m = 0; m < 5; m++) {
              dtemp[m] = ce[m][0] +
                xi  *(ce[m][1] + xi  *(ce[m][4] + xi  *(ce[m][7] + xi  *ce[m][10]))) +
                eta *(ce[m][2] + eta *(ce[m][5] + eta *(ce[m][8] + eta *ce[m][11]))) +
                zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + zeta*ce[m][12])));
            }

          for (m = 0; m < 5; m++)
              ue_loc[Z(z)][m] = dtemp[m];

          dtpp = 1.0 / dtemp[0];
          for (m = 1; m < 5; m++)
              buf_loc[Z(z)][m] = dtpp * dtemp[m];

          cuf_loc[Z(z)] = buf_loc[Z(z)][1] * buf_loc[Z(z)][1];

          buf_loc[Z(z)][0] = cuf_loc[Z(z)] +
                        buf_loc[Z(z)][2]*buf_loc[Z(z)][2] +
                        buf_loc[Z(z)][3]*buf_loc[Z(z)][3];


          q_loc[Z(z)] = 0.5*( buf_loc[Z(z)][1]*ue_loc[Z(z)][1] +
                          buf_loc[Z(z)][2]*ue_loc[Z(z)][2] +
                          buf_loc[Z(z)][3]*ue_loc[Z(z)][3] );
        }

        forcing[k][j][i][0] = forcing[k][j][i][0]
          - tx2    * ( ue_loc[Z( 1)][1] - ue_loc[Z(-1)][1] )
          + dx1tx1 * ( ue_loc[Z( 1)][0]
                  - 2.0*ue_loc[Z( 0)][0]
                  +      ue_loc[Z(-1)][0] );

        forcing[k][j][i][1] = forcing[k][j][i][1] - tx2 * (
            ( ue_loc[Z( 1)][1]*buf_loc[Z( 1)][1]
              + c2*(ue_loc[Z( 1)][4] - q_loc[Z( 1)]) ) -
            ( ue_loc[Z(-1)][1]*buf_loc[Z(-1)][1]
              + c2*(ue_loc[Z(-1)][4] - q_loc[Z(-1)]) ) )
          + xxcon1 * ( buf_loc[Z( 1)][1]
                    - 2.0*buf_loc[Z( 0)][1]
                    +      buf_loc[Z(-1)][1] )
          + dx2tx1 * ( ue_loc[Z( 1)][1]
                    - 2.0*ue_loc[Z( 0)][1]
                    +      ue_loc[Z(-1)][1] );

        forcing[k][j][i][2] = forcing[k][j][i][2] - tx2 * (
            ue_loc[Z( 1)][2]*buf_loc[Z( 1)][1]
          - ue_loc[Z(-1)][2]*buf_loc[Z(-1)][1] )
          + xxcon2 * ( buf_loc[Z( 1)][2]
                    - 2.0*buf_loc[Z( 0)][2]
                    +      buf_loc[Z(-1)][2] )
          + dx3tx1 * ( ue_loc[Z( 1)][2]
                    - 2.0*ue_loc[Z( 0)][2]
                    +      ue_loc[Z(-1)][2] );

        forcing[k][j][i][3] = forcing[k][j][i][3] - tx2 * (
            ue_loc[Z( 1)][3]*buf_loc[Z( 1)][1]
          - ue_loc[Z(-1)][3]*buf_loc[Z(-1)][1] )
          + xxcon2 * ( buf_loc[Z( 1)][3]
                    - 2.0*buf_loc[Z( 0)][3]
                    +      buf_loc[Z(-1)][3] )
          + dx4tx1 * ( ue_loc[Z( 1)][3]
                    - 2.0*ue_loc[Z( 0)][3]
                    +      ue_loc[Z(-1)][3] );

        forcing[k][j][i][4] = forcing[k][j][i][4] - tx2 * (
            buf_loc[Z( 1)][1]*(c1*ue_loc[Z( 1)][4] - c2*q_loc[Z( 1)]) -
            buf_loc[Z(-1)][1]*(c1*ue_loc[Z(-1)][4] - c2*q_loc[Z(-1)]) )
          + 0.5*xxcon3 * ( buf_loc[Z( 1)][0]
                        - 2.0*buf_loc[Z( 0)][0]
                        +      buf_loc[Z(-1)][0] )
          + xxcon4 * ( cuf_loc[Z( 1)]
                    - 2.0*cuf_loc[Z( 0)]
                    +      cuf_loc[Z(-1)] )
          + xxcon5 * ( buf_loc[Z( 1)][4]
                    - 2.0*buf_loc[Z( 0)][4]
                    +      buf_loc[Z(-1)][4] )
          + dx5tx1 * ( ue_loc[Z( 1)][4]
                    - 2.0*ue_loc[Z( 0)][4]
                    +      ue_loc[Z(-1)][4] );

        //---------------------------------------------------------------------
        // Fourth-order dissipation                         
        //---------------------------------------------------------------------
        for (m = 0; m < 5; m++) {
          if (i == 1)
            forcing[k][j][i][m] -= dssp *
              ( 5.0*ue_loc[Z( 0)][m]
              - 4.0*ue_loc[Z( 1)][m]
              +      ue_loc[Z( 2)][m] );
          else if (i == 2)
            forcing[k][j][i][m] -= dssp *
            ( -4.0*ue_loc[Z(-1)][m]
              + 6.0*ue_loc[Z( 0)][m]
              - 4.0*ue_loc[Z( 1)][m]
              +     ue_loc[Z( 2)][m] );
          else if (i == grid_points[0]-3)
            forcing[k][j][i][m] -= dssp *
            ( ue_loc[Z(-2)][m]
              - 4.0*ue_loc[Z(-1)][m]
              + 6.0*ue_loc[Z( 0)][m]
              - 4.0*ue_loc[Z( 1)][m] );
          else if (i == grid_points[0]-2)
            forcing[k][j][i][m] -= dssp *
            ( ue_loc[Z(-2)][m]
              - 4.0*ue_loc[Z(-1)][m]
              + 5.0*ue_loc[Z( 0)][m] );
          else
            forcing[k][j][i][m] -= dssp *
            ( ue_loc[Z(-2)][m]
              - 4.0*ue_loc[Z(-1)][m]
              + 6.0*ue_loc[Z( 0)][m]
              - 4.0*ue_loc[Z( 1)][m]
              +     ue_loc[Z( 2)][m] );
          }
        }
      }
    }
    //---------------------------------------------------------------
    // eta-direction flux differences             
    //---------------------------------------------------------------------
    #pragma dvm parallel([k][j][i] on forcing[k][j][i][]) private(z,m,zeta,eta,xi,dtpp,dtemp,ue_loc,buf_loc,cuf_loc,q_loc)
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        for (i = 1; i <= grid_points[0]-2; i++) {
          zeta = (double)k * dnzm1;
          xi = (double)i * dnxm1;

          for(z = -2; z < 3 ; z++) {
            eta = (double)(j+z) * dnym1;

            for (m = 0; m < 5; m++) {
              dtemp[m] = ce[m][0] +
                xi  *(ce[m][1] + xi  *(ce[m][4] + xi  *(ce[m][7] + xi  *ce[m][10]))) +
                eta *(ce[m][2] + eta *(ce[m][5] + eta *(ce[m][8] + eta *ce[m][11]))) +
                zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + zeta*ce[m][12])));
            }

            for (m = 0; m < 5; m++) 
              ue_loc[Z(z)][m] = dtemp[m];
            
            dtpp = 1.0/dtemp[0];

            for (m = 1; m < 5; m++)
              buf_loc[Z(z)][m] = dtpp * dtemp[m];
            

            cuf_loc[Z(z)]    = buf_loc[Z(z)][2] * buf_loc[Z(z)][2];
            buf_loc[Z(z)][0] = cuf_loc[Z(z)] + buf_loc[Z(z)][1] * buf_loc[Z(z)][1] + buf_loc[Z(z)][3] * buf_loc[Z(z)][3];
            q_loc[Z(z)] = 0.5*(buf_loc[Z(z)][1]*ue_loc[Z(z)][1] + buf_loc[Z(z)][2]*ue_loc[Z(z)][2] +
                        buf_loc[Z(z)][3]*ue_loc[Z(z)][3]);
          }

          forcing[k][j][i][0] = forcing[k][j][i][0]
            - ty2    * ( ue_loc[Z( 1)][2] - ue_loc[Z(-1)][2] )
            + dy1ty1 * ( ue_loc[Z( 1)][0]
                      - 2.0*ue_loc[Z( 0)][0]
                      +      ue_loc[Z(-1)][0] );

          forcing[k][j][i][1] = forcing[k][j][i][1] - ty2 * (
              ue_loc[Z( 1)][1]*buf_loc[Z( 1)][2]
            - ue_loc[Z(-1)][1]*buf_loc[Z(-1)][2] )
            + yycon2 * ( buf_loc[Z( 1)][1]
                      - 2.0*buf_loc[Z( 0)][1]
                      +      buf_loc[Z(-1)][1] )
            + dy2ty1 * ( ue_loc[Z( 1)][1]
                      - 2.0*ue_loc[Z( 0)][1]
                      +      ue_loc[Z(-1)][1] );

          forcing[k][j][i][2] = forcing[k][j][i][2] - ty2 * (
              ( ue_loc[Z( 1)][2]*buf_loc[Z( 1)][2]
                + c2*(ue_loc[Z( 1)][4] - q_loc[Z( 1)]) ) -
              ( ue_loc[Z(-1)][2]*buf_loc[Z(-1)][2]
                + c2*(ue_loc[Z(-1)][4] - q_loc[Z(-1)]) ) )
            + yycon1 * ( buf_loc[Z( 1)][2]
                      - 2.0*buf_loc[Z( 0)][2]
                      +      buf_loc[Z(-1)][2] )
            + dy3ty1 * ( ue_loc[Z( 1)][2]
                      - 2.0*ue_loc[Z( 0)][2]
                      +      ue_loc[Z(-1)][2] );

          forcing[k][j][i][3] = forcing[k][j][i][3] - ty2 * (
              ue_loc[Z( 1)][3]*buf_loc[Z( 1)][2]
            - ue_loc[Z(-1)][3]*buf_loc[Z(-1)][2] )
            + yycon2 * ( buf_loc[Z( 1)][3]
                      - 2.0*buf_loc[Z( 0)][3]
                      +      buf_loc[Z(-1)][3] )
            + dy4ty1 * ( ue_loc[Z( 1)][3]
                      - 2.0*ue_loc[Z( 0)][3]
                      +      ue_loc[Z(-1)][3] );

          forcing[k][j][i][4] = forcing[k][j][i][4] - ty2 * (
              buf_loc[Z( 1)][2]*(c1*ue_loc[Z( 1)][4] - c2*q_loc[Z( 1)]) -
              buf_loc[Z(-1)][2]*(c1*ue_loc[Z(-1)][4] - c2*q_loc[Z(-1)]) )
            + 0.5*yycon3 * ( buf_loc[Z( 1)][0]
                          - 2.0*buf_loc[Z( 0)][0]
                          +      buf_loc[Z(-1)][0] )
            + yycon4 * ( cuf_loc[Z( 1)]
                      - 2.0*cuf_loc[Z( 0)]
                      +      cuf_loc[Z(-1)] )
            + yycon5 * ( buf_loc[Z( 1)][4]
                      - 2.0*buf_loc[Z( 0)][4]
                      +      buf_loc[Z(-1)][4] )
            + dy5ty1 * ( ue_loc[Z( 1)][4]
                      - 2.0*ue_loc[Z( 0)][4]
                      +      ue_loc[Z(-1)][4] );

          //-----------------------------------------------------------------
          // Fourth-order dissipation in eta (y) direction
          //-----------------------------------------------------------------
          for (m = 0; m < 5; m++) {
            if (j == 1)
              forcing[k][j][i][m] -= dssp *
                ( 5.0*ue_loc[Z( 0)][m]
                - 4.0*ue_loc[Z( 1)][m]
                +      ue_loc[Z( 2)][m] );
            else if (j == 2)
              forcing[k][j][i][m] -= dssp *
                ( -4.0*ue_loc[Z(-1)][m]
                +  6.0*ue_loc[Z( 0)][m]
                -  4.0*ue_loc[Z( 1)][m]
                +      ue_loc[Z( 2)][m] );
            else if (j == grid_points[1]-3)
              forcing[k][j][i][m] -= dssp *
                ( ue_loc[Z(-2)][m]
                - 4.0*ue_loc[Z(-1)][m]
                + 6.0*ue_loc[Z( 0)][m]
                - 4.0*ue_loc[Z( 1)][m] );
            else if (j == grid_points[1]-2)
              forcing[k][j][i][m] -= dssp *
                ( ue_loc[Z(-2)][m]
                - 4.0*ue_loc[Z(-1)][m]
                + 5.0*ue_loc[Z( 0)][m] );
            else
              forcing[k][j][i][m] -= dssp *
                ( ue_loc[Z(-2)][m]
                - 4.0*ue_loc[Z(-1)][m]
                + 6.0*ue_loc[Z( 0)][m]
                - 4.0*ue_loc[Z( 1)][m]
                +      ue_loc[Z( 2)][m] );
          }
        }
      }
    }


    //---------------------------------------------------------------------
    // zeta-direction flux differences                      
    //---------------------------------------------------------------------
    #pragma dvm parallel([k][j][i] on forcing[k][j][i][]) private(z,m,zeta,eta,xi,dtpp,dtemp,ue_loc,buf_loc,cuf_loc,q_loc)
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        for (i = 1; i <= grid_points[0]-2; i++) {
          eta = (double)j * dnym1;
          xi  = (double)i * dnxm1;

          for (z = -2; z <= 2; z++) {
            zeta = (double)(k+z) * dnzm1;

            for (m = 0; m < 5; m++) {
              dtemp[m] = ce[m][0] +
                  xi  *(ce[m][1] + xi  *(ce[m][4] + xi  *(ce[m][7] + xi  *ce[m][10]))) +
                  eta *(ce[m][2] + eta *(ce[m][5] + eta *(ce[m][8] + eta *ce[m][11]))) +
                  zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + zeta*ce[m][12])));
            }

            for (m = 0; m < 5; m++)
              ue_loc[Z(z)][m] = dtemp[m];

            dtpp = 1.0 / dtemp[0];
            for (m = 1; m < 5; m++)
              buf_loc[Z(z)][m] = dtpp * dtemp[m];

            cuf_loc[Z(z)] = buf_loc[Z(z)][3] * buf_loc[Z(z)][3];

            buf_loc[Z(z)][0] = cuf_loc[Z(z)]
                            + buf_loc[Z(z)][1]*buf_loc[Z(z)][1]
                            + buf_loc[Z(z)][2]*buf_loc[Z(z)][2];

            q_loc[Z(z)] = 0.5 * ( buf_loc[Z(z)][1]*ue_loc[Z(z)][1]
                                + buf_loc[Z(z)][2]*ue_loc[Z(z)][2]
                                + buf_loc[Z(z)][3]*ue_loc[Z(z)][3] );
          }

          forcing[k][j][i][0] = forcing[k][j][i][0]
            - tz2    * ( ue_loc[Z( 1)][3] - ue_loc[Z(-1)][3] )
            + dz1tz1 * ( ue_loc[Z( 1)][0]
                      - 2.0*ue_loc[Z( 0)][0]
                      +      ue_loc[Z(-1)][0] );

          forcing[k][j][i][1] = forcing[k][j][i][1] - tz2 * (
              ue_loc[Z( 1)][1]*buf_loc[Z( 1)][3]
            - ue_loc[Z(-1)][1]*buf_loc[Z(-1)][3] )
            + zzcon2 * ( buf_loc[Z( 1)][1]
                      - 2.0*buf_loc[Z( 0)][1]
                      +      buf_loc[Z(-1)][1] )
            + dz2tz1 * ( ue_loc[Z( 1)][1]
                      - 2.0*ue_loc[Z( 0)][1]
                      +      ue_loc[Z(-1)][1] );

          forcing[k][j][i][2] = forcing[k][j][i][2] - tz2 * (
              ue_loc[Z( 1)][2]*buf_loc[Z( 1)][3]
            - ue_loc[Z(-1)][2]*buf_loc[Z(-1)][3] )
            + zzcon2 * ( buf_loc[Z( 1)][2]
                      - 2.0*buf_loc[Z( 0)][2]
                      +      buf_loc[Z(-1)][2] )
            + dz3tz1 * ( ue_loc[Z( 1)][2]
                      - 2.0*ue_loc[Z( 0)][2]
                      +      ue_loc[Z(-1)][2] );

          forcing[k][j][i][3] = forcing[k][j][i][3] - tz2 * (
              ( ue_loc[Z( 1)][3]*buf_loc[Z( 1)][3]
                + c2*(ue_loc[Z( 1)][4] - q_loc[Z( 1)]) ) -
              ( ue_loc[Z(-1)][3]*buf_loc[Z(-1)][3]
                + c2*(ue_loc[Z(-1)][4] - q_loc[Z(-1)]) ) )
            + zzcon1 * ( buf_loc[Z( 1)][3]
                      - 2.0*buf_loc[Z( 0)][3]
                      +      buf_loc[Z(-1)][3] )
            + dz4tz1 * ( ue_loc[Z( 1)][3]
                      - 2.0*ue_loc[Z( 0)][3]
                      +      ue_loc[Z(-1)][3] );

          forcing[k][j][i][4] = forcing[k][j][i][4] - tz2 * (
              buf_loc[Z( 1)][3]*(c1*ue_loc[Z( 1)][4] - c2*q_loc[Z( 1)]) -
              buf_loc[Z(-1)][3]*(c1*ue_loc[Z(-1)][4] - c2*q_loc[Z(-1)]) )
            + 0.5*zzcon3 * ( buf_loc[Z( 1)][0]
                          - 2.0*buf_loc[Z( 0)][0]
                          +      buf_loc[Z(-1)][0] )
            + zzcon4 * ( cuf_loc[Z( 1)]
                      - 2.0*cuf_loc[Z( 0)]
                      +      cuf_loc[Z(-1)] )
            + zzcon5 * ( buf_loc[Z( 1)][4]
                      - 2.0*buf_loc[Z( 0)][4]
                      +      buf_loc[Z(-1)][4] )
            + dz5tz1 * ( ue_loc[Z( 1)][4]
                      - 2.0*ue_loc[Z( 0)][4]
                      +      ue_loc[Z(-1)][4] );

          //-----------------------------------------------------------------
          // Fourth-order dissipation in zeta (z) direction
          //-----------------------------------------------------------------
          for (m = 0; m < 5; m++) {
            if (k == 1)
              forcing[k][j][i][m] -= dssp *
                ( 5.0*ue_loc[Z( 0)][m]
                - 4.0*ue_loc[Z( 1)][m]
                +      ue_loc[Z( 2)][m] );
            else if (k == 2)
              forcing[k][j][i][m] -= dssp *
                ( -4.0*ue_loc[Z(-1)][m]
                +  6.0*ue_loc[Z( 0)][m]
                -  4.0*ue_loc[Z( 1)][m]
                +      ue_loc[Z( 2)][m] );
            else if (k == grid_points[2]-3)
              forcing[k][j][i][m] -= dssp *
                ( ue_loc[Z(-2)][m]
                - 4.0*ue_loc[Z(-1)][m]
                + 6.0*ue_loc[Z( 0)][m]
                - 4.0*ue_loc[Z( 1)][m] );
            else if (k == grid_points[2]-2)
              forcing[k][j][i][m] -= dssp *
                ( ue_loc[Z(-2)][m]
                - 4.0*ue_loc[Z(-1)][m]
                + 5.0*ue_loc[Z( 0)][m] );
            else
              forcing[k][j][i][m] -= dssp *
                ( ue_loc[Z(-2)][m]
                - 4.0*ue_loc[Z(-1)][m]
                + 6.0*ue_loc[Z( 0)][m]
                - 4.0*ue_loc[Z( 1)][m]
                +      ue_loc[Z( 2)][m] );
          }
        }
      }
    }
  //---------------------------------------------------------------------
  // now change the sign of the forcing function, 
  //---------------------------------------------------------------------
    #pragma dvm parallel([k][j][i] on forcing[k][j][i][]) private(m)
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        for (i = 1; i <= grid_points[0]-2; i++) {
          for (m = 0; m < 5; m++) {
            forcing[k][j][i][m] = -1.0 * forcing[k][j][i][m];
          }
        }
      }
    }
  }
}

