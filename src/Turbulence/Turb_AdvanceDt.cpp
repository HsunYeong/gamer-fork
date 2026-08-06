#include "GAMER.h"

#ifdef TURBULENCE




//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_AdvanceDt
// Description :  Update momentum and total energy by turbulence acceleration table
//
// Note        :  1. Invoked by EvolveLevel()
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach
//                dt           : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//                SaveSg       : Sandglass to store the updated data
//-------------------------------------------------------------------------------------------------------
void Turb_AdvanceDt( const int lv, const double TimeNew, const double dt, const int SaveSg )
{
   const double dh = amr->dh[lv];
   const double _dh_table = TURB_TABLE_SIZE/BOX_SIZE;
   const long NPoint = TURB_TABLE_SIZE + 1;
   const long didx_x = 1;
   const long didx_y = NPoint;
   const long didx_z = SQR( NPoint );
   const double tfrac = (TimeNew - Turb->TimeLast) / Turb->dt;
   const double one_m_tfrac = 1 - tfrac;

   if ( Turb == NULL )
      Aux_Error( ERROR_INFO, "Turb == NULL at rank %d!!\n", MPI_Rank );

#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + 0.5*dh;
      for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + 0.5*dh;
      for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + 0.5*dh;
         real dx    = real( x * _dh_table );
         real dy    = real( y * _dh_table );
         real dz    = real( z * _dh_table );
         int  idx_x = int( dx );
         int  idx_y = int( dy );
         int  idx_z = int( dz );
         long idx0  = long( idx_x*didx_x + idx_y*didx_y ) + (long)idx_z*didx_z;

//       compute acc by trilinear interpolation
         real weight_xL, weight_yL, weight_zL;
         real weight_xR, weight_yR, weight_zR;
         real Acc[2][3] = {{0}};

         weight_xR = dx - (real)idx_x;
         weight_yR = dy - (real)idx_y;
         weight_zR = dz - (real)idx_z;
         weight_xL = 1.0 - weight_xR;
         weight_yL = 1.0 - weight_yR;
         weight_zL = 1.0 - weight_zR;

         for (int t=0; t<2; t++)
         for (int d=0; d<3; d++)
         {
            Acc[t][d] = Turb->AccTable[t][ idx0                            ][d] * weight_xL * weight_yL * weight_zL +
                        Turb->AccTable[t][ idx0 + didx_x                   ][d] * weight_xR * weight_yL * weight_zL +
                        Turb->AccTable[t][ idx0          + didx_y          ][d] * weight_xL * weight_yR * weight_zL +
                        Turb->AccTable[t][ idx0                   + didx_z ][d] * weight_xL * weight_yL * weight_zR +
                        Turb->AccTable[t][ idx0 + didx_x + didx_y          ][d] * weight_xR * weight_yR * weight_zL +
                        Turb->AccTable[t][ idx0          + didx_y + didx_z ][d] * weight_xL * weight_yR * weight_zR +
                        Turb->AccTable[t][ idx0 + didx_x          + didx_z ][d] * weight_xR * weight_yL * weight_zR +
                        Turb->AccTable[t][ idx0 + didx_x + didx_y + didx_z ][d] * weight_xR * weight_yR * weight_zR;
         }

         const double dens = amr->patch[SaveSg][lv][PID]->fluid[DENS][k][j][i];
         const double velx = amr->patch[SaveSg][lv][PID]->fluid[MOMX][k][j][i] / dens;
         const double vely = amr->patch[SaveSg][lv][PID]->fluid[MOMY][k][j][i] / dens;
         const double velz = amr->patch[SaveSg][lv][PID]->fluid[MOMZ][k][j][i] / dens;

         const double dMomX = dens * dt * ( one_m_tfrac*Acc[0][0] + tfrac*Acc[1][0] );
         const double dMomY = dens * dt * ( one_m_tfrac*Acc[0][1] + tfrac*Acc[1][1] );
         const double dMomZ = dens * dt * ( one_m_tfrac*Acc[0][2] + tfrac*Acc[1][2] );
         const double dE    = velx*dMomX + vely*dMomY + velz*dMomZ  + (dMomX*dMomX + dMomY*dMomY + dMomZ*dMomZ)/(2.0*dens);

         amr->patch[SaveSg][lv][PID]->fluid[MOMX][k][j][i] += (real)dMomX;
         amr->patch[SaveSg][lv][PID]->fluid[MOMY][k][j][i] += (real)dMomY;
         amr->patch[SaveSg][lv][PID]->fluid[MOMZ][k][j][i] += (real)dMomZ;
         amr->patch[SaveSg][lv][PID]->fluid[ENGY][k][j][i] += (real)dE;

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Turb_AdvanceDt



#endif // #ifdef TURBULENCE
