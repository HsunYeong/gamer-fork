#include "GAMER.h"

#ifdef TURBULENCE


//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_FillinTable
// Description :  Fillin turbulence acc table
//
// Note        :  1. Invoked by Turb_Init(), Turb_CheckUpdate()
//
// Parameter   :  time_idx : 0 -> Fill in both table 0 and table 1
//                           1 -> Fill in only table 1
//
// Return      :  Turb->AccTable
//-------------------------------------------------------------------------------------------------------
void Turb_FillinTable(int time_idx)
{
   if ( !TURB_ACTIVATE )   return;

   const long NPoint = TURB_TABLE_SIZE + 1;

#  pragma omp parallel for schedule( runtime )
   for (int k = 0; k < NPoint; k++)  {
   for (int j = 0; j < NPoint; j++)  {
   for (int i = 0; i < NPoint; i++)  {
      double Acc[2][3] = {{0}};
      long  idx = IDX321( i, j, k, NPoint, NPoint );

      for (int n = 0; n < Turb->NMode; n++)
      {
         double sinx = Turb->Sin[0][ n*NPoint + i ];
         double cosx = Turb->Cos[0][ n*NPoint + i ];
         double siny = Turb->Sin[1][ n*NPoint + j ];
         double cosy = Turb->Cos[1][ n*NPoint + j ];
         double sinz = Turb->Sin[2][ n*NPoint + k ];
         double cosz = Turb->Cos[2][ n*NPoint + k ];
         double amp  = Turb->Amplitude[n];

         double real = ( cosx*cosy - sinx*siny ) * cosz - ( sinx*cosy + cosx*siny ) * sinz;
         double imag = ( cosy*sinz + siny*cosz ) * cosx + ( cosy*cosz - siny*sinz ) * sinx;

         for (int t=time_idx; t<2; t++)
         for (int d=0; d<3; d++)
            Acc[t][d] += amp*( Turb->OUphase[t][2*3*n+2*d]*real - Turb->OUphase[t][2*3*n+2*d+1]*imag );

      } // for (int n = 0; n < Turb->NMode; n++)

      for (int t=time_idx; t<2; t++)
      for (int d=0; d<3; d++)
         Turb->AccTable[t][idx][d] = (real)Acc[t][d]*TURB_AMPL_FACTOR;

   }}} // for i, j, k

}
#endif // #ifdef TURBULENCE
