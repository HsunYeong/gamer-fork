#include "GAMER.h"


//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_FillinTable
// Description :  Fillin h_Turb_AccTable host array
//
// Note        :  1. Invoked by Turb_Init(), Turb_CheckUpdate()
//
// Parameter   :  IdxTable : Turbulence table index-> IdxLast or IdxNext
//
// Return      :  h_Turb_AccTable
//-------------------------------------------------------------------------------------------------------
void Turb_FillinTable( int IdxTable )
{
   const long NPoint = TURB_TABLE_SIZE + 1;

#  pragma omp parallel for schedule( runtime )
   for (int k = 0; k < NPoint; k++)  {
   for (int j = 0; j < NPoint; j++)  {
   for (int i = 0; i < NPoint; i++)  {
      double Acc[3] = {0};
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

         for (int d=0; d<3; d++)
            Acc[d] += amp*( Turb->OUphase[IdxTable][2*3*n+2*d]*real - Turb->OUphase[IdxTable][2*3*n+2*d+1]*imag );

      } // for (int n = 0; n < Turb->NMode; n++)

      for (int d=0; d<3; d++)
         h_Turb_AccTable[IdxTable][3*idx + d] = (real)Acc[d];

   }}} // for i, j, k

}

