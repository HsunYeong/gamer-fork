#include "GAMER.h"

#ifdef TURBULENCE

void Turb_Aux_CheckUpdate()
{
   bool hasUpdate = false;

   while ( Time[0] > Turb->Time )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "Time ( %13.7e ) > Turbulence time ( %13.7e ): Update turbulence pattern ...\n", Time[0], Turb->Time );

      double dt  = Turb->dt;
      double tau = Turb->Tdecay;
      double var = Turb->OUvar;
      double coeff1 = exp(-dt/tau);
      double coeff2 = var * sqrt( 1 - SQR(coeff1) );

//    Construct OU phase vector
      for (int n = 0; n < Turb->NMode; ++n)
      {
         double kk       = 0;
         double k_dot_Nr = 0;
         double k_dot_Ni = 0;
         double Nr[3], Ni[3];

         for (int d = 0; d < 3; ++d)
         {
//          get random number Nr and Ni
            Turb_GetRNG( Nr[d], Ni[d], Turb->RSeed );

            kk       += SQR( Turb->Kmode[d][n] );
            k_dot_Nr += Turb->Kmode[d][n]*Nr[d];
            k_dot_Ni += Turb->Kmode[d][n]*Ni[d];
         }

         for (int d = 0; d < 3; ++d)
         {
//          Helmholtz decomposition
            Nr[d] = TURB_ZETA*Nr[d] + (1 - 2*TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Nr/kk;
            Ni[d] = TURB_ZETA*Ni[d] + (1 - 2*TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Ni/kk;

//          Update OU phases
            Turb->OUphase[2*3*n+2*d  ] = coeff1 * Turb->OUphase[2*3*n+2*d  ] + coeff2 * Nr[d];
            Turb->OUphase[2*3*n+2*d+1] = coeff1 * Turb->OUphase[2*3*n+2*d+1] + coeff2 * Ni[d];
         }

      } // for (int n = 0; n < Turb->NMode; ++n)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, " done\n" );

//    update turbulence time
      Turb->Time += dt;
      hasUpdate = true;

   } // while ( Time[0] > Turb->Time )

// update turb acc table
   if ( hasUpdate )
   {
      const long NPoint = TURB_TABLE_SIZE + 1;
      const double dh   = BOX_SIZE/TURB_TABLE_SIZE;

#     pragma omp parallel for schedule( runtime )
      for (int k = 0; k < NPoint; k++)  {
      for (int j = 0; j < NPoint; j++)  {
      for (int i = 0; i < NPoint; i++)  {
         double Acc[3] = {0, 0, 0};
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

            for (int d = 0; d < 3; d++)
            {
               Acc[d] += amp*( Turb->OUphase[2*3*n+2*d]*real - Turb->OUphase[2*3*n+2*d+1]*imag );
            }
         } // for (int n = 0; n < Turb->NMode; n++)
         for (int d = 0; d < 3; d++)
         {
            TurbAccTable[d][idx] = Acc[d]*TURB_AMPL_FACTOR;
         }
      }}} // for i, j, k

   } // if ( hasUpdate )

} // FUNCTION : Turb_Aux_CheckUpdate


#endif // #ifdef TURBULENCE
