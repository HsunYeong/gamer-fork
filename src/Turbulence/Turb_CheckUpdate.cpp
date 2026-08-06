#include "GAMER.h"
#include <algorithm>

#ifdef TURBULENCE

//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_CheckUpdate
// Description :  Check TimeNew at lv 0, and update OUphase and AccTable
//
// Note        :  1. Invoked by EvolveLevel()
//
// Parameter   :  TimeNew : TimeNew for base level
//
// Return      :  Turb->TimeLast, Turb->TimeNext, Turb->OUphase, Turb->AccTable
//-------------------------------------------------------------------------------------------------------
void Turb_CheckUpdate( const double TimeNew )
{
   if ( !TURB_ACTIVATE )   return;

   int hasUpdate = 0;

   while ( TimeNew > Turb->TimeNext )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "TimeNew ( %13.7e ) > Turbulence TimeNext ( %13.7e ): Update turbulence pattern ...", TimeNew, Turb->TimeNext );

      double dt  = Turb->dt;
      double tau = Turb->Tdecay;
      double var = Turb->OUvar;
      double coeff1 = exp(-dt/tau);
      double coeff2 = sqrt( 1 - SQR(coeff1) );

//    swap last and next OUphases
      std::swap( Turb->OUphase[0], Turb->OUphase[1] );

//    construct OU phase vector
      for (int n = 0; n < Turb->NMode; ++n)
      {
         double kk       = 0;
         double k_dot_Nr = 0;
         double k_dot_Ni = 0;
         double Nr[3], Ni[3];

         for (int d = 0; d < 3; ++d)
         {
//          get random number Nr and Ni
            Turb_GetRNG( Nr[d], Ni[d], Turb->RSeed, var );

            kk       += SQR( Turb->Kmode[d][n] );
            k_dot_Nr += Turb->Kmode[d][n]*Nr[d];
            k_dot_Ni += Turb->Kmode[d][n]*Ni[d];
         }

         for (int d = 0; d < 3; ++d)
         {
//          Helmholtz decomposition
            Nr[d] = TURB_ZETA*Nr[d] + (1 - 2*TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Nr/kk;
            Ni[d] = TURB_ZETA*Ni[d] + (1 - 2*TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Ni/kk;

//          Update OU phases to time_new
            Turb->OUphase[1][2*3*n+2*d  ] = coeff1 * Turb->OUphase[0][2*3*n+2*d  ] + coeff2 * Nr[d];
            Turb->OUphase[1][2*3*n+2*d+1] = coeff1 * Turb->OUphase[0][2*3*n+2*d+1] + coeff2 * Ni[d];
         }

      } // for (int n = 0; n < Turb->NMode; ++n)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, " done\n" );

//    update turbulence time
      Turb->TimeLast =  Turb->TimeNext;
      Turb->TimeNext += Turb->dt;

      hasUpdate += 1;

   } // while ( Time[0] > Turb->Time )

// check
   if ( TimeNew > Turb->TimeNext || TimeNew < Turb->TimeLast )
      Aux_Error( ERROR_INFO, "TimeNew of lv 0 outside turbulence time range ( TimeNew %24.17e, Turb->TimeLast %24.17e, Turb->TimeNext %24.17e ) !!\n",
                              TimeNew, Turb->TimeLast, Turb->TimeNext );

// update turb acc table
   if ( hasUpdate > 0 )
   {
//    if there is only one OU update
      if ( hasUpdate == 1 )
      {
//       swap last and next AccTable
         std::swap( Turb->AccTable[0], Turb->AccTable[1] );

//       update new table
         Turb_FillinTable(1);

      } // if ( hasUpdate == 1 )
      else
      {
//       fillin both tables
//       this should be prevented in general by choosing a large enough turbulence dt
         Turb_FillinTable(0);
      } // else
   } // if ( hasUpdate )

} // FUNCTION : Turb_CheckUpdate


#endif // #ifdef TURBULENCE
