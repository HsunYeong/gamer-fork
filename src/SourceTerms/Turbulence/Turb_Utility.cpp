#include "GAMER.h"

#if ( MODEL == HYDRO )
void   Turb_Init();
void   Turb_End();
void   Turb_GetRNG( double& a, double& b, int& Seed, const double OUvar );
void   Turb_FillinTable( int IdxTable );
double Turb_ran1s( int& Seed );

/********************************************************************************************************
Turbulence structure:

1. Stores turbulence Fourier modes, amplitude, random phases, and sin, cos Fourier basis

2. Use Helmholtz decomposition to separate compressive and solenoial components

3. Update random phases by Ornstein-Uhlenbeck process

4. References: Federrath et al. (2010), A&A 512, A81 (https://doi.org/10.1051/0004-6361/200912437)
               TurbGen (https://github.com/chfeder/turbulence_generator)

*********************************************************************************************************/



//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_Init
// Description :  Initialize turbulence data structure
//
// Note        :  1. Invoked by Src_Init_Turbulence()
//                2. Don't fill in AccTable here since global arrays are not yet initialized.
//                3. When restart, load OU phases, times, and random seed.
//
// Parameter   :  None
//
// Return      :  Turb
//-------------------------------------------------------------------------------------------------------
void Turb_Init()
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// initialize turbulence field
   Turb = new Turbulence_t;

// assign values to structure members
   Turb->Tdecay   = BOX_SIZE/TURB_KDRIV/TURB_VEL;
   Turb->dt       = Turb->Tdecay/TURB_UPDATE_STEP;
   Turb->RSeed    = TURB_RSEED_INIT;

   const double ZetaNorm = sqrt(3.0)/sqrt(1.0 - 2.0*TURB_ZETA + 3.0*SQR(TURB_ZETA));
   const double EnergyInputRate = CUBE(TURB_AMPL_FACTOR*0.15*TURB_VEL)/BOX_SIZE;

// OUvar ~ a_rms
   Turb->OUvar = sqrt( EnergyInputRate/Turb->Tdecay );

// initialize k modes
   double kmin   = (TURB_KMIN - __DBL_EPSILON__) * 2*M_PI / BOX_SIZE;
   double kmax   = (TURB_KMAX + __DBL_EPSILON__) * 2*M_PI / BOX_SIZE;

   if ( kmax < kmin )
      Aux_Error( ERROR_INFO, "Turbulence: kmax ( %13.7e ) < kmin ( %13.7e )!!\n", kmax, kmin );

   double kmid   = 0.5*(kmin + kmax);
   int    Nmax   = 2*TURB_KMAX + 1;
   int    nmodes = 0;
   double kmodes[Nmax];
   for (int i = 0; i < Nmax; ++ i)
      kmodes[i] = 2*M_PI/BOX_SIZE * (i - (Nmax - 1)/2.0);

   for (int k = 0; k < Nmax; ++k) {
   for (int j = 0; j < Nmax; ++j) {
   for (int i = 0; i < Nmax; ++i) {

      double kmag = sqrt( SQR(kmodes[i]) + SQR(kmodes[j]) + SQR(kmodes[k]) );
      if ( kmag >= kmin && kmag <= kmax ) nmodes++;

   }}}

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   initialize %d turbulence modes\n", nmodes );

   Turb->NMode = nmodes;

   if ( Turb->Amplitude == NULL)  Turb->Amplitude = new double [nmodes];

   for (int i = 0; i < 3; ++i)
   {
      if ( Turb->Kmode[i] == NULL ) Turb->Kmode[i] = new double [nmodes];
   }

// get amplitude of each mode
   int n = 0;
   for (int k = 0; k < Nmax; ++k)
   for (int j = 0; j < Nmax; ++j)
   for (int i = 0; i < Nmax; ++i)
   {
      double amp = 0;
      double kmag = sqrt( SQR(kmodes[i]) + SQR(kmodes[j]) + SQR(kmodes[k]) );
      if ( kmag >= kmin && kmag <= kmax ) {
//       constant
         if ( TURB_SPEC_FORM == 0)
            amp = 1.0*kmin/kmag;
//       parabolic
         else if ( TURB_SPEC_FORM == 1 )
            amp = sqrt( fabs(-4 * SQR( (kmag - kmid)/(kmax - kmin) ) + 1) )*kmid/kmag;
//       power law
         else if ( TURB_SPEC_FORM == 2 )
            amp = sqrt( pow(kmag/kmin, TURB_POW) )*kmin/kmag;
         else
            Aux_Error( ERROR_INFO, "Unknown TURB_SPEC_FORM = %d!!\n", TURB_SPEC_FORM );

         Turb->Kmode[0][n] = kmodes[i];
         Turb->Kmode[1][n] = kmodes[j];
         Turb->Kmode[2][n] = kmodes[k];

         Turb->Amplitude[n] = amp*2*ZetaNorm;
         n++;
      } // if ( kmag >= kmin && kmag <= kmax )
   } // for i, j, k

// print turbulence information
   if ( MPI_Rank == 0 )
   {
       Aux_Message( stdout, "Turbulence parameters:\n" );
       Aux_Message( stdout, "   velocity dispersion    = %13.7e:\n", TURB_VEL         );
       Aux_Message( stdout, "   amplitute factor       = %13.7e:\n", TURB_AMPL_FACTOR );
       Aux_Message( stdout, "   energy injection rate  = %13.7e:\n", EnergyInputRate  );
       Aux_Message( stdout, "   kmin                   = %13.7e:\n", kmin             );
       Aux_Message( stdout, "   kmax                   = %13.7e:\n", kmax             );
       Aux_Message( stdout, "   correlation time       = %13.7e:\n", Turb->Tdecay     );
       Aux_Message( stdout, "   update pattern dt      = %13.7e:\n", Turb->dt         );
       Aux_Message( stdout, "   OU variance            = %13.7e:\n", Turb->OUvar      );
       Aux_Message( stdout, "   solenoidal weight norm = %13.7e:\n", ZetaNorm         );
       Aux_Message( stdout, "\n");

      for (int n = 0; n < Turb->NMode; ++n)
      {
         Aux_Message( stdout, "    mode = %3d, amplitude = %13.7e\n", n, Turb->Amplitude[n] );
      }
   }

// initialize OU noise
   for (int t = 0; t < 2; t++)
      if ( Turb->OUphase[t] == NULL ) Turb->OUphase[t] = new double [6*nmodes];

// when restart, load turbulence field
   if ( OPT__INIT == INIT_BY_RESTART && !OPT__RESTART_RESET && !TURB_RESET )
   {
//    read random seed, step, OUarr size from restart file

//    read OUphase, check size

      Aux_Error( ERROR_INFO, "Turbulence from restart is not supported yet !!\n" );
   }
   else
   {
//    loop through two sets
      for (int t = 0; t < 2 ; t++)
      {
//       construct OU phase vector
         for (int n = 0; n < Turb->NMode; ++n)
         {
            double kk       = 0;
            double k_dot_Nr = 0;
            double k_dot_Ni = 0;
            double Nr[3], Ni[3];
            for (int d = 0; d < 3; ++d)
            {
//             get random number Nr and Ni
               Turb_GetRNG( Nr[d], Ni[d], Turb->RSeed, Turb->OUvar );

               kk       += SQR( Turb->Kmode[d][n] );
               k_dot_Nr += Turb->Kmode[d][n]*Nr[d];
               k_dot_Ni += Turb->Kmode[d][n]*Ni[d];
            }

//          Helmholtz decomposition
            for (int d = 0; d < 3; ++d)
            {
               Turb->OUphase[t][2*3*n+2*d  ] = TURB_ZETA*Nr[d] + (1 - 2*TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Nr/kk;
               Turb->OUphase[t][2*3*n+2*d+1] = TURB_ZETA*Ni[d] + (1 - 2*TURB_ZETA)*Turb->Kmode[d][n]*k_dot_Ni/kk;
            }
         } // for (int n = 0; n < Turb->NMode; ++n)
      } // for t

//    perform Ornstein-Uhlenbeck process to update OUphase[1]
      double coeff1 = exp( -Turb->dt/Turb->Tdecay );
      double coeff2 = sqrt( 1 - SQR(coeff1) );
      for (int n = 0; n < Turb->NMode; ++n)
      {
         for (int d = 0; d < 3; ++d)
         {
            Turb->OUphase[1][2*3*n+2*d  ] = coeff1 * Turb->OUphase[0][2*3*n+2*d  ] + coeff2 * Turb->OUphase[1][2*3*n+2*d  ];
            Turb->OUphase[1][2*3*n+2*d+1] = coeff1 * Turb->OUphase[0][2*3*n+2*d+1] + coeff2 * Turb->OUphase[1][2*3*n+2*d+1];
         }
      }
   }

// initialize acc table, store values on box corner
   const long NPoint = TURB_TABLE_SIZE + 1;
   const double dh   = BOX_SIZE/TURB_TABLE_SIZE;

   for (int d = 0; d < 3; ++d)
   {
      Turb->Sin[d] = new double [ NPoint*Turb->NMode ];
      Turb->Cos[d] = new double [ NPoint*Turb->NMode ];
   }

#  pragma omp parallel for schedule( runtime )
   for (int n = 0; n < Turb->NMode; n++)
   {
      for (int i = 0; i < TURB_TABLE_SIZE; i++)
      {
         Turb->Sin[0][ n*NPoint + i ] = sin( Turb->Kmode[0][n]*i*dh );
         Turb->Cos[0][ n*NPoint + i ] = cos( Turb->Kmode[0][n]*i*dh );
      }
      for (int j = 0; j < TURB_TABLE_SIZE; j++)
      {
         Turb->Sin[1][ n*NPoint + j ] = sin( Turb->Kmode[1][n]*j*dh );
         Turb->Cos[1][ n*NPoint + j ] = cos( Turb->Kmode[1][n]*j*dh );
      }
      for (int k = 0; k < TURB_TABLE_SIZE; k++)
      {
         Turb->Sin[2][ n*NPoint + k ] = sin( Turb->Kmode[2][n]*k*dh );
         Turb->Cos[2][ n*NPoint + k ] = cos( Turb->Kmode[2][n]*k*dh );
      }
//    apply periodicity
      for (int d = 0; d < 3; d++)
      {
         Turb->Sin[d][ n*NPoint + TURB_TABLE_SIZE ] = Turb->Sin[d][ n*NPoint ];
         Turb->Cos[d][ n*NPoint + TURB_TABLE_SIZE ] = Turb->Cos[d][ n*NPoint ];
      }
   }

// set next update time
   Turb->TimeLast = Time[0];
   Turb->TimeNext = Time[0] + Turb->dt;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Turb_Init


//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_End
// Description :  Free memories
//-------------------------------------------------------------------------------------------------------
void Turb_End()
{
   if ( Turb != NULL ) delete Turb;

} // FUNCTION : Turb_End


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

} // FUNCTION : Turb_FillinTable


//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_GetRNG
// Description :  Get random number using Box–Muller transformation based on RSeed, and multiply by OUvar
//
// Return      :  Gaussian random number pair
//-------------------------------------------------------------------------------------------------------
void Turb_GetRNG( double& a, double& b, int& Seed, const double OUvar )
{
   double r1 = Turb_ran1s(Seed);
   double r2 = Turb_ran1s(Seed);
   a = OUvar*sqrt(2.0*log(1.0/r1))*cos(2*M_PI*r2);
   b = OUvar*sqrt(2.0*log(1.0/r1))*sin(2*M_PI*r2);
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Turb_ran1s
// Description :  Park–Miller random number generator, and update random seed
//
// Return      :  uniformly distributed random number in [0,1[
//-------------------------------------------------------------------------------------------------------
double Turb_ran1s(int& Seed)
{
   static const int IA=16807, IM=2147483647, IQ=127773, IR=2836;
   static const double AM=1.0/IM, RNMX=1.0-1.2e-7;
   if (Seed <= 0) Seed = MAX(-Seed, 1);
   int k = Seed/IQ;
   Seed = IA*(Seed-k*IQ)-IR*k;
   if (Seed < 0) Seed = Seed+IM;
   int iy = Seed;

   return MIN(AM*iy, RNMX);
}


#endif // if ( MODEL == HYDRO )
