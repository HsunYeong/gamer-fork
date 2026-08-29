#include "GAMER.h"



// =======================================================================================
// CR_Streaming_Test : initial-condition geometry selector (Jiang & Oh 2018, Section 4)
// --> the CR *physics* parameters (CR_SIGMA, CR_SIGMA_PERP, CR_STREAM, CR_SOURCE, CR_VMAX, ...)
//     are taken from Input__Parameter, NOT set here.  Each example directory carries the
//     Input__Parameter/Input__TestProb appropriate for one paper test.
// =======================================================================================
#define CR_TEST_TRIANGULAR_1D     0     // Sec 4.1.1 : 1D triangular Ec profile      (streaming)
#define CR_TEST_GAUSSIAN_1D       1     // Sec 4.1.1 / 4.1.4 : 1D Gaussian Ec         (streaming or diffusion)
#define CR_TEST_GAUSSIAN_2D       2     // Sec 4.1.2 / 4.1.5 : 2D Gaussian Ec, diagonal B (streaming or anisotropic diffusion)
#define CR_TEST_CIRCLE_2D         3     // Sec 4.1.5 : 2D ring Ec, circular B         (anisotropic diffusion)
#define CR_TEST_BOTTLENECK_1D     4     // Sec 4.1.3 : 1D density cloud, CR injected from the -x boundary
#define CR_TEST_WAVE_1D           5     // Sec 4.2.1 : 1D CR-driven wave              (full MHD + CR coupling)
#define CR_TEST_BLAST_2D          6     // Sec 4.2.3 : 2D CR-driven blast wave        (full MHD + CR coupling)
// -- classic single-moment CR module comparison tests (port of example/test_problem/Hydro/CR_{Diffusion,ShockTube,SoundWave}) --
// -- compare the two-moment CR_E against the classic-module CRay for the same setup --
#define CR_TEST_CLASSIC_DIFFUSION 7     // classic CR_Diffusion : Gaussian Ec, uniform B||x, frozen gas (pure diffusion)
#define CR_TEST_CLASSIC_SHOCKTUBE 8     // classic CR_ShockTube : CR-hydro shock, CR locked to gas (high opacity + CR_SOURCE)
#define CR_TEST_CLASSIC_SOUNDWAVE 9     // classic CR_SoundWave : CR-modified acoustic wave, CR locked to gas
// -- the paper's own CR-modified shock test (Jiang & Oh 2018, Sec 4.2.2, Fig 11) --
#define CR_TEST_PAPER_SHOCK      10     // Sec 4.2.2 : colliding-flow CR-modified shock (streaming + diffusion + advection)
// -- 3D generalization of the CR-driven blast (Sec 4.2.3) with a configurable B direction --
#define CR_TEST_BLAST_3D         11     // 3D CR-driven blast wave; B direction set by CR_Streaming_B_theta/B_phi


// problem-specific global variables
// =======================================================================================
static int    CR_Streaming_Test;        // test selector (one of the CR_TEST_* macros above)
static int    CR_Streaming_Dir;          // streaming/diffusion axis for the 1D tests (0/1/2 = x/y/z)
static double CR_Streaming_FlowV;        // uniform background flow velocity along CR_Streaming_Dir
                                         // (e.g. the moving-fluid diffusion test, Sec 4.1.4);
                                         // also the collision half-speed for the shock test (Sec 4.2.2, |v|=10)
static double CR_Streaming_Ec0;          // uniform CR energy density for the shock test (Sec 4.2.2: 1, 50, or 200)
static int    CR_Streaming_FcInit;       // shock-test initial CR flux: 1 = advective Fc=(4/3)v*Ec (paper),
                                         // 0 = zero (avoids the sharp t=0 flux discontinuity at the collision;
                                         //     the flux relaxes to the advective value on ~1/(Vm*sigma))
static double CR_Streaming_B_theta;      // polar angle (deg, from +z axis) of the uniform B field for the 3D blast
static double CR_Streaming_B_phi;        // azimuthal angle (deg, in the x-y plane from +x axis) of the uniform B
                                         // field for the 3D blast --> B = (sinT cosP, sinT sinP, cosT), |B| = 1;
                                         // the default (theta=90, phi=0) reproduces the 2D-blast B||x
                                         // (only used by CR_TEST_BLAST_3D)
static int    CR_Streaming_GradOutflowBC; // 1 = gradient-preserving CR outflow BC on both x faces (streaming
                                         // tests): linear-extrapolate Ec (positive-clamped) + copy Fc so the
                                         // boundary Pc gradient keeps full strength and the profile tracks the
                                         // Jiang & Oh (2018) analytic solution to the domain edge (Fig. 3);
                                         // 0 = default public-Athena outflow copy (Ec ghost = last active cell)
// =======================================================================================


// function prototypes shared within this file
#if ( MODEL == HYDRO  &&  defined CR_STREAMING )
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );
void ShockBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
              const int GhostSize, const int idx[], const double pos[], const double Time,
              const int lv, const int TFluVarIdxList[], double AuxArray[] );
void GradOutflowBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                    const int GhostSize, const int idx[], const double pos[], const double Time,
                    const int lv, const int TFluVarIdxList[], double AuxArray[] );
bool Flag_CR_Streaming( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef CR_STREAMING
   Aux_Error( ERROR_INFO, "CR_STREAMING must be enabled !!\n" );
#  endif

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled (CR streaming/diffusion is field-aligned) !!\n" );
#  endif

// the two-moment CR_STREAMING module is standalone: the gas may use a pure gamma-law EoS
// (recommended, fully decoupled from CRAY) or the COSMIC_RAY EoS during the transition period
#  if ( EOS != EOS_GAMMA  &&  EOS != EOS_COSMIC_RAY )
   Aux_Error( ERROR_INFO, "EOS must be EOS_GAMMA (standalone) or EOS_COSMIC_RAY !!\n" );
#  endif

#  ifdef CR_DIFFUSION
   Aux_Error( ERROR_INFO, "CR_DIFFUSION must be disabled !!\n" );
#  endif


// warnings


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == HYDRO  &&  defined CR_STREAMING )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//                3. Must call EoS_Init() before calling any other EoS routine
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "CR_Streaming_Test",  &CR_Streaming_Test,      0,           0,                11                );
   ReadPara->Add( "CR_Streaming_Dir",   &CR_Streaming_Dir,       0,           0,                3                 );
   ReadPara->Add( "CR_Streaming_FlowV", &CR_Streaming_FlowV,     0.0,         NoMin_double,     NoMax_double      );
   ReadPara->Add( "CR_Streaming_Ec0",   &CR_Streaming_Ec0,       1.0,         Eps_double,       NoMax_double      );
   ReadPara->Add( "CR_Streaming_FcInit", &CR_Streaming_FcInit,    1,           0,                1                 );
   ReadPara->Add( "CR_Streaming_GradOutflowBC", &CR_Streaming_GradOutflowBC, 0, 0,               1                 );
   ReadPara->Add( "CR_Streaming_B_theta", &CR_Streaming_B_theta,  90.0,        NoMin_double,     NoMax_double      );
   ReadPara->Add( "CR_Streaming_B_phi",   &CR_Streaming_B_phi,     0.0,        NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
// the configurable B direction is only wired up for the 3D blast; guard against silently
// ignoring a non-default angle for any other test
   if ( CR_Streaming_Test != CR_TEST_BLAST_3D  &&
        ( CR_Streaming_B_theta != 90.0  ||  CR_Streaming_B_phi != 0.0 ) )
      Aux_Error( ERROR_INFO, "CR_Streaming_B_theta/B_phi are only supported for CR_Streaming_Test = %d (3D blast) !!\n",
                 CR_TEST_BLAST_3D );


// (2) set the problem-specific derived parameters
   const bool is_1D = ( CR_Streaming_Test == CR_TEST_TRIANGULAR_1D     ||
                        CR_Streaming_Test == CR_TEST_GAUSSIAN_1D       ||
                        CR_Streaming_Test == CR_TEST_BOTTLENECK_1D     ||
                        CR_Streaming_Test == CR_TEST_WAVE_1D           ||
                        CR_Streaming_Test == CR_TEST_CLASSIC_DIFFUSION ||
                        CR_Streaming_Test == CR_TEST_CLASSIC_SHOCKTUBE ||
                        CR_Streaming_Test == CR_TEST_PAPER_SHOCK         );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// for the 1D tests, output a line along the streaming direction
   if ( is_1D )
   {
      const int target = ( CR_Streaming_Dir == 0 ) ? OUTPUT_X : ( CR_Streaming_Dir == 1 ) ? OUTPUT_Y : OUTPUT_Z;
      if ( OPT__OUTPUT_PART != target )
      {
         OPT__OUTPUT_PART = target;
         PRINT_RESET_PARA( OPT__OUTPUT_PART, FORMAT_INT, "" );
      }
   }

// the gradient-preserving CR outflow BC needs the -x/+x fluid faces set to the user BC (=4);
// reset them here so enabling the flag is a single switch (the user BC is enrolled in Init_TestProb)
   if ( CR_Streaming_GradOutflowBC )
   {
      if ( OPT__BC_FLU[0] != BC_FLU_USER )
      {
         OPT__BC_FLU[0] = BC_FLU_USER;
         PRINT_RESET_PARA( OPT__BC_FLU[0], FORMAT_INT, "(gradient-preserving CR outflow BC)" );
      }
      if ( OPT__BC_FLU[1] != BC_FLU_USER )
      {
         OPT__BC_FLU[1] = BC_FLU_USER;
         PRINT_RESET_PARA( OPT__BC_FLU[1], FORMAT_INT, "(gradient-preserving CR outflow BC)" );
      }
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  CR_Streaming_Test     = %d\n",     CR_Streaming_Test  );
      Aux_Message( stdout, "  CR_Streaming_Dir      = %d\n",     CR_Streaming_Dir   );
      Aux_Message( stdout, "  CR_Streaming_FlowV    = %14.7e\n", CR_Streaming_FlowV );
      Aux_Message( stdout, "  CR_Streaming_Ec0      = %14.7e\n", CR_Streaming_Ec0   );
      Aux_Message( stdout, "  CR_Streaming_FcInit   = %d\n",     CR_Streaming_FcInit );
      Aux_Message( stdout, "  CR_Streaming_GradOutflowBC = %d\n", CR_Streaming_GradOutflowBC );
      Aux_Message( stdout, "  CR_Streaming_B_theta  = %14.7e\n", CR_Streaming_B_theta );
      Aux_Message( stdout, "  CR_Streaming_B_phi    = %14.7e\n", CR_Streaming_B_phi   );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  CR_Bottleneck_Dens
// Description :  Density profile of the cold cloud used in the bottleneck test (Sec 4.1.3, Eq. 24)
//
// Note        :  1. Uses ABSOLUTE coordinates (the cloud is centered at x0 = 200 in the domain x in [0,1000])
//                2. rho = rho_h + (rho_c - rho_h)*[1+tanh((x-x0)/dx0)]*[1+tanh((x0-x)/dx0)]
//                   --> a dense cold cloud (rho_c) embedded in a diffuse hot background (rho_h)
//
// Parameter   :  x : Position along the streaming direction (absolute coordinate)
//
// Return      :  gas density
//-------------------------------------------------------------------------------------------------------
static double CR_Bottleneck_Dens( const double x )
{
   const double rho_c = 1.0, rho_h = 0.1, x0 = 200.0, dx0 = 25.0;
   const double t1 = 1.0 + std::tanh( (x  - x0)/dx0 );
   const double t2 = 1.0 + std::tanh( (x0 - x )/dx0 );
   return rho_h + (rho_c - rho_h)*t1*t2;
} // FUNCTION : CR_Bottleneck_Dens



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//                5. The ADV_* fields (ADV_SIGMA/ADV_VX/ADV_VY/ADV_VZ) are recomputed by CR_UpdateOpacity()
//                   at the start of every fluid solve, so the IC values set here are overwritten before use
//                   --> we therefore set safe no-streaming defaults
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double xc = amr->BoxCenter[0];
   const double yc = amr->BoxCenter[1];

   double Dens = 1.0, vx = 0.0, vy = 0.0, vz = 0.0, Pgas = 1.0, cr_E = 0.0;
   double cr_F1 = 0.0;   // initial CR flux along x (zero for every test except the Sec 4.2.2 shock)

   switch ( CR_Streaming_Test )
   {
      case CR_TEST_TRIANGULAR_1D :
      {
//       Sec 4.1.1: Ec = 2 - |r| for |r| < 1, else 1   (r relative to the box center)
         const double r = ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z ) - amr->BoxCenter[CR_Streaming_Dir];
         const double d = std::fabs( r );
         cr_E = ( d < 1.0 ) ? (2.0 - d) : 1.0;
         vx   = ( CR_Streaming_Dir==0 ) ? CR_Streaming_FlowV : 0.0;
         vy   = ( CR_Streaming_Dir==1 ) ? CR_Streaming_FlowV : 0.0;
         vz   = ( CR_Streaming_Dir==2 ) ? CR_Streaming_FlowV : 0.0;
         break;
      }

      case CR_TEST_GAUSSIAN_1D :
      {
//       Sec 4.1.1 / 4.1.4: Ec = exp(-40 r^2)   (r relative to the box center)
         const double r = ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z ) - amr->BoxCenter[CR_Streaming_Dir];
         cr_E = std::exp( -40.0 * r*r );
         vx   = ( CR_Streaming_Dir==0 ) ? CR_Streaming_FlowV : 0.0;
         vy   = ( CR_Streaming_Dir==1 ) ? CR_Streaming_FlowV : 0.0;
         vz   = ( CR_Streaming_Dir==2 ) ? CR_Streaming_FlowV : 0.0;
         break;
      }

      case CR_TEST_GAUSSIAN_2D :
      {
//       Sec 4.1.2 / 4.1.5: Ec = exp(-40 (dx^2 + dy^2)),  B along the diagonal (set in SetBFieldIC)
         const double dx = x - xc, dy = y - yc;
         cr_E = std::exp( -40.0 * (dx*dx + dy*dy) );
         break;
      }

      case CR_TEST_CIRCLE_2D :
      {
//       Sec 4.1.5: a ring of Ec; B is circular (set in SetBFieldIC)
//       Ec = 12 in 0.5 < r < 0.7 and |phi| < pi/12 (consistent with the analytic solution Eq. 28), else 10
//       The paper keeps "all the MHD variables fixed": the circular field B = (-y/r, x/r) has |B| = 1 but an
//       unbalanced hoop stress (tension force -B^2/r r^) that cannot be pressure-balanced (Pgas ~ -ln r diverges).
//       With rho = 1 the gas implodes at the field vortex (rho_center 1 -> 2.7 by t = 0.26) and adiabatically
//       compresses the CRs into a spurious central peak (Ec ~ 35 vs background 10).  Emulate the frozen fluid
//       with a very inert gas instead: a = B^2/(r*rho) --> rho = 1e8 keeps the center displacement < 1e-5 cells
//       by t = 0.26.  With streaming off the CR module never uses rho, so the CR test itself is unaffected.
         const double dx = x - xc, dy = y - yc;
         const double r   = std::sqrt( dx*dx + dy*dy );
         const double phi = std::atan2( dy, dx );
         cr_E = ( r > 0.5  &&  r < 0.7  &&  std::fabs(phi) < M_PI/12.0 ) ? 12.0 : 10.0;
         Dens = 1.0e8;                           // quasi-static gas (see the note above)
         break;
      }

      case CR_TEST_BOTTLENECK_1D :
      {
//       Sec 4.1.3: a cold dense cloud; CRs are injected from the -x boundary (see BottleneckBC)
//       --> the simulation domain starts essentially empty of CRs
         Dens = CR_Bottleneck_Dens( x );        // absolute coordinate (cloud at x0 = 200)
         cr_E = 1.0e-6;
         break;
      }

      case CR_TEST_WAVE_1D :
      {
//       Sec 4.2.1: Ec = 20 + 10 sin(pi*(x-xc)); uniform gas; full MHD + CR coupling (CR_SOURCE=1)
         cr_E = 20.0 + 10.0*std::sin( M_PI*(x - xc) );
         break;
      }

      case CR_TEST_BLAST_2D :
      {
//       Sec 4.2.3: Ec = 100 inside r < 0.02, else 0.1; uniform background (rho=1, Eint=2.5)
         const double dx = x - xc, dy = y - yc;
         const double r  = std::sqrt( dx*dx + dy*dy );
         cr_E = ( r < 0.02 ) ? 100.0 : 0.1;
         Pgas = (GAMMA - 1.0)*2.5;               // background internal energy = 2.5
         break;
      }

      case CR_TEST_BLAST_3D :
      {
//       3D generalization of Sec 4.2.3: identical to CR_TEST_BLAST_2D but with a spherical
//       (rather than cylindrical) CR overpressure and a configurable B direction (set in
//       SetBFieldIC via CR_Streaming_B_theta/B_phi):
//       Ec = 100 inside r < 0.02, else 0.1; uniform background (rho=1, Eint=2.5)
         const double dx = x - xc, dy = y - yc, dz = z - amr->BoxCenter[2];
         const double r  = std::sqrt( dx*dx + dy*dy + dz*dz );
         cr_E = ( r < 0.02 ) ? 100.0 : 0.1;
         Pgas = (GAMMA - 1.0)*2.5;               // background internal energy = 2.5
         break;
      }

      case CR_TEST_CLASSIC_DIFFUSION :
      {
//       Port of the classic CR_Diffusion test (Gaussian ball, Type 0, Mag_Type 0) with the EXACT
//       classic main-branch setup for a direct CR_E-vs-CRay comparison:
//         Ec(t=0) = E0*exp(-R02*r^2) + BG  (3D ball at the box center),  uniform B||x,
//         rho = 1, Pgas = 5/3 (classic CR_Diffusion_PGas0 default).
//       The classic parallel coefficient kappa_para is matched in Input__Parameter via the two-moment
//       inverse-diffusion convention kappa = 1/(3*CR_SIGMA)  -->  CR_SIGMA = 1/(3*kappa_classic);
//       kappa_perp = 0 is emulated with a large-but-stable CR_SIGMA_PERP (see the example input).
//       CR_SOURCE selects the physics: 1 = CR back-reacts on the gas exactly like the classic
//       module's EoS coupling (matched live-gas comparison); 0 = no coupling --> the gas stays
//       static (uniform Pgas) and the center line follows the pure-diffusion analytic (Eq. 25).
         const double E0 = 1.0, R02 = 40.0, BG = 0.1;
         const double dxc = x - xc, dyc = y - yc, dzc = z - amr->BoxCenter[2];
         cr_E = E0*std::exp( -R02*( dxc*dxc + dyc*dyc + dzc*dzc ) ) + BG;
         Pgas = 1.666666666667;                  // classic CR_Diffusion_PGas0
         break;
      }

      case CR_TEST_CLASSIC_SHOCKTUBE :
      {
//       Port of the classic CR_ShockTube test (CR coupled to gas, gamma_cr=4/3, no diffusion).
//       Emulated as single-fluid CR-hydro: CR locked to the gas (CR_SIGMA large, CR_STREAM=0) with
//       CR_SOURCE=1 providing the CR-pressure back-reaction.  Ec = Pcr/(gamma_cr-1) = 3*Pcr.
//       L (x<xc): rho=1.0, Pgas=6.7e4, Pcr=1.3e5 ;  R (x>xc): rho=0.2, Pgas=2.4e2, Pcr=2.4e2.
         const double Pcr_L = 1.3e5, Pcr_R = 2.4e2;
         const double _gcrm1 = 1.0/(GAMMA_CR - 1.0);
         if      ( x < xc ) { Dens = 1.0;            Pgas = 6.7e4;             cr_E = Pcr_L*_gcrm1;             }
         else if ( x > xc ) { Dens = 0.2;            Pgas = 2.4e2;             cr_E = Pcr_R*_gcrm1;             }
         else               { Dens = 0.5*(1.0+0.2); Pgas = 0.5*(6.7e4+2.4e2); cr_E = 0.5*(Pcr_L+Pcr_R)*_gcrm1; }
         break;
      }

      case CR_TEST_CLASSIC_SOUNDWAVE :
      {
//       Port of the classic CR_SoundWave test (CR-modified acoustic wave, gamma_cr=4/3).
//       Emulated as single-fluid CR-hydro: CR locked to the gas (CR_SIGMA large, CR_STREAM=0) with
//       CR_SOURCE=1.  The wave speed includes the CR pressure: cs^2 = (gamma*Pg0 + gamma_cr*Pcr0)/rho0.
//       CR_Streaming_Dir = 0/1/2 -> grid-aligned 1D wave along x/y/z (B along same axis);
//       CR_Streaming_Dir = 3     -> diagonal wave along (1,1,1) with B = (1,1,1)/sqrt(3).
         const double Rho0 = 1.0, Pg0 = 1.0, Pcr0 = 1.0, Delta = 1.0e-6, Sign = 1.0;
         const double cs       = std::sqrt( GAMMA*Pg0 + GAMMA_CR*Pcr0 );   // rho0 = 1
         const double delta_cs = Delta/cs;
         const bool    diag    = ( CR_Streaming_Dir == 3 );
         const double WaveL    = diag ? amr->BoxSize[0]/std::sqrt(3.0) : amr->BoxSize[CR_Streaming_Dir];
         const double WaveK    = 2.0*M_PI/WaveL;
         const double r        = diag ? ( x + y + z )/std::sqrt(3.0)
                                      : ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z );
         const double wave     = std::sin( WaveK*r );                      // t = 0
         const double vel      = Sign*Delta*wave;                          // velocity perturbation
         Dens = ( 1.0 + delta_cs*wave )*Rho0;
         if ( diag ) {
            vx = vy = vz = vel/std::sqrt(3.0);
         } else {
            vx = ( CR_Streaming_Dir==0 ) ? vel : 0.0;
            vy = ( CR_Streaming_Dir==1 ) ? vel : 0.0;
            vz = ( CR_Streaming_Dir==2 ) ? vel : 0.0;
         }
         Pgas = ( 1.0 + delta_cs*wave*GAMMA    )*Pg0;
         cr_E = ( 1.0 + delta_cs*wave*GAMMA_CR )*Pcr0/(GAMMA_CR - 1.0);    // Ec = 3*Pcr
         break;
      }

      case CR_TEST_PAPER_SHOCK :
      {
//       Sec 4.2.2 (Fig 11): two symmetric gas streams collide at the box center and drive a
//       CR-modified shock outward to each side.  This is the paper's OWN two-moment CR shock
//       (full streaming + diffusion + advection), NOT the classic-module port (mode 8).
//         rho = 1, Pgas = 1 uniform;   v = +FlowV for x<xc, -FlowV for x>xc   (paper: FlowV=10)
//         Ec  = CR_Streaming_Ec0 uniform (paper: 1, 50, 200 -> Pc = Ec/3 = 1/3, 50/3, 200/3)
//         initial CR flux is the advective flux Fc = (4/3) v Ec (CRs co-moving with the gas)
//       Streaming uses a constant v_A = 1 (uniform B||x with rho=1 upstream); diffusion sigma'=10
//       and Vm=100 come from Input__Parameter.  The x ghost zones are fixed to this IC (ShockBC).
//       NOTE: the stored CR_F field is the REDUCED flux Fc/Vm (same convention as Athena++'s
//       u_cr(CRF*); see the transport flux vmax*Fc and the source equilibrium Fc -> (4/3)v*Ec/Vm
//       in CPU_CR_TwoMoment.cpp), so the advective IC must be divided by CR_VMAX here.
         const double v = ( x < xc ) ? CR_Streaming_FlowV : -CR_Streaming_FlowV;
         Dens  = 1.0;
         Pgas  = 1.0;
         vx    = v;
         cr_E  = CR_Streaming_Ec0;
         cr_F1 = ( CR_Streaming_FcInit == 0 ) ? 0.0 : (4.0/3.0)*v*cr_E/CR_VMAX;
         break;
      }

      default :
         Aux_Error( ERROR_INFO, "unsupported CR_Streaming_Test (%d) !!\n", CR_Streaming_Test );
   } // switch ( CR_Streaming_Test )


// momentum
   const double MomX = Dens*vx;
   const double MomY = Dens*vy;
   const double MomZ = Dens*vz;

// gas pressure for the EoS
// --> standalone two-moment build (EOS_GAMMA, no COSMIC_RAY): the gas is a pure gamma-law fluid;
//     ALL CR back-reaction is handled by the two-moment source terms (CR_E/CR_F), so no CR pressure
//     is folded into the gas energy here
   double Pres = Pgas;
#  ifdef COSMIC_RAY
// COSMIC_RAY EoS path (transition only): the EoS needs a CRAY value, so add a constant background
// CR pressure. This couples the gas EoS to CRAY and is exactly what the standalone build removes.
   const double P_cr = 1.0;
   Pres       += P_cr;
   fluid[CRAY] = P_cr / (GAMMA_CR - 1.0);
#  endif

   const double Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt,
                                                 EoS_AuxArray_Int, h_EoS_Table );
   const double Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );   // do NOT include magnetic energy here

// two-moment CR fields (Ec and the three flux components; flux = 0 except the Sec 4.2.2 shock,
// which initializes the advective flux Fc = (4/3) v Ec)
   fluid[CR_E ] = cr_E;
   fluid[CR_F1] = cr_F1;
   fluid[CR_F2] = 0.0;
   fluid[CR_F3] = 0.0;

// streaming opacity/velocity: recomputed by CR_UpdateOpacity() each step --> safe no-streaming defaults
   fluid[ADV_SIGMA] = CR_MAX_OPACITY;
   fluid[ADV_VX   ] = 0.0;
   fluid[ADV_VY   ] = 0.0;
   fluid[ADV_VZ   ] = 0.0;

// gas fields
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                2. The Alfven speed used by the streaming module is v_A = |B|/sqrt(rho); the field
//                   geometry also sets the direction of field-aligned streaming/diffusion
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   switch ( CR_Streaming_Test )
   {
      case CR_TEST_TRIANGULAR_1D :
      case CR_TEST_GAUSSIAN_1D :
//       uniform field along the streaming direction, |B| = 1  --> v_A = 1
         magnetic[MAGX] = ( CR_Streaming_Dir==0 ) ? 1.0 : 0.0;
         magnetic[MAGY] = ( CR_Streaming_Dir==1 ) ? 1.0 : 0.0;
         magnetic[MAGZ] = ( CR_Streaming_Dir==2 ) ? 1.0 : 0.0;
         break;

      case CR_TEST_GAUSSIAN_2D :
//       uniform field along the x-y diagonal, |B| = 1  --> v_A = 1   (Sec 4.1.2 / 4.1.5)
         magnetic[MAGX] = 1.0/std::sqrt(2.0);
         magnetic[MAGY] = 1.0/std::sqrt(2.0);
         magnetic[MAGZ] = 0.0;
         break;

      case CR_TEST_CIRCLE_2D :
      {
//       circular field B = ( -(y-yc), (x-xc) ) / r, |B| = 1   (Sec 4.1.5)
//       --> set B as the DISCRETE CURL of the vector potential Az(X,Y) = -sqrt((X-xc)^2+(Y-yc)^2):
//              Bx = dAz/dy ,  By = -dAz/dx
//           evaluated by one-cell finite differences over the cell faces.
//       --> this makes B divergence-free to MACHINE PRECISION on the GAMER grid and reproduces
//           Athena++'s setup (src/pgen/cr_diffusion.cpp) cell-for-cell.
//       --> the previous pointwise field B=(-dy/r,dx/r) is only *analytically* div-free; on the
//           grid it leaves a large discrete div(B) (~1e-2 max, every cell failing the 1e-11
//           tolerance), so the field-aligned diffusion direction b(x)b(x) differs from Athena's
//           and the GAMER<->Athena error grows with the diffusion time (the uniform-field
//           anisotropic-diffusion test, which IS discretely div-free, stays at ~3e-8 instead).
         const double xc = amr->BoxCenter[0];
         const double yc = amr->BoxCenter[1];
         const double dh = amr->dh[lv];                 // SetBFieldIC is called per component at the
                                                        // correct face center, so a single-cell curl
                                                        // here matches GAMER's face-centered B storage
#        define CR_CIRCLE_AZ( X, Y )  ( -std::sqrt( SQR((X)-xc) + SQR((Y)-yc) ) )
         magnetic[MAGX] =  ( CR_CIRCLE_AZ( x, y+0.5*dh ) - CR_CIRCLE_AZ( x, y-0.5*dh ) ) / dh;
         magnetic[MAGY] = -( CR_CIRCLE_AZ( x+0.5*dh, y ) - CR_CIRCLE_AZ( x-0.5*dh, y ) ) / dh;
         magnetic[MAGZ] = 0.0;
#        undef CR_CIRCLE_AZ
         break;
      }

      case CR_TEST_CLASSIC_SOUNDWAVE :
//       B aligned with the wave so the field-aligned CR source has no perpendicular cross-terms:
//       CR_Streaming_Dir = 0/1/2 -> B along x/y/z (grid-aligned); = 3 -> B = (1,1,1)/sqrt(3)
         if ( CR_Streaming_Dir == 3 ) {
            magnetic[MAGX] = 1.0/std::sqrt(3.0);
            magnetic[MAGY] = 1.0/std::sqrt(3.0);
            magnetic[MAGZ] = 1.0/std::sqrt(3.0);
         } else {
            magnetic[MAGX] = ( CR_Streaming_Dir==0 ) ? 1.0 : 0.0;
            magnetic[MAGY] = ( CR_Streaming_Dir==1 ) ? 1.0 : 0.0;
            magnetic[MAGZ] = ( CR_Streaming_Dir==2 ) ? 1.0 : 0.0;
         }
         break;

      case CR_TEST_BLAST_3D :
      {
//       uniform field in the direction (theta, phi), |B| = 1  --> v_A = 1  (rho = 1):
//          B = ( sinT cosP, sinT sinP, cosT ),  theta = polar angle from +z, phi = azimuth from +x
//       the angles are given in degrees; the default (theta=90, phi=0) recovers B||x (the 2D blast)
         const double theta = CR_Streaming_B_theta * M_PI/180.0;
         const double phi   = CR_Streaming_B_phi   * M_PI/180.0;
         magnetic[MAGX] = std::sin(theta)*std::cos(phi);
         magnetic[MAGY] = std::sin(theta)*std::sin(phi);
         magnetic[MAGZ] = std::cos(theta);
         break;
      }

      case CR_TEST_BOTTLENECK_1D :
      case CR_TEST_WAVE_1D :
      case CR_TEST_BLAST_2D :
      case CR_TEST_CLASSIC_DIFFUSION :   // B||x: classic diffusion is along x (kappa_para)
      case CR_TEST_CLASSIC_SHOCKTUBE :   // B||x (parallel to the 1D shock normal -> no magnetic force)
      case CR_TEST_PAPER_SHOCK :         // B||x, |B|=1 -> v_A = 1 upstream (rho=1), as in Sec 4.2.2
      default :
//       uniform field along x, |B| = 1
         magnetic[MAGX] = 1.0;
         magnetic[MAGY] = 0.0;
         magnetic[MAGZ] = 0.0;
         break;
   } // switch ( CR_Streaming_Test )

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  BottleneckBC
// Description :  User boundary condition for the bottleneck test (Sec 4.1.3): inject CRs from the
//                -x boundary by fixing Ec = 3 while keeping the background gas profile
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//                2. Only the -x face is set to the user BC (OPT__BC_FLU_XM = 4); the +x face uses outflow
//                3. Following Jiang & Oh (2018, Sec 4.1.3), the boundary CR flux is "reflecting": the
//                   ghost-zone Fc is set to the sign-flipped value of the first active zone
//                   (Fc_ghost = -Fc_active), while Ec is fixed to 3.  This matches the Athena++ setup
//                   (src/pgen/cr_diffusion.cpp : BottleneckCRInnerX1) so the boundary injects exactly
//                   the same CR energy flux.
//                   --> a simpler Fc_ghost = 0 also drives the bottleneck, but the HLLE boundary flux
//                       then leaves the upstream Ec plateau only ~half as far below the Ec=3 reservoir,
//                       i.e. ~2% high relative to Athena++ (the only place the two codes disagreed)
//                4. The first active zone's Fc is read from the prepared "Array" (its interior is
//                   already filled before the ghost zones; same layout/pattern as the JetICMWall user BC)
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                ArraySize      : Size of Array including the ghost zones on each side
//                fluid          : Fluid fields to be set
//                NVar_Flu       : Number of fluid variables to be prepared
//                GhostSize      : Number of ghost zones
//                idx            : Array indices
//                pos            : Physical coordinates
//                Time           : Physical time
//                lv             : Refinement level
//                TFluVarIdxList : List recording the target fluid variable indices
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BottleneckBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                   const int GhostSize, const int idx[], const double pos[], const double Time,
                   const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// start from the background IC (gas density profile, Ec = 1e-6, Fc = 0, ADV defaults)
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

// fix the injected CR energy density; keep the transverse CR fluxes at zero (B is along x here)
   fluid[CR_E ] = 3.0;
   fluid[CR_F2] = 0.0;
   fluid[CR_F3] = 0.0;

// reflecting CR flux along the streaming (x) direction: Fc_ghost = -Fc(first active zone)
// --> the prepared "Array" already holds the interior data when the ghost zones are filled, so we
//     read CR_F1 of the first active cell; slot v in Array corresponds to field TFluVarIdxList[v]
   typedef real (*vla)[ ArraySize[2] ][ ArraySize[1] ][ ArraySize[0] ];
   vla Array3D = ( vla )Array;

   const int i_ref = GhostSize;   // first active cell along +x (the -x ghost zones are idx[0] < GhostSize)
   const int jg    = idx[1];
   const int kg    = idx[2];

   real CRF1_active = 0.0;
   for (int v=0; v<NVar_Flu; v++)
      if ( TFluVarIdxList[v] == CR_F1 )   CRF1_active = Array3D[v][kg][jg][i_ref];

   fluid[CR_F1] = -CRF1_active;

} // FUNCTION : BottleneckBC



//-------------------------------------------------------------------------------------------------------
// Function    :  ShockBC
// Description :  User boundary condition for the Sec 4.2.2 CR-modified shock test: the x ghost zones are
//                fixed to the initial condition ("variables in the left and right ghost zones are fixed
//                to be these initial values", Jiang & Oh 2018)
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr" for both the -x and +x faces
//                2. The shock forms at the box center and does not reach the boundary within the run
//                   time, so this simply keeps the inflowing streams (v = +/-FlowV) steady at the edges
//
// Parameter   :  (same as BottleneckBC; only "fluid" and "pos" are used here)
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void ShockBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
              const int GhostSize, const int idx[], const double pos[], const double Time,
              const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// fixed ghost zones = initial condition
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

} // FUNCTION : ShockBC



//-------------------------------------------------------------------------------------------------------
// Function    :  GradOutflowBC
// Description :  Gradient-preserving CR outflow boundary condition on the -x and +x faces (streaming tests)
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr" for BOTH the -x and +x faces
//                   (OPT__BC_FLU_XM/XP are reset to user=4 in SetParameter when CR_Streaming_GradOutflowBC=1)
//                2. The gas (DENS/MOM/ENGY) keeps standard outflow behavior; for the streaming tests the gas
//                   is frozen (CR_SOURCE=0, uniform IC), so the SetGridIC background equals the outflow-copy
//                   state exactly.  The CR fluxes Fc1/2/3 are copied from the last active cell (outflow).
//                3. The CR energy density is LINEARLY EXTRAPOLATED along the face normal from the two nearest
//                   active cells (positive-clamped): Ec_ghost(n) = max( Ec_last + n*(Ec_last - Ec_lastm1), TINY ).
//                   --> the enrolled opacity central-difference b.grad(Pc) and the transport flux both see the
//                       interior slope through the boundary, so v_adv and sigma_adv keep full strength and the
//                       CRs stream out at the physical rate (4/3) v_A Ec instead of jamming (Jiang & Oh 2018 Fig 3).
//                   --> the public-Athena copy BC (default) instead halves b.grad(Pc) at the last active cell
//                       (Ec_ghost = Ec_last), which flattens Ec over the last ~12 cells.
//                4. The ghost streaming fields (ADV_SIGMA/ADV_VX/VY/VZ) do NOT need to be filled here: they are
//                   recomputed from the extrapolated Ec by CR_UpdateOpacity() over indices [1, FLU_NXT-2] before
//                   the flux uses them (see CPU_FluidSolver_MHM.cpp).  We leave the SetGridIC defaults.
//                5. The interior (active) cells are already filled in the prepared "Array" when the domain-boundary
//                   ghost zones are set, so we read the last two active cells directly (same mechanism as BottleneckBC).
//
// Parameter   :  (same as BottleneckBC / ShockBC)
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void GradOutflowBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                    const int GhostSize, const int idx[], const double pos[], const double Time,
                    const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// start from the frozen background (correct gas state for CR_SOURCE=0 + safe CR/ADV defaults)
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

// map the prepared 1D "Array" to 4D [NVar_Flu][k][j][i]
   typedef real (*vla)[ ArraySize[2] ][ ArraySize[1] ][ ArraySize[0] ];
   vla Array3D = ( vla )Array;

   const int ig = idx[0];
   const int jg = idx[1];
   const int kg = idx[2];
   const int Nx = ArraySize[0];

// identify the x-face and the two nearest active cells along the outward normal
   int i_last, i_lastm1, ndepth;
   if ( ig < GhostSize )                     // -x ghost zone
   {
      i_last   = GhostSize;                  // last active cell (adjacent to the -x boundary)
      i_lastm1 = GhostSize + 1;              // one cell further into the interior
      ndepth   = GhostSize - ig;             // 1, 2, ... outward from the boundary
   }
   else                                      // +x ghost zone ( ig >= Nx - GhostSize )
   {
      i_last   = Nx - GhostSize - 1;         // last active cell (adjacent to the +x boundary)
      i_lastm1 = Nx - GhostSize - 2;
      ndepth   = ig - i_last;                // 1, 2, ...
   }

// read Ec (two active cells) and Fc1/2/3 (last active cell) from the prepared interior
   real Ec_last = 0.0, Ec_lastm1 = 0.0, Fc1 = 0.0, Fc2 = 0.0, Fc3 = 0.0;
   for (int v=0; v<NVar_Flu; v++)
   {
      const int f = TFluVarIdxList[v];
      if      ( f == CR_E  ) { Ec_last = Array3D[v][kg][jg][i_last]; Ec_lastm1 = Array3D[v][kg][jg][i_lastm1]; }
      else if ( f == CR_F1 )   Fc1     = Array3D[v][kg][jg][i_last];
      else if ( f == CR_F2 )   Fc2     = Array3D[v][kg][jg][i_last];
      else if ( f == CR_F3 )   Fc3     = Array3D[v][kg][jg][i_last];
   }

// gradient-preserving Ec (linear extrapolation, positive-clamped) + outflow-copy Fc
   const real CR_EC_TINY = (real)1.0e-30;
   const real ec_extrap  = Ec_last + (real)ndepth*( Ec_last - Ec_lastm1 );
   fluid[CR_E ] = FMAX( ec_extrap, CR_EC_TINY );
   fluid[CR_F1] = Fc1;
   fluid[CR_F2] = Fc2;
   fluid[CR_F3] = Fc3;

} // FUNCTION : GradOutflowBC
#endif // #if ( MODEL == HYDRO  &&  defined CR_STREAMING )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CR_Streaming
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CR_Streaming()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == HYDRO  &&  defined CR_STREAMING )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr            = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr     = SetBFieldIC;
#  endif

// static, hollow spherical-shell refinement (enabled by OPT__FLAG_USER with OPT__FLAG_USER_NUM=2;
// the inner/outer shell radii per level are given in "Input__Flag_User")
   Flag_User_Ptr                     = Flag_CR_Streaming;

// the bottleneck test injects CRs through a user boundary condition on the -x face
   if ( CR_Streaming_Test == CR_TEST_BOTTLENECK_1D )
   {
      BC_User_Ptr                    = BottleneckBC;
//    with MHD on, a user fluid BC also requires a matching magnetic-field user BC; the field is
//    uniform (Bx = 1) in this test, so SetBFieldIC already provides the correct ghost-zone values
#     ifdef MHD
      BC_BField_User_Ptr             = SetBFieldIC;
#     endif
   }

// the Sec 4.2.2 shock test fixes the x ghost zones to the initial (inflow) state
   if ( CR_Streaming_Test == CR_TEST_PAPER_SHOCK )
   {
      BC_User_Ptr                    = ShockBC;
#     ifdef MHD
      BC_BField_User_Ptr             = SetBFieldIC;   // uniform Bx = 1 in the ghost zones
#     endif
   }

// the gradient-preserving CR outflow BC (streaming tests) uses the user fluid BC on both x faces
// (OPT__BC_FLU_XM/XP were reset to user in SetParameter); the field is uniform, so SetBFieldIC
// supplies the (copy) B ghost values
   if ( CR_Streaming_GradOutflowBC )
   {
      BC_User_Ptr                    = GradOutflowBC;
#     ifdef MHD
      BC_BField_User_Ptr             = SetBFieldIC;
#     endif
   }
#  endif // #if ( MODEL == HYDRO  &&  defined CR_STREAMING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CR_Streaming
