#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double WLMDwarfGalaxy_PEHeatingScaleLength;  // scale length of the photoelectric heating rate distribution
static double WLMDwarfGalaxy_PEHeatingScaleHeight;  // scale height of the photoelectric heating rate distribution
static double WLMDwarfGalaxy_PEHeatingCoreLength;   // core radius of the photoelectric heating rate distribution
static double WLMDwarfGalaxy_PEHeatingCoreHeight;   // core height of the photoelectric heating rate distribution
static double WLMDwarfGalaxy_PEHeatingRate0;        // central photoelectric heating rate (in erg/cm^3/s/n_H)
static double WLMDwarfGalaxy_PEHeatingRateBg;       // background photoelectric heating rate (in erg/cm^3/s/n_H)
static double WLMDwarfGalaxy_CenterOfMass_X = -1.0; // x coordinate of the center of mass of the system
static double WLMDwarfGalaxy_CenterOfMass_Y = -1.0; // y coordinate of the center of mass of the system
static double WLMDwarfGalaxy_CenterOfMass_Z = -1.0; // z coordinate of the center of mass of the system
// =======================================================================================


// problem-specific function prototypes
bool Flag_WLMDwarfGalaxy( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );




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


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifndef SUPPORT_GRACKLE
   Aux_Error( ERROR_INFO, "SUPPORT_GRACKLE must be enabled !!\n" );
#  endif

#  ifndef STAR_FORMATION
   Aux_Error( ERROR_INFO, "STAR_FORMATION must be enabled !!\n" );
#  endif

#  ifndef FEEDBACK
   Aux_Error( ERROR_INFO, "FEEDBACK must be enabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_FLU* = 1) for this test !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_POT = 1) for this test !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FILE  &&  amr->Par->Init != PAR_INIT_BY_FILE )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 3 (by FILE) !!\n" );
#  endif

#  ifdef FEEDBACK
   if ( !FB_RESOLVED_SNEII )
      Aux_Error( ERROR_INFO, "FB_RESOLVED_SNEII must be enabled !!\n" );
#  endif

#  ifdef SUPPORT_GRACKLE
   if ( !GRACKLE_ACTIVATE )
      Aux_Error( ERROR_INFO, "GRACKLE_ACTIVATE must be enabled for this test !!\n" );

   if ( !GRACKLE_METAL )
      Aux_Error( ERROR_INFO, "GRACKLE_METAL must be enabled for this test !!\n" );

   if ( GRACKLE_PE_HEATING  &&  GRACKLE_PE_HEATING_RATE > 0.0 )
      Aux_Error( ERROR_INFO, "GRACKLE_PE_HEATING_RATE must be 0.0 for this test !!\n" );

   if ( !GRACKLE_USE_V_HEATING_RATE )
      Aux_Error( ERROR_INFO, "GRACKLE_USE_V_HEATING_RATE must be enabled for this test !!\n" );
#  endif // #ifdef SUPPORT_GRACKLE

   if ( !OPT__RECORD_USER )
      Aux_Error( ERROR_INFO, "OPT__RECORD_USER must be enabled for this test !!\n" );

   if ( OPT__RECORD_CENTER )
      Aux_Error( ERROR_INFO, "OPT__RECORD_CENTER is redundant and must be disabled as OPT__RECORD_USER is enabled for this test !!\n" );

   if ( MPI_Rank == 0 )
   {
#     if ( FLU_SCHEME != MHM )
         Aux_Message( stderr, "WARNING : it's recommended to adopt the MHM scheme for this test !!\n" );
#     endif

#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif

#     if ( !defined MHD  &&  RSOLVER != HLLC )
         Aux_Message( stderr, "WARNING : it's recommended to adopt the HLLC Riemann solver for this test !!\n" );
#     endif

#     if ( MODEL == HYDRO )
      if ( MINMOD_COEFF > 1.5 )
         Aux_Message( stderr, "WARNING : it's recommended to set MINMOD_COEFF <= 1.5 for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadInputTestProb
// Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
//
// Note        :  1. Invoked by SetParameter() to read parameters
//                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
//                3. If there is no problem-specific runtime parameter to load, add at least one parameter
//                   to prevent an empty structure in HDF5_Output_t
//                   --> Example:
//                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
//
// Parameter   :  load_mode      : Mode for loading parameters
//                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
//                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
//                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
//                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
{

#  ifndef SUPPORT_HDF5
   if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
#  endif

   if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
   if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// --> LOAD_PARA() is defined in "include/TestProb.h"
// ********************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",                       &VARIABLE,                               DEFAULT,      MIN,              MAX               );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "WLMDwarfGalaxy_PEHeatingScaleLength",   &WLMDwarfGalaxy_PEHeatingScaleLength,    0.4,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "WLMDwarfGalaxy_PEHeatingScaleHeight",   &WLMDwarfGalaxy_PEHeatingScaleHeight,    0.3,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "WLMDwarfGalaxy_PEHeatingCoreLength",    &WLMDwarfGalaxy_PEHeatingCoreLength,     0.4,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "WLMDwarfGalaxy_PEHeatingCoreHeight",    &WLMDwarfGalaxy_PEHeatingCoreHeight,     0.0,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "WLMDwarfGalaxy_PEHeatingRate0",         &WLMDwarfGalaxy_PEHeatingRate0,          2.6e-27,      0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "WLMDwarfGalaxy_PEHeatingRateBg",        &WLMDwarfGalaxy_PEHeatingRateBg,         2.1e-29,      0.0,              NoMax_double      );

} // FUNCITON : LoadInputTestProb



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
// (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// check the runtime parameters

// convert to code units
   WLMDwarfGalaxy_PEHeatingScaleLength *= Const_kpc / UNIT_L;
   WLMDwarfGalaxy_PEHeatingScaleHeight *= Const_kpc / UNIT_L;
   WLMDwarfGalaxy_PEHeatingCoreLength  *= Const_kpc / UNIT_L;
   WLMDwarfGalaxy_PEHeatingCoreHeight  *= Const_kpc / UNIT_L;


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  1000.0*Const_Myr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                     = %d\n",                    TESTPROB_ID                                            );
      Aux_Message( stdout, "  WLMDwarfGalaxy_PEHeatingScaleLength = %13.7e kpc\n",            WLMDwarfGalaxy_PEHeatingScaleLength * UNIT_L/Const_kpc );
      Aux_Message( stdout, "  WLMDwarfGalaxy_PEHeatingScaleHeight = %13.7e kpc\n",            WLMDwarfGalaxy_PEHeatingScaleHeight * UNIT_L/Const_kpc );
      Aux_Message( stdout, "  WLMDwarfGalaxy_PEHeatingCoreLength  = %13.7e kpc\n",            WLMDwarfGalaxy_PEHeatingCoreLength  * UNIT_L/Const_kpc );
      Aux_Message( stdout, "  WLMDwarfGalaxy_PEHeatingCoreHeight  = %13.7e kpc\n",            WLMDwarfGalaxy_PEHeatingCoreHeight  * UNIT_L/Const_kpc );
      Aux_Message( stdout, "  WLMDwarfGalaxy_PEHeatingRate0       = %13.7e erg/cm^3/s/n_H\n", WLMDwarfGalaxy_PEHeatingRate0                          );
      Aux_Message( stdout, "  WLMDwarfGalaxy_PEHeatingRateBg      = %13.7e erg/cm^3/s/n_H\n", WLMDwarfGalaxy_PEHeatingRateBg                         );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Initialize the fluid for the WLM dwarf galaxy test
//
// Note        :  1. NOT needed beacuse the fluid is initialized by file
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
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

   Aux_Error( ERROR_INFO, "OPT__INIT = 1 is not supported for this test problem !!\n" );

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_WLMDwarfGalaxy
// Description :  Add the problem-specific grid fields
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_WLMDwarfGalaxy()
{

// add the metallicity field only if it has not been done
// --> since Grackle may already add this field automatically when GRACKLE_METAL is enabled
// --> also note that "Idx_Metal" has been predefined in Field.h
   if ( Idx_Metal == Idx_Undefined )
      Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_NO, INTERP_FRAC_YES );

   if ( Idx_Dust == Idx_Undefined )
      Idx_Dust = AddField( "Dust", FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_NO, INTERP_FRAC_YES );

} // FUNCTION : AddNewField_WLMDwarfGalaxy



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_WLMDwarfGalaxy
// Description :  Add the problem-specific particle attributes
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
//                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
//                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
//                3. Pre-declared attribute indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_WLMDwarfGalaxy()
{

// "Idx_ParMetalFrac" has been predefined in Field.h
   if ( Idx_ParMetalFrac == Idx_Undefined )
      Idx_ParMetalFrac = AddParticleAttributeFlt( "ParMetalFrac" );

} // FUNCTION : AddNewParticleAttribute_WLMDwarfGalaxy



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_WLMDwarfGalaxy
// Description :  Record the center for WLMDwarfGalaxy
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//                4. This is the same as Aux_Record_Center
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_WLMDwarfGalaxy()
{

   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__Center", OUTPUT_DIR );


// 1. Maximum fluid density in HYDRO/ELBDM
   Extrema_t Max_Dens;
   Max_Dens.Field     = _DENS;
   Max_Dens.Radius    = __FLT_MAX__; // entire domain
   Max_Dens.Center[0] = amr->BoxCenter[0];
   Max_Dens.Center[1] = amr->BoxCenter[1];
   Max_Dens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_Dens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );


#  ifdef PARTICLE
// 2. Maximum particle density
   Extrema_t Max_ParDens;
   Max_ParDens.Field     = _PAR_DENS;
   Max_ParDens.Radius    = __FLT_MAX__; // entire domain
   Max_ParDens.Center[0] = amr->BoxCenter[0];
   Max_ParDens.Center[1] = amr->BoxCenter[1];
   Max_ParDens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_ParDens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );


// 3. Maximum total density including fluid density and particle density
   Extrema_t Max_TotDens;
   Max_TotDens.Field     = _TOTAL_DENS;
   Max_TotDens.Radius    = __FLT_MAX__; // entire domain
   Max_TotDens.Center[0] = amr->BoxCenter[0];
   Max_TotDens.Center[1] = amr->BoxCenter[1];
   Max_TotDens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_TotDens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );
#  endif


#  ifdef GRAVITY
// 4. Minimum gravitational potential
   Extrema_t Min_Pote;
   Min_Pote.Field     = _POTE;
   Min_Pote.Radius    = __FLT_MAX__; // entire domain
   Min_Pote.Center[0] = amr->BoxCenter[0];
   Min_Pote.Center[1] = amr->BoxCenter[1];
   Min_Pote.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Min_Pote, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );
#  endif


// 5. Center of mass for the total density field
// find the center of mass for the total density field
   double CoM_Coord[3];
   double FinaldR;
   int    FinalNIter;
   Aux_FindWeightedAverageCenter( CoM_Coord, amr->BoxCenter, __FLT_MAX__, 0.0, _TOTAL_DENS,  __FLT_MAX__, 1, &FinaldR, &FinalNIter );


// 6. Center of mass for the gas density field
// find the center of mass for the gas density field
   double GasCoM_Coord[3];
   double GasFinaldR;
   int    GasFinalNIter;
   Aux_FindWeightedAverageCenter( GasCoM_Coord, amr->BoxCenter, __FLT_MAX__, 0.0, _DENS,  __FLT_MAX__, 1, &GasFinaldR, &GasFinalNIter );


// Output the center to file
   if ( MPI_Rank == 0 )
   {
//    Output the header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );
         else
         {
            FILE *File = fopen( FileName, "w" );
            fprintf( File, "#%19s  %10s  %14s  %14s  %14s  %14s",
                           "Time", "Step", "MaxDens", "MaxDens_x", "MaxDens_y", "MaxDens_z" );

#           ifdef PARTICLE
            fprintf( File, "  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s",
                           "MaxParDens", "MaxParDens_x", "MaxParDens_y", "MaxParDens_z",
                           "MaxTotalDens", "MaxTotalDens_x", "MaxTotalDens_y", "MaxTotalDens_z" );
#           endif

#           ifdef GRAVITY
            fprintf( File, "  %14s  %14s  %14s  %14s",
                           "MinPote", "MinPote_x", "MinPote_y", "MinPote_z" );
#           endif

            fprintf( File, "  %14s  %14s  %14s  %14s  %14s",
                           "Final_NIter", "Final_dR", "CoM_x", "CoM_y", "CoM_z" );

            fprintf( File, "  %14s  %14s  %14s  %14s  %14s",
                           "GasFl_NIter", "GasFl_dR", "GasCoM_x", "GasCoM_y", "GasCoM_z" );

            fprintf( File, "\n" );
            fclose( File );
         }

         FirstTime = false;
      } // if ( FirstTime )

      FILE *File = fopen( FileName, "a" );
      fprintf( File, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e",
                     Time[0], Step, Max_Dens.Value, Max_Dens.Coord[0], Max_Dens.Coord[1], Max_Dens.Coord[2] );

#     ifdef PARTICLE
      fprintf( File, "  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
                     Max_ParDens.Value, Max_ParDens.Coord[0], Max_ParDens.Coord[1], Max_ParDens.Coord[2],
                     Max_TotDens.Value, Max_TotDens.Coord[0], Max_TotDens.Coord[1], Max_TotDens.Coord[2] );
#     endif

#     ifdef GRAVITY
      fprintf( File, "  %14.7e  %14.7e  %14.7e  %14.7e",
                     Min_Pote.Value, Min_Pote.Coord[0], Min_Pote.Coord[1], Min_Pote.Coord[2] );
#     endif

      fprintf( File, "  %14d  %14.7e  %14.7e  %14.7e  %14.7e",
                     FinalNIter, FinaldR, CoM_Coord[0], CoM_Coord[1], CoM_Coord[2] );

      fprintf( File, "  %14d  %14.7e  %14.7e  %14.7e  %14.7e",
                     GasFinalNIter, GasFinaldR, GasCoM_Coord[0], GasCoM_Coord[1], GasCoM_Coord[2] );

      fprintf( File, "\n" );
      fclose( File );
   } // if ( MPI_Rank == 0 )

// pass the values to set the volumetric heating rate
   WLMDwarfGalaxy_CenterOfMass_X = GasCoM_Coord[0];
   WLMDwarfGalaxy_CenterOfMass_Y = GasCoM_Coord[1];
   WLMDwarfGalaxy_CenterOfMass_Z = GasCoM_Coord[2];

} // FUNCTION : Aux_Record_WLMDwarfGalaxy



#  ifdef SUPPORT_GRACKLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_vHeatingRate_WLMDwarfGalaxy
// Description :  Function to set Grackle's volumetric heating rate for WLMDwarfGalaxy
//
// Note        :  1. Invoked by Grackle_Prepare() using the function pointer
//                   "Grackle_vHeatingRate_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Returned rate should be in the unit of erg s^-1 cm^-3
//
// Parameter   :  x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                n_H      : Hydrogen number density, should be in the unit of cm^-3
//
// Return      :  volumetric_heating_rate
//-------------------------------------------------------------------------------------------------------
real_che Grackle_vHeatingRate_WLMDwarfGalaxy( const double x, const double y, const double z, const double Time, const double n_H )
{
   if ( WLMDwarfGalaxy_CenterOfMass_X > amr->BoxEdgeR[0]  ||  WLMDwarfGalaxy_CenterOfMass_X < amr->BoxEdgeL[0]  ||
        WLMDwarfGalaxy_CenterOfMass_Y > amr->BoxEdgeR[1]  ||  WLMDwarfGalaxy_CenterOfMass_Y < amr->BoxEdgeL[1]  ||
        WLMDwarfGalaxy_CenterOfMass_Z > amr->BoxEdgeR[2]  ||  WLMDwarfGalaxy_CenterOfMass_Z < amr->BoxEdgeL[2] )
      Aux_Error( ERROR_INFO, "Incorrect CenterOfMass [%14.7e, %14.7e, %14.7e] !!\n",
                 WLMDwarfGalaxy_CenterOfMass_X, WLMDwarfGalaxy_CenterOfMass_Y, WLMDwarfGalaxy_CenterOfMass_Z );

   const double   R_to_cen = sqrt( SQR( x - WLMDwarfGalaxy_CenterOfMass_X ) + SQR( y - WLMDwarfGalaxy_CenterOfMass_Y ) );
   const double   H_to_cen = fabs(      z - WLMDwarfGalaxy_CenterOfMass_Z                                              );

   const double   R_normed = ( R_to_cen - WLMDwarfGalaxy_PEHeatingCoreLength ) / WLMDwarfGalaxy_PEHeatingScaleLength;
   const double   H_normed = ( H_to_cen - WLMDwarfGalaxy_PEHeatingCoreHeight ) / WLMDwarfGalaxy_PEHeatingScaleHeight;

   const double   factor_R = exp( - fmax( R_normed, 0.0 ) );
   const double   factor_H = exp( - fmax( H_normed, 0.0 ) );

   const double   PEHeatingRate_at_R_H    = WLMDwarfGalaxy_PEHeatingRate0 * factor_R * factor_H;
   const real_che volumetric_heating_rate = n_H * fmax( PEHeatingRate_at_R_H, WLMDwarfGalaxy_PEHeatingRateBg );

   return volumetric_heating_rate;

} // FUNCTION : Grackle_vHeatingRate_WLMDwarfGalaxy
#endif // #ifdef SUPPORT_GRACKLE
#endif // #if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_WLMDwarfGalaxy_IsolatedGalaxy
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_WLMDwarfGalaxy()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   Init_Field_User_Ptr           = AddNewField_WLMDwarfGalaxy;
   Par_Init_Attribute_User_Ptr   = AddNewParticleAttribute_WLMDwarfGalaxy;
   Flag_User_Ptr                 = Flag_WLMDwarfGalaxy;
   Aux_Record_User_Ptr           = Aux_Record_WLMDwarfGalaxy;
#  ifdef SUPPORT_GRACKLE
   Grackle_vHeatingRate_User_Ptr = Grackle_vHeatingRate_WLMDwarfGalaxy;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_WLMDwarfGalaxy
