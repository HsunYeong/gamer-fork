#include "GAMER.h"
#include "TestProb.h"

static void GetCenterOfMass( const double CM_Old[], double CM_New[], const double CM_MaxR );
static void Record_Max( void );
static void Init_Load_StepTable( void );
static void AddNewField_ELBDM_UniformGranule( void );
static void Init_User_ELBDM_UniformGranule( void );
// problem-specific global variables
// =======================================================================================
static FieldIdx_t Idx_Dens0 = Idx_Undefined;  // field index for storing the **initial** density
static double   CM_MaxR;      // maximum radius for determining CM
static double   CM_TolErrR;   // maximum allowed errors for determining CM
static double   Center[3];                    // use CoM coordinate of the whole halo as center
static double   dr_min_prof;                  // bin size of correlation function statistics (minimum size if logarithic bin) (profile)
static double   LogBinRatio_prof;             // ratio of bin size growing rate for logarithmic bin (profile)
static double   RadiusMax_prof;               // maximum radius for correlation function statistics (profile)
static double   dr_min_corr;                  // bin size of correlation function statistics (minimum size if logarithic bin) (correlation)
static double   LogBinRatio_corr;             // ratio of bin size growing rate for logarithmic bin (correlation)
static double   RadiusMax_corr;               // maximum radius for correlation function statistics (correlation)
//static double   PrepTime;                     // time for doing statistics
static bool     ComputeCorrelation;           // flag for compute correlation
static bool     ReComputeCorrelation;         // flag for recompute correlation for restart; use the simulation time of RESTART as initial time for computing time correlation; only available for RESTART
static bool     LogBin_prof;                  // logarithmic bin or not (profile)
static bool     RemoveEmpty_prof;             // remove 0 sample bins; false: Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0 (profile)
static bool     LogBin_corr;                  // logarithmic bin or not (correlation)
static bool     RemoveEmpty_corr;             // remove 0 sample bins; false: Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0 (correlation)
static int      MinLv;                        // do statistics from MinLv to MaxLv
static int      MaxLv;                        // do statistics from MinLv to MaxLv
static int      OutputCorrelationMode;        // output correlation function mode=> 0: constant interval 1: by table
static int      StepInitial;                  // inital step for recording correlation function (OutputCorrelationMode = 0)
static int      StepInterval;                 // interval for recording correlation function (OutputCorrelationMode = 0)
static int      *StepTable;                   // step index table for output correlation function (OutputCorrelationMode = 1)
static bool     Fluid_Periodic_BC_Flag;       // flag for checking the fluid boundary condtion is setup to periodic (0: user defined; 1: periodic)
static char     FilePath_corr[MAX_STRING];    // output path for correlation function text files

static int step_counter;                             // counter for caching consumed step indices
static Profile_with_Sigma_t Prof_Dens_initial;                      // pointer to save initial density profile
static Profile_with_Sigma_t *Prof[] = { &Prof_Dens_initial };
static Profile_t            Correlation_Dens;                       // pointer to save density correlation function
static Profile_t            *Correlation[] = { &Correlation_Dens };
// =======================================================================================

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
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

// only accept OPT__INIT == INIT_BY_RESTART or OPT__INIT == INIT_BY_FILE
   if ( OPT__INIT != INIT_BY_RESTART && OPT__INIT != INIT_BY_FILE )
      Aux_Error( ERROR_INFO, "enforced to accept only OPT__INIT == INIT_BY_RESTART or OPT__INIT == INIT_BY_FILE !!\n" );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == ELBDM  && defined GRAVITY )
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
//                3. Must NOT call any EoS routine here since it hasn't been initialized at this point
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
// ReadPara->Add( "KEY_IN_THE_FILE",          &VARIABLE,                DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "CM_MaxR",                  &CM_MaxR,                 NoMax_double,  Eps_double,       NoMax_double      );
   ReadPara->Add( "CM_TolErrR",               &CM_TolErrR,              0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ComputeCorrelation",       &ComputeCorrelation,      false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Fluid_Periodic_BC_Flag",   &Fluid_Periodic_BC_Flag,  true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "dr_min_corr",              &dr_min_corr,             Eps_double,    Eps_double,       NoMax_double      );
   ReadPara->Add( "LogBinRatio_corr",         &LogBinRatio_corr,        1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "RadiusMax_corr",           &RadiusMax_corr,          Eps_double,    Eps_double,       NoMax_double      );
   ReadPara->Add( "LogBin_corr",              &LogBin_corr,             false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "RemoveEmpty_corr",         &RemoveEmpty_corr,        false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "dr_min_prof",              &dr_min_prof,             Eps_double,    Eps_double,       NoMax_double      );
   ReadPara->Add( "MinLv",                    &MinLv,                   0,             0,                MAX_LEVEL         );
   ReadPara->Add( "MaxLv",                    &MaxLv,                   MAX_LEVEL,     0,                MAX_LEVEL         );
   ReadPara->Add( "OutputCorrelationMode",    &OutputCorrelationMode,   0,             0,                1                 );
   ReadPara->Add( "StepInitial",              &StepInitial,             0,             0,                NoMax_int         );
   ReadPara->Add( "StepInterval",             &StepInterval,            1,             1,                NoMax_int         );
   ReadPara->Add( "FilePath_corr",            FilePath_corr,            Useless_str,   Useless_str,      Useless_str       );
   if ( OPT__INIT == INIT_BY_RESTART )
      ReadPara->Add( "ReComputeCorrelation", &ReComputeCorrelation,     false,         Useless_bool,     Useless_bool      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   if ( CM_TolErrR < 0.0 )   CM_TolErrR = 1.0*amr->dh[MAX_LEVEL];

   if (ComputeCorrelation)
   {
       if ( dr_min_corr <=Eps_double )          dr_min_corr = 1e-3*0.5*amr->BoxSize[0];
       if ( RadiusMax_corr<=Eps_double )        RadiusMax_corr = 0.5*amr->BoxSize[0];
       if ( LogBinRatio_corr<=1.0 )             LogBinRatio_corr = 2.0;

       if ( dr_min_prof <=Eps_double )          dr_min_prof = dr_min_corr;
       RadiusMax_prof                           = RadiusMax_corr * 1.05;   // assigned by Test Problem
       LogBinRatio_prof                         = 1.0;                     // hard-coded by Test Problem (no effect)
       LogBin_prof                              = false;                   // hard-coded by Test Problem
       RemoveEmpty_prof                         = false;                   // hard-coded by Test Problem

       if ( MinLv < 0 ) MinLv = 0;
       if ( MaxLv <= MinLv ) MaxLv = MAX_LEVEL;
       if ( FilePath_corr == "\0" )  sprintf( FilePath_corr, "./" );
       else
       {
          FILE *file_checker = fopen(FilePath_corr, "r");
          if (!file_checker)
             Aux_Error( ERROR_INFO, "File path %s for saving correlation function text files does not exist!! Please create!!\n", FilePath_corr );
          else
             fclose(file_checker);
       }
   }

// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "OPT__INIT = 1 is not supported for this test problem !!\n" );
   // check whether fluid boundary condition in Input__Parameter is set properly
   if ( Fluid_Periodic_BC_Flag )  // use periodic boundary condition
   {
      for ( int direction = 0; direction < 6; direction++ )
      {
         if ( OPT__BC_FLU[direction] != BC_FLU_PERIODIC )
            Aux_Error( ERROR_INFO, "must set periodic BC for fluid --> reset OPT__BC_FLU[%d] to 1 !!\n", direction );
      }
   }
   else  // use user define boundary condition
   {
      for ( int direction = 0; direction < 6; direction++ )
      {
         if ( OPT__BC_FLU[direction] != BC_FLU_USER )
            Aux_Error( ERROR_INFO, "must adopt user defined BC for fluid --> reset OPT__BC_FLU[%d] to 4 !!\n", direction );
      }
   }


// (2) set the problem-specific derived parameters

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__; //2.5

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
      Aux_Message( stdout, "=================================================================================\n" );
      Aux_Message( stdout, "  test problem ID                             = %d\n",         TESTPROB_ID           );
      Aux_Message( stdout, "  CM_MaxR                                     = %13.7e\n",     CM_MaxR               );
      Aux_Message( stdout, "  CM_TolErrR                                  = %13.7e\n",     CM_TolErrR            );
      Aux_Message( stdout, "  fluid periodic boundary condition flag      = %d\n"    , Fluid_Periodic_BC_Flag    );
      Aux_Message( stdout, "  compute correlation                         = %d\n"    , ComputeCorrelation        );
      if (ComputeCorrelation)
      {
         Aux_Message( stdout, "  histogram bin size  (correlation)           = %13.6e\n", dr_min_corr            );
         Aux_Message( stdout, "  log bin ratio       (correlation)           = %13.6e\n", LogBinRatio_corr       );
         Aux_Message( stdout, "  radius maximum      (correlation)           = %13.6e\n", RadiusMax_corr         );
         Aux_Message( stdout, "  use logarithmic bin (correlation)           = %d\n"    , LogBin_corr            );
         Aux_Message( stdout, "  remove empty bin    (correlation)           = %d\n"    , RemoveEmpty_corr       );
         Aux_Message( stdout, "  histogram bin size  (profile)               = %13.6e\n", dr_min_prof            );
         Aux_Message( stdout, "  log bin ratio       (profile, no effect)    = %13.6e\n", LogBinRatio_prof       );
         Aux_Message( stdout, "  radius maximum      (profile, assigned)     = %13.6e\n", RadiusMax_prof         );
         Aux_Message( stdout, "  use logarithmic bin (profile, assigned)     = %d\n"    , LogBin_prof            );
         Aux_Message( stdout, "  remove empty bin    (profile, assigned)     = %d\n"    , RemoveEmpty_prof       );
//         Aux_Message( stdout, "  prepare time                                = %13.6e\n", PrepTime               );
         Aux_Message( stdout, "  minimum level                               = %d\n"    , MinLv                  );
         Aux_Message( stdout, "  maximum level                               = %d\n"    , MaxLv                  );
         Aux_Message( stdout, "  output correlation function mode            = %d\n"    , OutputCorrelationMode  );
         Aux_Message( stdout, "  file path for correlation text file         = %s\n"    , FilePath_corr          );
         if ( OPT__INIT == INIT_BY_RESTART )
         Aux_Message( stdout, "  re-compute correlation using restart time as initial time = %d\n", ReComputeCorrelation );
      }
      Aux_Message( stdout, "=================================================================================\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_StepTable
// Description :  Load the dump table from the file "Input__StepTable"
//-------------------------------------------------------------------------------------------------------
static void Init_Load_StepTable()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_StepTable ...\n" );


   const char FileName[] = "Input__StepTable";

   if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "r" );

   const int MaxLine = 10000;
   char *input_line = NULL;
   size_t len = 0;
   int Trash, line, n;


// allocate the step table
   StepTable = new int [MaxLine];


// skip the header
   getline( &input_line, &len, File );

// begin to read
   for (line=0; line<MaxLine; line++)
   {
      n = getline( &input_line, &len, File );

//    check
      if ( n <= 1 )
         Aux_Error( ERROR_INFO, "incorrect reading at line %d of the file <%s> !!\n", line+2, FileName );

      sscanf( input_line, "%d%d", &Trash, &StepTable[line] );

//    stop the reading
      if ( input_line[0] == 42 )                   // '*' == 42
      {

//       ensure that at least one step index is loaded
         if ( line == 0 )
            Aux_Error( ERROR_INFO, "please provide at least one step index in the step table !!\n" );

         int StepTable_NDump   = line;             // record the number of step indices
         StepTable[line]   = __INT_MAX__;          // set the next step as an extremely large number

         if ( StepTable[line-1] < END_STEP )
         {
            END_STEP          = StepTable[line-1];    // reset the ending time as the time of the last step

            if ( MPI_Rank == 0 )
               Aux_Message( stdout, "NOTE : the END_STEP is reset to the time of the last step index = %de\n",
                            END_STEP );
         }


//       verify the loaded dump table
         for (int t=1; t<=line; t++)
         {
            if ( StepTable[t] < StepTable[t-1] )
               Aux_Error( ERROR_INFO, "values recorded in \"%s\" must be monotonically increasing !!\n",
                          FileName );
         }

         break;

      } // if ( input_line[0] == 42 )
   } // for (line=0; line<MaxLine; line++)


   if ( line == MaxLine )
      Aux_Error( ERROR_INFO, "please prepare a symbol * in the end of the file <%s> !!\n", FileName );


   fclose( File );

   if ( input_line != NULL )     free( input_line );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_StepTable ... done\n" );

} // FUNCTION : Init_Load_StepTable

//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_ELBDM_UniformGranule
// Description :  Add the problem-specific fields
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
static void AddNewField_ELBDM_UniformGranule(void)
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   Idx_Dens0 = AddField( "Dens0", NORMALIZE_NO, INTERP_FRAC_NO );
//   if ( MPI_Rank == 0 )   printf("Idx_Dens0 = %d \n", Idx_Dens0);
#  endif

} // FUNCTION : AddNewField_ELBDM_Halo_Stability_Test

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_ELBDM_UniformGranule
// Description :  Store the initial density
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Init_User_ELBDM_UniformGranule(void)
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
//    store the initial density in both Sg so that we don't have to worry about which Sg to be used
//    a. for restart and ReComputeCorrelation disabled, the initial density has already been loaded and we just need to copy the data to another Sg
      if ( ( OPT__INIT == INIT_BY_RESTART ) && ( !ReComputeCorrelation ) ) {
         const real Dens0 = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i];

         amr->patch[ 1-amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
      }

//    b. for starting a new simulation, we must copy the initial density to both Sg
      else {
         const real Dens0 = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];

         amr->patch[   amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
         amr->patch[ 1-amr->FluSg[lv] ][lv][PID]->fluid[Idx_Dens0][k][j][i] = Dens0;
      }
   }

   if (ComputeCorrelation)
   {
       step_counter = 0;
       const double InitialTime = Time[0];
       if ( OutputCorrelationMode == 0)
       {
          if ( MPI_Rank==0 ) Aux_Message( stdout, "StepInitial = %d ; StepInterval = %d \n", StepInitial, StepInterval);
       }
       else if ( OutputCorrelationMode == 1)
          Init_Load_StepTable();
       if ( MPI_Rank == 0 )  Aux_Message( stdout, "InitialTime = %13.6e \n", InitialTime );

       // compute the enter position for passive field
       if ( MPI_Rank == 0 )  Aux_Message( stdout, "Calculate halo center for passive field:\n");

       Record_Max();
       if ( MPI_Rank == 0 )  Aux_Message( stdout, "Center of passive field is ( %14.11e,%14.11e,%14.11e )\n", Center[0], Center[1], Center[2] );
       // commpute density profile for passive field;
       if ( MPI_Rank == 0 )  Aux_Message( stdout, "Calculate density profile for passive field:\n");

       const long TVar[] = {BIDX(Idx_Dens0)};
       Aux_ComputeProfile_with_Sigma( Prof, Center, RadiusMax_prof, dr_min_prof, LogBin_prof, LogBinRatio_prof, RemoveEmpty_prof, TVar, 1, MinLv, MaxLv, PATCH_LEAF, InitialTime );

       if ( MPI_Rank == 0 )
       {
          char Filename[MAX_STRING];
          sprintf( Filename, "%s/initial_profile_with_Sigma.txt", FilePath_corr );
          FILE *output_initial_prof = fopen(Filename, "w");
          fprintf( output_initial_prof, "#%19s  %21s  %21s  %21s  %11s\n", "Radius", "Dens", "Dens_Sigma" , "Weighting", "Cell_Number");
          for (int b=0; b<Prof[0]->NBin; b++)
             fprintf( output_initial_prof, "%20.14e  %21.14e  %21.14e  %21.14e  %11ld\n",
                       Prof[0]->Radius[b], Prof[0]->Data[b], Prof[0]->Data_Sigma[b], Prof[0]->Weight[b], Prof[0]->NCell[b] );
          fclose(output_initial_prof);
       }
   }
#  endif

} // FUNCTION : Init_User_ELBDM_UniformGranule

//-------------------------------------------------------------------------------------------------------
// Function    :  GetCenterOfMass
// Description :  Record the center of mass (CM)
//
// Note        :  1. Invoked by Record_EridanusII() recursively
//                2. Only include cells within CM_MaxR from CM_Old[] when updating CM
//
// Parameter   :  CM_Old[] : Previous CM
//                CM_New[] : New CM to be returned
//                CM_MaxR  : Maximum radius to compute CM
//
// Return      :  CM_New[]
//-------------------------------------------------------------------------------------------------------
void GetCenterOfMass( const double CM_Old[], double CM_New[], const double CM_MaxR )
{

   const double CM_MaxR2          = SQR( CM_MaxR );
   const double HalfBox[3]        = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic          = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC );
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const real   MinEntr_No        = -1.0;
   const bool   DE_Consistency_No = false;

   int   *PID0List = NULL;
   double M_ThisRank, MR_ThisRank[3], M_AllRank, MR_AllRank[3];
   real (*TotalDens)[PS1][PS1][PS1];

   M_ThisRank = 0.0;
   for (int d=0; d<3; d++)    MR_ThisRank[d] = 0.0;


   for (int lv=0; lv<NLEVEL; lv++)
   {
//    get the total density on grids
      TotalDens = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      PID0List  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, _DENS, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

      delete [] PID0List;

//    calculate the center of mass
      const double dh = amr->dh[lv];
      const double dv = CUBE( dh );

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         double x, y, z, dx, dy, dz;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;  dz = z - CM_Old[2];
                                       if ( Periodic ) {
                                          if      ( dz > +HalfBox[2] )  {  z -= amr->BoxSize[2];  dz -= amr->BoxSize[2];  }
                                          else if ( dz < -HalfBox[2] )  {  z += amr->BoxSize[2];  dz += amr->BoxSize[2];  }
                                       }
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;  dy = y - CM_Old[1];
                                       if ( Periodic ) {
                                          if      ( dy > +HalfBox[1] )  {  y -= amr->BoxSize[1];  dy -= amr->BoxSize[1];  }
                                          else if ( dy < -HalfBox[1] )  {  y += amr->BoxSize[1];  dy += amr->BoxSize[1];  }
                                       }
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;  dx = x - CM_Old[0];
                                       if ( Periodic ) {
                                          if      ( dx > +HalfBox[0] )  {  x -= amr->BoxSize[0];  dx -= amr->BoxSize[0];  }
                                          else if ( dx < -HalfBox[0] )  {  x += amr->BoxSize[0];  dx += amr->BoxSize[0];  }
                                       }

//          only include cells within CM_MaxR
            const double R2 = SQR(dx) + SQR(dy) + SQR(dz);
            if ( R2 < CM_MaxR2 )
            {
               const double dm = TotalDens[PID][k][j][i]*dv;

               M_ThisRank     += dm;
               MR_ThisRank[0] += dm*x;
               MR_ThisRank[1] += dm*y;
               MR_ThisRank[2] += dm*z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// collect data from all ranks to calculate the CM
// --> note that all ranks will get CM_New[]
   MPI_Allreduce( &M_ThisRank, &M_AllRank, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( MR_ThisRank, MR_AllRank, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)    CM_New[d] = MR_AllRank[d] / M_AllRank;

// map the new CM back to the simulation domain
   if ( Periodic )
   for (int d=0; d<3; d++)
   {
      if      ( CM_New[d] >= amr->BoxSize[d] )  CM_New[d] -= amr->BoxSize[d];
      else if ( CM_New[d] < 0.0              )  CM_New[d] += amr->BoxSize[d];

   }

   for (int d=0; d<3; d++)
      if ( CM_New[d] >= amr->BoxSize[d]  ||  CM_New[d] < 0.0 )
         Aux_Error( ERROR_INFO, "CM_New[%d] = %14.7e lies outside the domain !!\n", d, CM_New[d] );

} // FUNCTION : GetCenterOfMass


//-------------------------------------------------------------------------------------------------------
// Function    :  Record_Max
// Description :  Record the maximum density and center coordinates
//
// Note        :  1. It will also record the real and imaginary parts associated with the maximum density
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filenames are fixed to "Record__MaxDens" and "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Record_Max( )
{
   const char filename_max_dens[] = "Record__MaxDens";
   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 10;

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3], real_loc, imag_loc;
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];


// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID]->son != -1 )  continue;

      for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
      for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
      for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

         dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
         pote = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];

         if ( dens > max_dens_loc )
         {
            max_dens_loc        = dens;
            real_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i];
            imag_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i];
            max_dens_pos_loc[0] = x;
            max_dens_pos_loc[1] = y;
            max_dens_pos_loc[2] = z;
         }

         if ( pote < min_pote_loc )
         {
            min_pote_loc        = pote;
            min_pote_pos_loc[0] = x;
            min_pote_pos_loc[1] = y;
            min_pote_pos_loc[2] = z;
         }
      }}}
   }


// gather data to the root rank
   send[0] = max_dens_loc;
   send[1] = real_loc;
   send[2] = imag_loc;
   send[3] = max_dens_pos_loc[0];
   send[4] = max_dens_pos_loc[1];
   send[5] = max_dens_pos_loc[2];
   send[6] = min_pote_loc;
   send[7] = min_pote_pos_loc[0];
   send[8] = min_pote_pos_loc[1];
   send[9] = min_pote_pos_loc[2];

   MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );


// record the maximum density and center coordinates
   double max_dens      = -__DBL_MAX__;
   double min_pote      = +__DBL_MAX__;
   int    max_dens_rank = -1;
   int    min_pote_rank = -1;

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         if ( recv[r][0] > max_dens )
         {
            max_dens      = recv[r][0];
            max_dens_rank = r;
         }

         if ( recv[r][6] < min_pote )
         {
            min_pote      = recv[r][6];
            min_pote_rank = r;
         }
      }

      if ( max_dens_rank < 0  ||  max_dens_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect max_dens_rank (%d) !!\n", max_dens_rank );

      if ( min_pote_rank < 0  ||  min_pote_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect min_pote_rank (%d) !!\n", min_pote_rank );

      static bool FirstTime = true;

      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_max_dens) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_max_dens );
         else
         {
            FILE *file_max_dens = fopen( filename_max_dens, "w" );
            fprintf( file_max_dens, "#%19s   %10s   %14s   %14s   %14s\n", "Time", "Step", "Dens", "Real", "Imag" );
            fclose( file_max_dens );
         }

         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter", "CM_x", "CM_y", "CM_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_max_dens = fopen( filename_max_dens, "a" );
      fprintf( file_max_dens, "%20.14e   %10ld   %14.7e   %14.7e   %14.7e\n",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2] );
      fclose( file_max_dens );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][3], recv[max_dens_rank][4], recv[max_dens_rank][5],
                              recv[min_pote_rank][6], recv[min_pote_rank][7], recv[min_pote_rank][8], recv[min_pote_rank][9] );
      fclose( file_center );
   } // if ( MPI_Rank == 0 )


// compute the center of mass until convergence
   const double TolErrR2 = SQR( CM_TolErrR );
   const int    NIterMax = 10;

   double dR2, CM_Old[3], CM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    CM_Old[d] = recv[max_dens_rank][3+d];

   MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   while ( true )
   {
      GetCenterOfMass( CM_Old, CM_New, CM_MaxR );

      dR2 = SQR( CM_Old[0] - CM_New[0] )
          + SQR( CM_Old[1] - CM_New[1] )
          + SQR( CM_Old[2] - CM_New[2] );
      NIter ++;

      if ( dR2 <= TolErrR2  ||  NIter >= NIterMax )
         break;
      else
         memcpy( CM_Old, CM_New, sizeof(double)*3 );
   }

   if ( MPI_Rank == 0 )
   {
      if ( dR2 > TolErrR2 )
         Aux_Message( stderr, "WARNING : dR (%13.7e) > CM_TolErrR (%13.7e) !!\n", sqrt(dR2), CM_TolErrR );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n", NIter, CM_New[0], CM_New[1], CM_New[2] );
      fclose( file_center );
   }

   delete [] recv;

} // FUNCTION : Record_Max



#endif // #if ( MODEL == ELBDM )


static void Do_COM_and_CF( void )
{
   Record_Max();

// Compute correlation if ComputeCorrelation flag is true
   if (ComputeCorrelation)
   {
      if ( ((OutputCorrelationMode==1) && (Step==StepTable[step_counter])) || ((OutputCorrelationMode==0) && (Step>=StepInitial) && (((Step-StepInitial)%StepInterval)==0)) )
      {
         const long TVar[] = {_DENS};
         Aux_ComputeCorrelation( Correlation, (const Profile_with_Sigma_t**)Prof, Center, RadiusMax_corr, dr_min_corr, LogBin_corr, LogBinRatio_corr,
                                 RemoveEmpty_corr, TVar, 1, MinLv, MaxLv, PATCH_LEAF, Time[0], dr_min_prof);

         if ( MPI_Rank == 0 )
         {
            char Filename[MAX_STRING];
            sprintf( Filename, "%s/correlation_function_t=%.4e.txt", FilePath_corr, Time[0] );
            FILE *output_correlation = fopen(Filename, "w");
            fprintf( output_correlation, "#%19s  %21s  %21s  %11s\n", "Radius", "Correlation_Function", "Weighting", "Cell_Number");
            for (int b=0; b<Correlation[0]->NBin; b++)
                fprintf( output_correlation, "%20.14e  %21.14e  %21.14e  %11ld\n",
                         Correlation[0]->Radius[b], Correlation[0]->Data[b], Correlation[0]->Weight[b], Correlation[0]->NCell[b] );
            fclose(output_correlation);
         }
         // accumulate the step counter
         step_counter ++;
      }
   }  // end of if ComputeCorrelation
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Template
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_UniformGranule()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-// 1. define a user-specified function (example functions are given below)
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
//   Init_Function_User_Ptr = SetGridIC;
   Init_Field_User_Ptr    = AddNewField_ELBDM_UniformGranule;
   Init_User_Ptr          = Init_User_ELBDM_UniformGranule;
   Aux_Record_User_Ptr    = Do_COM_and_CF;
//   BC_User_Ptr            = BC;
#  ifdef MHD
//   Init_Function_BField_User_Ptr  = SetBFieldIC;
#  endif
// comment out Init_ByFile_User_Ptr to use the default
//   Init_ByFile_User_Ptr           = NULL; // option: OPT__INIT=3;             example: Init/Init_ByFile.cpp -> Init_ByFile_Default()
//   Init_Field_User_Ptr            = NULL; // set NCOMP_PASSIVE_USER;          example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
//   Flag_User_Ptr                  = NULL; // option: OPT__FLAG_USER;          example: Refine/Flag_User.cpp
//   Mis_GetTimeStep_User_Ptr       = NULL; // option: OPT__DT_USER;            example: Miscellaneous/Mis_GetTimeStep_User.cpp
//   BC_User_Ptr                    = NULL; // option: OPT__BC_FLU_*=4;         example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
#  ifdef MHD
//   BC_BField_User_Ptr             = NULL; // option: OPT__BC_FLU_*=4;
#  endif
//   Flu_ResetByUser_Func_Ptr       = NULL; // option: OPT__RESET_FLUID;        example: Fluid/Flu_ResetByUser.cpp
//   Init_DerivedField_User_Ptr     = NULL; // option: OPT__OUTPUT_USER_FIELD;  example: Fluid/Flu_DerivedField_User.cpp
//   Output_User_Ptr                = NULL; // option: OPT__OUTPUT_USER;        example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
//   Aux_Record_User_Ptr            = NULL; // option: OPT__RECORD_USER;        example: Auxiliary/Aux_Record_User.cpp
//   Init_User_Ptr                  = NULL; // option: none;                    example: none
//   End_User_Ptr                   = NULL; // option: none;                    example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
//   Init_ExtAcc_Ptr                = NULL; // option: OPT__EXT_ACC;            example: SelfGravity/CPU_Gravity/CPU_ExtAcc_PointMass.cpp
//   End_ExtAcc_Ptr                 = NULL;
//   Init_ExtPot_Ptr                = NULL; // option: OPT__EXT_POT;            example: SelfGravity/CPU_Poisson/CPU_ExtPot_PointMass.cpp
//   End_ExtPot_Ptr                 = NULL;
//   Poi_AddExtraMassForGravity_Ptr = NULL; // option: OPT__GRAVITY_EXTRA_MASS; example: none
//   Poi_UserWorkBeforePoisson_Ptr  = NULL; // option: none;                    example: SelfGravity/Poi_UserWorkBeforePoisson.cpp
#  endif
#  if ( EOS == EOS_USER )
//   EoS_Init_Ptr                   = NULL; // option: EOS in the Makefile;     example: EoS/User_Template/CPU_EoS_User_Template.cpp
//   EoS_End_Ptr                    = NULL;
#  endif
#  endif // #if ( MODEL == HYDRO )
//   Src_Init_User_Ptr              = NULL; // option: SRC_USER;                example: SourceTerms/User_Template/CPU_Src_User_Template.cpp



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_DiskHeating
