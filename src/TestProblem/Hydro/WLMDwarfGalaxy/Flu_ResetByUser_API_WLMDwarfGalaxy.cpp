#include "GAMER.h"

extern bool WLMDwarfGalaxy_UniformB;

void MHD_Init_BField_ByVecPot_File_WLMDwarfGalaxy( const int B_lv );

//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API_WLMDwarfGalaxy
// Description :  Default API for resetting the fluid array
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked by EvolveLevel() using the function pointer "Flu_ResetByUser_API_Ptr"
//                   --> This function pointer may be reset by a test problem initializer, in which case
//                       this funtion will become useless
//                3. Currently NOT applied to the input uniform array
//                   --> Init_ByFile() does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//
// Parameter   :  lv      : Target refinement level
//                FluSg   : Target fluid sandglass
//                MagSg   : Target B field sandglass
//                TimeNew : Current physical time (system has been updated from TimeOld to TimeNew in EvolveLevel())
//                dt      : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API_WLMDwarfGalaxy( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt )
{

// check
#  ifdef MHD
   if ( Flu_ResetByUser_Func_Ptr == NULL  &&  MHD_ResetByUser_BField_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_Func_Ptr == NULL and MHD_ResetByUser_BField_Ptr == NULL for OPT__RESET_FLUID !!\n" );
#  else
   if ( Flu_ResetByUser_Func_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_Func_Ptr == NULL for OPT__RESET_FLUID !!\n" );
#  endif


   const bool   ResetFlu  = ( Flu_ResetByUser_Func_Ptr   != NULL );
#  ifdef MHD
   const bool   ResetMag  = ( MHD_ResetByUser_BField_Ptr != NULL );
   const bool   UseVecPot = ( MHD_ResetByUser_VecPot_Ptr != NULL );
#  endif
   const double dh   = amr->dh[lv];
   const double dh_2 = 0.5*dh;
#  ifdef OPENMP
   const int    NT   = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int    NT   = 1;
#  endif
   const int    Der_NP              = 8;
   const bool   IntPhase_No         = false;
   const bool   DE_Consistency_No   = false;
   const real   MinDens_No          = -1.0;
   const real   MinPres_No          = -1.0;
   const real   MinTemp_No          = -1.0;
   const real   MinEntr_No          = -1.0;

// 1. get sibling data
#  ifdef LOAD_BALANCE
   Buf_GetBufferData( lv, FluSg, MagSg, NULL_INT, DATA_GENERAL, _DENS, _NONE, Flu_ParaBuf, USELB_YES );
#  endif

// 1-0. reset the magnetic field
#  ifdef MHD
   if ( ResetMag ) {

// 1-1.
   if ( WLMDwarfGalaxy_UniformB ) {
// 1-2. allocate memory
   real (*Ax)[ CUBE(PS1+1) ] = NULL;
   real (*Ay)[ CUBE(PS1+1) ] = NULL;
   real (*Az)[ CUBE(PS1+1) ] = NULL;

   if ( UseVecPot )
   {
      Ax = new real [NT][ CUBE(PS1+1) ];
      Ay = new real [NT][ CUBE(PS1+1) ];
      Az = new real [NT][ CUBE(PS1+1) ];
   }


#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    1-3. compute vector potential
//         --> compute the entire patch at once to avoid redundant calculations
      real *AxTID = ( UseVecPot ) ? Ax[TID] : NULL;
      real *AyTID = ( UseVecPot ) ? Ay[TID] : NULL;
      real *AzTID = ( UseVecPot ) ? Az[TID] : NULL;

      if ( UseVecPot )
      {
         int idx = 0;

         for (int k=0; k<PS1+1; k++) {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + k*dh;
         for (int j=0; j<PS1+1; j++) {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + j*dh;
         for (int i=0; i<PS1+1; i++) {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + i*dh;

            AxTID[idx] = (real)0.0;
            AyTID[idx] = (real)0.0;
            AzTID[idx] = (real)0.0;

            if ( i != PS1 )   AxTID[idx] = MHD_ResetByUser_VecPot_Ptr( x+dh_2, y,      z,      TimeNew, dt, lv, 'x', NULL );
            if ( j != PS1 )   AyTID[idx] = MHD_ResetByUser_VecPot_Ptr( x,      y+dh_2, z,      TimeNew, dt, lv, 'y', NULL );
            if ( k != PS1 )   AzTID[idx] = MHD_ResetByUser_VecPot_Ptr( x,      y,      z+dh_2, TimeNew, dt, lv, 'z', NULL );

            idx ++;
         }}} // i,j,k
      } // if ( UseVecPot )


//    1-4. reset B field
//         --> set one component at a time since different components are defined at different cell faces
      for (int v=0; v<NCOMP_MAG; v++)
      {
         int    ijk_end[3], idx=0;
         double dxyz0[3];

         for (int d=0; d<3; d++)
         {
            ijk_end[d] = ( d == v ) ? PS1+1 : PS1;
            dxyz0  [d] = ( d == v ) ? 0.0   : dh_2;
         }

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + dxyz0[0];
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + dxyz0[1];
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + dxyz0[2];

         for (int k=0; k<ijk_end[2]; k++)    {  const double z = z0 + k*dh;
         for (int j=0; j<ijk_end[1]; j++)    {  const double y = y0 + j*dh;
         for (int i=0; i<ijk_end[0]; i++)    {  const double x = x0 + i*dh;

            real B_in, B_out;

            B_in  = amr->patch[MagSg][lv][PID]->magnetic[v][idx];
            B_out = MHD_ResetByUser_BField_Ptr( x, y, z, TimeNew, dt, lv, 'x'+v, NULL, B_in,
                                                UseVecPot, AxTID, AyTID, AzTID, i, j, k );

            amr->patch[MagSg][lv][PID]->magnetic[v][ idx ++ ] = B_out;
         }}} // i,j,k
      } // for (int v=0; v<NCOMP_MAG; v++)


//    1-5. update the total energy density
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         const real Emag_new = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
         amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i] += Emag_new;
      }
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// 1-6. free memory
   if ( UseVecPot )
   {
      delete [] Ax;
      delete [] Ay;
      delete [] Az;
   }

   } // ( WLMDwarfGalaxy_UniformB )
   else
   {
      MHD_Init_BField_ByVecPot_File_WLMDwarfGalaxy( lv );

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
            amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i] += Emag;
         }
      }
   }

   } // if ( ResetMag )
#  endif // #ifdef MHD

} // FUNCTION : Flu_ResetByUser_API_WLMDwarfGalaxy
