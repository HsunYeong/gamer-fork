#include "GAMER.h"



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

// 1-1. reset the magnetic field
#  ifdef MHD
   if ( ResetMag ) {

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

   real (*DensIn)[ CUBE(PS1+2) ] = new real [Der_NP][ CUBE(PS1+2) ];

   typedef real (*Dens3D)[PS1+2][PS1+2][PS1+2];

// loop over patch groups
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {

//    get patch data with one ghost size
      Prepare_PatchData( lv, Time[lv], DensIn[0], NULL, 1, 1, &PID0,
                         _DENS, _NONE, OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26,
                         IntPhase_No, OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No,
                         DE_Consistency_No );

      Dens3D DensIn3D = ( Dens3D )DensIn;

#     pragma omp parallel for schedule( runtime )
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID  = PID0 + LocalID;

#        ifdef OPENMP
         const int TID = omp_get_thread_num();
#        else
         const int TID = 0;
#        endif

   //    1-3. compute vector potential
   //         --> compute the entire patch at once to avoid redundant calculations
         real *AxTID = ( UseVecPot ) ? Ax[TID] : NULL;
         real *AyTID = ( UseVecPot ) ? Ay[TID] : NULL;
         real *AzTID = ( UseVecPot ) ? Az[TID] : NULL;

         if ( UseVecPot )
         {
            int idx = 0;

            for (int k=0; k<PS1+1; k++) {  const double z  = amr->patch[0][lv][PID]->EdgeL[2] + k*dh;
                                           const int    ki = k + 1;
            for (int j=0; j<PS1+1; j++) {  const double y  = amr->patch[0][lv][PID]->EdgeL[1] + j*dh;
                                           const int    ji = j + 1;
            for (int i=0; i<PS1+1; i++) {  const double x  = amr->patch[0][lv][PID]->EdgeL[0] + i*dh;
                                           const int    ii = i + 1;

               AxTID[idx] = (real)0.0;
               AyTID[idx] = (real)0.0;
               AzTID[idx] = (real)0.0;

               double DensArray[3];
//             compute edge-centered density
               // dens( x+dh_2, y,      z     )
               DensArray[0] = 0.25*(DensIn3D[LocalID][ki][ji][ii] + DensIn3D[LocalID][k ][ji][ii] + DensIn3D[LocalID][ki][j ][ii] + DensIn3D[LocalID][k ][j ][ii]);
               // dens( x,      y+dh_2, z     )
               DensArray[1] = 0.25*(DensIn3D[LocalID][ki][ji][ii] + DensIn3D[LocalID][ki][ji][i ] + DensIn3D[LocalID][k ][ji][ii] + DensIn3D[LocalID][k ][ji][i ]);
               // dens( x,      y,      z+dh_2)
               DensArray[2] = 0.25*(DensIn3D[LocalID][ki][ji][ii] + DensIn3D[LocalID][ki][j ][ii] + DensIn3D[LocalID][ki][ji][i ] + DensIn3D[LocalID][ki][j ][i ]);

               if ( i != PS1 )   AxTID[idx] = MHD_ResetByUser_VecPot_Ptr( x+dh_2, y,      z,      TimeNew, dt, lv, 'x', DensArray );
               if ( j != PS1 )   AyTID[idx] = MHD_ResetByUser_VecPot_Ptr( x,      y+dh_2, z,      TimeNew, dt, lv, 'y', DensArray );
               if ( k != PS1 )   AzTID[idx] = MHD_ResetByUser_VecPot_Ptr( x,      y,      z+dh_2, TimeNew, dt, lv, 'z', DensArray );

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
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


// 1-6. free memory
   if ( UseVecPot )
   {
      delete [] Ax;
      delete [] Ay;
      delete [] Az;
   }
   delete [] DensIn;

   } // if ( ResetMag )
#  endif // #ifdef MHD

} // FUNCTION : Flu_ResetByUser_API_WLMDwarfGalaxy
