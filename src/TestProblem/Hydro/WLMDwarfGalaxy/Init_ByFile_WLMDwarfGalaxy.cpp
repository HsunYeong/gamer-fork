#include "GAMER.h"


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFile_WLMDwarfGalaxy
// Description :  Function to actually set the fluid field from the input uniform-mesh array
//
// Note        :  1. Invoked by Init_ByFile_AssignData() using the function pointer Init_ByFile_User_Ptr()
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Does not floor and normalize passive scalars
//                3. Calculate the dual-energy variable automatically instead of load it from the disk
//                   --> When adopting DUAL_ENERGY, the input uniform-mesh array must NOT include the dual-energy
//                       variable
//                4. ELBDM:
//                   ELBDM_SCHEME == ELBDM_WAVE:
//                   --> Calculate the density field automatically instead of loading it from the disk for ELBDM
//                   --> The input uniform-mesh array must NOT include the density field
//                   ELBDM_SCHEME == ELBDM_HYBRID:
//                   --> We will load the density and phase fields from the disk on all levels
//                   --> There is no need to separately calculate the density field
//                5. Assuming nvar_in (i.e., OPT__UM_IC_NVAR) == NCOMP_TOTAL
//                   --> Unless either DUAL_ENERGY or ELBDM is adopted, for which it assumes nvar_in == NCOMP_TOTAL-1
//
// Parameter   :  fluid_out : Fluid field to be set
//                fluid_in  : Fluid field loaded from the uniform-mesh array (UM_IC)
//                nvar_in   : Number of variables in fluid_in
//                x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                lv        : Target AMR level
//                AuxArray  : Auxiliary array
//
// Return      :  fluid_out
//-------------------------------------------------------------------------------------------------------
void Init_ByFile_WLMDwarfGalaxy( real fluid_out[], const real fluid_in[], const int nvar_in,
                                 const double x, const double y, const double z, const double Time,
                                 const int lv, double AuxArray[] )
{

#  ifdef GAMER_DEBUG
#  if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )
   if ( nvar_in != NCOMP_TOTAL-1 )
      Aux_Error( ERROR_INFO, "nvar_in (%d) != NCOMP_TOTAL-1 (%d) when enabling DUAL_ENERGY !!\n", nvar_in, NCOMP_TOTAL-1 );

#  elif ( MODEL == ELBDM )
   if ( nvar_in != NCOMP_TOTAL-1 )
      Aux_Error( ERROR_INFO, "nvar_in (%d) != NCOMP_TOTAL-1 (%d) for ELBDM !!\n", nvar_in, NCOMP_TOTAL-1 );

#  else
   if ( nvar_in != NCOMP_TOTAL )
      Aux_Error( ERROR_INFO, "nvar_in (%d) != NCOMP_TOTAL (%d) !!\n", nvar_in, NCOMP_TOTAL );
#  endif
#  endif // #ifdef GAMER_DEBUG

   for (int v_in=0, v_out=0; v_in<nvar_in; v_in++, v_out++)
   {
//    skip the dual-energy field for HYDRO
#     if   ( MODEL == HYDRO )
#     ifdef DUAL_ENERGY
      if ( v_out == DUAL )    v_out ++;
#     endif
#     endif // MODEL

      fluid_out[v_out] = fluid_in[v_in];
   }

// calculate the dual-energy field for HYDRO
#  if   ( MODEL == HYDRO )

#  ifdef MHD
   const real Emag = 0.0;
#  else
   const real Emag = NULL_REAL;
#  endif

#  ifdef DUAL_ENERGY
   fluid_out[DUAL] = Hydro_Con2Dual( fluid_in[DENS], fluid_in[MOMX], fluid_in[MOMY], fluid_in[MOMZ], fluid_in[ENGY], Emag,
                                     EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                     PassiveFloorMask );
#  endif

#  endif // MODEL

} // FUNCTION : Init_ByFile_WLMDwarfGalaxy



