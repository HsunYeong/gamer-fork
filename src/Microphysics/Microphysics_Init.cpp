#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Microphysics_Init
// Description :  Initialize the microphysics routines
//
// Note        :  1. Invoked by Init_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Microphysics_Init()
{

// check if microphysics has been initialized already
   static bool MicroPhy_Initialized = false;

   if ( MicroPhy_Initialized )  return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   MicroPhy.useless            = NULL_BOOL;  // to avoid the empty structure issue
#  ifdef CR_DIFFUSION
   MicroPhy.CR_safety          = DT__CR_DIFFUSION;
   MicroPhy.CR_diff_coeff_para = CR_DIFF_PARA;
   MicroPhy.CR_diff_coeff_perp = CR_DIFF_PERP;
   MicroPhy.CR_diff_min_b      = CR_DIFF_MIN_B;
#  endif // #ifdef CR_DIFFUSION

#  ifdef CR_STREAMING
   MicroPhy.CR_source          = CR_SOURCE;
   MicroPhy.CR_stream          = CR_STREAM;
   MicroPhy.CR_Ec_source       = CR_EC_SOURCE;
   MicroPhy.CR_vmax            = CR_VMAX;
// CR_SIGMA/CR_SIGMA_PERP are the physical inverse-diffusion coefficients sigma' of Jiang & Oh (2018),
// defined so that the field-aligned diffusivity is D = 1/(3*sigma') (Eq. 25, Sec. 4.1.4).
// Internally the two-moment source/transport terms use the "code-unit" opacity sigma_code = Vm*sigma'
// (the same convention as Athena++; see src/pgen/cr_diffusion.cpp where D = vmax/(3*sigma)). This also
// matches the streaming opacity sigma_adv, which already carries a factor of Vm via invlim=1/Vm.
// --> multiply by CR_VMAX here so that D = 1/(3*CR_SIGMA) as intended.
   MicroPhy.CR_sigma           = CR_SIGMA      * CR_VMAX;
   MicroPhy.CR_sigma_perp      = CR_SIGMA_PERP * CR_VMAX;
   MicroPhy.CR_max_opacity     = CR_MAX_OPACITY;
   MicroPhy.CR_tau_asym_lim    = CR_TAU_ASYM_LIM;
   MicroPhy.CR_taufact         = CR_TAUFACT;
   MicroPhy.CR_vel_flx_flag    = CR_VEL_FLX_FLAG;
   MicroPhy.CR_cfl             = CR_CFL;
#  endif // #ifdef CR_STREAMING

   MicroPhy_Initialized = true;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Microphysics_Init
