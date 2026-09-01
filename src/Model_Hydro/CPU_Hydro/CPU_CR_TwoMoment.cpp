#include "CUFLU.h"

#ifdef CR_STREAMING




//-------------------------------------------------------------------------------------------------------
// Function    : RotateVec
// Description : Rotate a vector from lab frame to B-aligned frame
//
// Note        : 1. B-aligned frame: x' along B, y' and z' perpendicular to B
//               2. The rotation is defined by two angles: theta (B to z-axis) and phi (Bxy to x-axis)
//
// Parameter   : sint : sin(theta) = |Bxy|/|B|
//               cost : cos(theta) =    Bz/|B|
//               sinp : sin(phi  ) =    By/|Bxy|
//               cosp : cos(phi  ) =    Bx/|Bxy|
//               v1, v2, v3 : input/output vector components
//
// Reference   : Athena++ src/utils/rotate_vectors.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void RotateVec( const real sint, const real cost,
                       const real sinp, const real cosp,
                       real &v1, real &v2, real &v3 )
{
// The two rotation matrices:
// R_1 (around z by phi):
// [cos_p  sin_p 0]
// [-sin_p cos_p 0]
// [0       0    1]
//
// R_2 (around y by theta):
// [sin_t  0 cos_t]
// [0      1     0]
// [-cos_t 0 sin_t]

// First apply R1, then apply R2
   real newv1 =  cosp * v1 + sinp * v2;
   v2         = -sinp * v1 + cosp * v2;

// now apply R2
   v1         =  sint * newv1 + cost * v3;
   real newv3 = -cost * newv1 + sint * v3;
   v3 = newv3;

} // FUNCTION : RotateVec



//-------------------------------------------------------------------------------------------------------
// Function    : InvRotateVec
// Description : Rotate a vector from B-aligned frame back to lab frame (inverse of RotateVec)
//
// Note        : 1. Applies R1^-1 * R2^-1
//
// Parameter   : sint, cost, sinp, cosp : same as RotateVec
//               v1, v2, v3 : input/output vector components
//
// Reference   : Athena++ src/utils/rotate_vectors.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void InvRotateVec( const real sint, const real cost,
                          const real sinp, const real cosp,
                          real &v1, real &v2, real &v3 )
{
// R_1^-1 (around z by -phi):
// [cos_p  -sin_p 0]
// [sin_p  cos_p  0]
// [0       0     1]
//
// R_2^-1 (around y by -theta):
// [sin_t  0 -cos_t]
// [0      1    0  ]
// [cos_t  0  sin_t]

// First apply R2^-1, then apply R1^-1
   real newv1 = sint * v1 - cost * v3;
   v3         = cost * v1 + sint * v3;

// now apply R1^-1
   v1         = cosp * newv1 - sinp * v2;
   real newv2 = sinp * newv1 + cosp * v2;
   v2 = newv2;

} // FUNCTION : InvRotateVec



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeBFieldAngles
// Description : Compute the rotation angles from B field direction
//
// Note        : 1. Returns sin(theta), cos(theta), sin(phi), cos(phi)
//               2. theta is the angle between B and z-axis
//               3. phi is the angle between B_xy projection and x-axis
//
// Parameter   : Bx, By, Bz : magnetic field components
//               sint, cost, sinp, cosp : output angles
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_ComputeBFieldAngles( const real Bx, const real By, const real Bz,
                                    real &sint, real &cost, real &sinp, real &cosp )
{
   const real bxby = SQRT( SQR(Bx) + SQR(By) );
   const real btot = SQRT( SQR(Bx) + SQR(By) + SQR(Bz) );

   if ( btot > TINY_NUMBER ) {
      sint = bxby / btot;       // sin(theta) = |Bxy|/|B|
      cost = Bz / btot;         // cos(theta) = Bz/|B|
   } else {
      sint = (real)1.0;
      cost = (real)0.0;
   }

   if ( bxby > TINY_NUMBER ) {
      sinp = By / bxby;         // sin(phi) = By/|Bxy|
      cosp = Bx / bxby;         // cos(phi) = Bx/|Bxy|
   } else {
      sinp = (real)0.0;
      cosp = (real)1.0;
   }

} // FUNCTION : CR_ComputeBFieldAngles



//-------------------------------------------------------------------------------------------------------
// Function    : CR_UpdateStreaming_OneCell
// Description : Update streaming opacity (sigma_adv) and streaming velocity (v_adv) from grad(Pc)
//               for a single cell
//
// Note        : 1. Based on Athena++ DefaultStreaming() in cr.cpp
//               2. sigma_adv[0] is parallel to B, sigma_adv[1,2] are perpendicular (set to max_opacity)
//               3. v_adv = -sign(B dot grad Pc) * v_Alfven
//               4. This function should be called after flux calculation with grad_pc computed
//                  from flux divergence
//
// Parameter   : Ec          : CR energy density
//               rho         : gas density
//               Bx, By, Bz  : magnetic field components
//               grad_pc[3]  : gradient of CR pressure (dPc/dx, dPc/dy, dPc/dz)
//               vmax        : effective speed of light
//               sigma_adv   : output streaming opacity (parallel component)
//               v_adv[3]    : output streaming velocity components
//               MicroPhy    : Microphysics object containing CR_max_opacity
//
// Reference   : Athena++ src/cr/cr.cpp DefaultStreaming()
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_UpdateStreaming_OneCell( const real Ec, const real rho,
                              const real Bx, const real By, const real Bz,
                              const real grad_pc[3], const real vmax,
                              real &sigma_adv, real v_adv[3], const MicroPhy_t *MicroPhy )
{
   const real invlim = (real)1.0 / vmax;

// compute B field magnitude and Alfven velocity
   const real bsq  = SQR(Bx) + SQR(By) + SQR(Bz);
   const real btot = SQRT( bsq );
   const real inv_sqrt_rho = (real)1.0 / SQRT( rho );
   const real va = btot * inv_sqrt_rho;

// compute B dot grad(Pc)
   const real b_grad_pc = Bx * grad_pc[0] + By * grad_pc[1] + Bz * grad_pc[2];

// determine sign of B dot grad(Pc)
   real dpc_sign = (real)0.0;
   if ( b_grad_pc > TINY_NUMBER )
      dpc_sign = (real)1.0;
   else if ( -b_grad_pc > TINY_NUMBER )
      dpc_sign = (real)-1.0;

// compute streaming velocity: v_adv = -sign(B dot grad Pc) * v_Alfven * b_hat
   const real va1 = Bx * inv_sqrt_rho;
   const real va2 = By * inv_sqrt_rho;
   const real va3 = Bz * inv_sqrt_rho;

   v_adv[0] = -va1 * dpc_sign;
   v_adv[1] = -va2 * dpc_sign;
   v_adv[2] = -va3 * dpc_sign;

// compute streaming opacity (parallel to B)
// sigma_adv = |B dot grad Pc| / (|B| * v_A * (4/3) * (1/vmax) * Ec)
   if ( va > TINY_NUMBER && Ec > TINY_NUMBER ) {
      sigma_adv = FABS(b_grad_pc) / ( btot * va * ((real)4.0/(real)3.0) * invlim * Ec );
   } else {
      sigma_adv = MicroPhy->CR_max_opacity;
   }

} // FUNCTION : CR_UpdateStreaming_OneCell



//-------------------------------------------------------------------------------------------------------
// Function    : CR_UpdateStreaming
// Description : Update streaming opacity (sigma_adv) and streaming velocity (v_adv) from flux divergence
//               for all interior cells
//
// Note        : 1. This function loops over interior cells and calls CR_UpdateStreaming_OneCell for each
//               2. Works for both half-step and full-step by using appropriate parameters
//               3. grad_pc[n] is the FULL flux divergence of the CR_F(n+1) equation, summing all
//                  three flux directions (9 terms total): grad_pc[n] = sum_d (F[d][CR_F(n+1),d+1/2]
//                  - F[d][CR_F(n+1),d-1/2]) / dx / vmax, matching Athena++'s grad_pc_
//               4. Based on Athena++ DefaultStreaming() - called AFTER flux calculation
//
// Parameter   : g_Output     : Array to store the updated opacity (ADV_SIGMA, ADV_VX, ADV_VY, ADV_VZ)
//               g_CellVar    : Array storing cell-centered variables (for reading rho, Ec)
//                              If NULL, read from g_Output
//               g_CC_B       : Array storing cell-centered B field [3][ CUBE(N) ]
//               g_Flux       : Array storing fluxes [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ]
//               NFlux        : Stride for accessing g_Flux[] (N_HF_FLUX or N_FL_FLUX)
//               NVar_Out     : Stride for accessing g_Output[] (FLU_NXT, N_HF_VAR, or PS2)
//               NVar_In      : Stride for accessing g_CellVar[] (FLU_NXT or N_HF_VAR)
//               NVar_B       : Stride for accessing g_CC_B[] (FLU_NXT or N_HF_VAR)
//               out_offset   : Offset between flux cell index and output index
//               in_offset    : Offset between flux cell index and input (g_CellVar) index
//               dh           : Cell size
//               MicroPhy     : Microphysics object containing CR_vmax and CR_max_opacity
//
// Return      : Modified ADV_SIGMA, ADV_VX, ADV_VY, ADV_VZ in g_Output
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_UpdateStreaming( real g_Output[][ CUBE(FLU_NXT) ],
                       const real g_CellVar[][ CUBE(FLU_NXT) ],
                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                       const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                       const int NFlux, const int NVar_Out, const int NVar_In, const int NVar_B,
                       const int out_offset, const int in_offset,
                       const real dh, const MicroPhy_t *MicroPhy )
{
   const real vmax = MicroPhy->CR_vmax;
   const real _dh  = (real)1.0 / dh;
   const int  cell_offset = 1;  // skip boundary cells
   const int  didx_flux[3] = { 1, NFlux, SQR(NFlux) };
   const int  CRF_v[3] = { CR_F1, CR_F2, CR_F3 };

// determine input array: use g_CellVar if provided, else g_Output
   const real (*g_Input)[CUBE(FLU_NXT)] = ( g_CellVar != NULL ) ? g_CellVar : g_Output;

// loop bounds
   const int cell_size_i  = NFlux - 2*cell_offset;
   const int cell_size_j  = NFlux - 2*cell_offset;
   const int cell_size_k  = NFlux - 2*cell_offset;
   const int cell_size_ij = cell_size_i * cell_size_j;

   CGPU_LOOP( idx_cell, cell_size_i*cell_size_j*cell_size_k )
   {
//    flux cell indices
      const int i_cell = idx_cell % cell_size_i + cell_offset;
      const int j_cell = idx_cell % cell_size_ij / cell_size_i + cell_offset;
      const int k_cell = idx_cell / cell_size_ij + cell_offset;
      const int idx_flux = IDX321( i_cell, j_cell, k_cell, NFlux, NFlux );

//    output indices (with out_offset)
      const int i_out = i_cell + out_offset;
      const int j_out = j_cell + out_offset;
      const int k_out = k_cell + out_offset;
      const int idx_out = IDX321( i_out, j_out, k_out, NVar_Out, NVar_Out );

//    input indices (with in_offset)
      const int i_in = i_cell + in_offset;
      const int j_in = j_cell + in_offset;
      const int k_in = k_cell + in_offset;
      const int idx_in = IDX321( i_in, j_in, k_in, NVar_In, NVar_In );

//    B field index (may differ from input index if NVar_B != NVar_In)
      const int idx_B = IDX321( i_in, j_in, k_in, NVar_B, NVar_B );

//    compute grad_pc from flux divergence
//    --> matching Athena++ (cr_transport.cpp:366-403), each component n is the FULL divergence of
//        the flux vector of the CR_F(n+1) equation, summing all three flux directions; the
//        transverse fluxes carry nonzero HLLE upwind/dissipation terms (-bm*Fc_L, -bp*Fc_R) even
//        though the off-diagonal Eddington factors vanish, so they must not be dropped
      real grad_pc[3];
      for (int n=0; n<3; n++) {
         grad_pc[n] = (real)0.0;
         for (int d=0; d<3; d++)
            grad_pc[n] += g_Flux[d][ CRF_v[n] ][idx_flux] - g_Flux[d][ CRF_v[n] ][ idx_flux - didx_flux[d] ];
         grad_pc[n] *= _dh / vmax;
      }

//    get cell-centered values for updating opacity
      const real Ec  = g_Input[CR_E ][idx_in];
      const real rho = g_Input[DENS][idx_in];

//    get B field from g_CC_B array
      const real Bx = g_CC_B[0][idx_B];
      const real By = g_CC_B[1][idx_B];
      const real Bz = g_CC_B[2][idx_B];

//    call CR_UpdateStreaming_OneCell to compute new sigma_adv and v_adv
      real sigma_adv_new;
      real v_adv_new[3];
      CR_UpdateStreaming_OneCell( Ec, rho, Bx, By, Bz, grad_pc, vmax, sigma_adv_new, v_adv_new, MicroPhy );

//    store updated values
      g_Output[ADV_SIGMA][idx_out] = sigma_adv_new;
      g_Output[ADV_VX   ][idx_out] = v_adv_new[0];
      g_Output[ADV_VY   ][idx_out] = v_adv_new[1];
      g_Output[ADV_VZ   ][idx_out] = v_adv_new[2];

   } // CGPU_LOOP

} // FUNCTION : CR_UpdateStreaming



//-------------------------------------------------------------------------------------------------------
// Function    : CR_UpdateOpacity
// Description : Update streaming opacity (sigma_adv) and streaming velocity (v_adv) using
//               CENTRAL DIFFERENCE gradient of CR pressure
//
// Note        : 1. Based on Athena++ DefaultOpacity() in cr.cpp
//               2. This function uses central differences to compute grad(Pc), NOT flux divergence
//               3. Should be called AFTER the source term (like Athena's DefaultOpacity is called
//                  after AddSourceTerms in the task list)
//               4. The sigma_adv and v_adv computed here will be used in the NEXT timestep's
//                  flux calculation
//               5. g_Output may hold either conserved data (g_Flu_Array_In) or primitive data
//                  (g_PriVar_Half, i.e. AFTER Con2Pri in Hydro_RiemannPredict): safe because only
//                  DENS, CR_E, and ADV_* are accessed, which are identical in the two
//                  representations (CR passive fields are not converted by Hydro_Con2Pri())
//
// Parameter   : g_Output     : Flat pointer to the output array for updated opacity
//               OutStride    : Stride between variables in g_Output (CUBE(FLU_NXT) or CUBE(PS2))
//               g_CC_B       : Array storing cell-centered B field [3][ CUBE(N) ]
//               NVar_Out     : Size for computing output indices (FLU_NXT, N_HF_VAR, or PS2)
//               NVar_In      : Size for computing input indices for neighbor access
//               NVar_B       : Stride for accessing g_CC_B[] (FLU_NXT or N_HF_VAR)
//               out_offset   : Offset for output index relative to interior
//               in_offset    : Offset for input index relative to interior
//               NSize        : Size of the interior region to update
//               dh           : Cell size
//               MicroPhy     : Microphysics object containing vmax
//
// Return      : Modified ADV_SIGMA, ADV_VX, ADV_VY, ADV_VZ in g_Output
//
// Reference   : Athena++ src/cr/cr.cpp DefaultOpacity()
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_UpdateOpacity( real *g_Output,
                       const int OutStride,
                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                       const int NVar_Out, const int NVar_In, const int NVar_B,
                       const int out_offset, const int in_offset, const int NSize,
                       const real dh, const MicroPhy_t *MicroPhy )
{
   const real vmax   = MicroPhy->CR_vmax;
   const real invlim = (real)1.0 / vmax;
   const real _2dh   = (real)0.5 / dh;

// loop bounds: NSize already accounts for boundary cells (caller provides NSize-2)
   const int size_ij = NSize * NSize;

   CGPU_LOOP( idx_cell, NSize*NSize*NSize )
   {
//    cell indices within the loop region
      const int i_cell = idx_cell % NSize;
      const int j_cell = idx_cell % size_ij / NSize;
      const int k_cell = idx_cell / size_ij;

//    output indices (with out_offset, caller adds +1 for boundary skip)
      const int i_out = i_cell + out_offset;
      const int j_out = j_cell + out_offset;
      const int k_out = k_cell + out_offset;
      const int idx_out = IDX321( i_out, j_out, k_out, NVar_Out, NVar_Out );

//    input indices for central difference (with in_offset)
      const int i_in = i_cell + in_offset;
      const int j_in = j_cell + in_offset;
      const int k_in = k_cell + in_offset;
      const int idx_in = IDX321( i_in, j_in, k_in, NVar_In, NVar_In );

//    neighbor indices for central difference
      const int idx_ip1 = IDX321( i_in+1, j_in,   k_in,   NVar_In, NVar_In );
      const int idx_im1 = IDX321( i_in-1, j_in,   k_in,   NVar_In, NVar_In );
      const int idx_jp1 = IDX321( i_in,   j_in+1, k_in,   NVar_In, NVar_In );
      const int idx_jm1 = IDX321( i_in,   j_in-1, k_in,   NVar_In, NVar_In );
      const int idx_kp1 = IDX321( i_in,   j_in,   k_in+1, NVar_In, NVar_In );
      const int idx_km1 = IDX321( i_in,   j_in,   k_in-1, NVar_In, NVar_In );

//    B field index
      const int idx_B = IDX321( i_in, j_in, k_in, NVar_B, NVar_B );

//    compute grad(Pc) using central differences
//    Pc = Ec / 3, so grad(Pc) = (1/3) * grad(Ec)
//    Following Athena++: dprdx = (Ec[i+1] - Ec[i-1]) / 3.0 / (2*dx)
      real grad_pc[3];
      grad_pc[0] = ( g_Output[CR_E*OutStride + idx_ip1] - g_Output[CR_E*OutStride + idx_im1] ) / (real)3.0 * _2dh;
      grad_pc[1] = ( g_Output[CR_E*OutStride + idx_jp1] - g_Output[CR_E*OutStride + idx_jm1] ) / (real)3.0 * _2dh;
      grad_pc[2] = ( g_Output[CR_E*OutStride + idx_kp1] - g_Output[CR_E*OutStride + idx_km1] ) / (real)3.0 * _2dh;

//    get cell-centered values
      const real Ec  = g_Output[CR_E*OutStride + idx_in];
      const real rho = g_Output[DENS*OutStride + idx_in];

//    get B field from g_CC_B array
      const real Bx = g_CC_B[0][idx_B];
      const real By = g_CC_B[1][idx_B];
      const real Bz = g_CC_B[2][idx_B];

//    compute B field magnitude and Alfven velocity
      const real bsq  = SQR(Bx) + SQR(By) + SQR(Bz);
      const real btot = SQRT( bsq );
      const real inv_sqrt_rho = (real)1.0 / SQRT( rho );
      const real va = btot * inv_sqrt_rho;

//    compute B dot grad(Pc)
      const real b_grad_pc = Bx * grad_pc[0] + By * grad_pc[1] + Bz * grad_pc[2];

//    determine sign of B dot grad(Pc)
      real dpc_sign = (real)0.0;
      if ( b_grad_pc > TINY_NUMBER )
         dpc_sign = (real)1.0;
      else if ( -b_grad_pc > TINY_NUMBER )
         dpc_sign = (real)-1.0;

//    compute streaming velocity: v_adv = -sign(B dot grad Pc) * v_Alfven * b_hat
      const real va1 = Bx * inv_sqrt_rho;
      const real va2 = By * inv_sqrt_rho;
      const real va3 = Bz * inv_sqrt_rho;

      real v_adv[3];
      v_adv[0] = -va1 * dpc_sign;
      v_adv[1] = -va2 * dpc_sign;
      v_adv[2] = -va3 * dpc_sign;

//    compute streaming opacity (parallel to B)
//    sigma_adv = |B dot grad Pc| / (|B| * v_A * (4/3) * (1/vmax) * Ec)
      real sigma_adv;
      if ( va > TINY_NUMBER && Ec > TINY_NUMBER ) {
         sigma_adv = FABS(b_grad_pc) / ( btot * va * ((real)4.0/(real)3.0) * invlim * Ec );
      } else {
         sigma_adv = MicroPhy->CR_max_opacity;
      }

//    store updated values
      g_Output[ADV_SIGMA*OutStride + idx_out] = sigma_adv;
      g_Output[ADV_VX   *OutStride + idx_out] = v_adv[0];
      g_Output[ADV_VY   *OutStride + idx_out] = v_adv[1];
      g_Output[ADV_VZ   *OutStride + idx_out] = v_adv[2];

   } // CGPU_LOOP

} // FUNCTION : CR_UpdateOpacity



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeVdiff
// Description : Compute v_diff (diffusion velocity) at a CELL CENTER along a given direction
//
// Note        : 1. Based on Athena++ cr_transport.cpp CalculateFluxes()
//               2. Computes v_diff in B-aligned frame, then rotates to lab frame
//               3. Returns v_diff component along the specified flux direction
//               4. This should be called separately for left and right cells at each interface
//
// Parameter   : sigma_adv   : streaming opacity (parallel to B)
//               Bx, By, Bz  : cell-centered magnetic field components
//               Ec          : CR energy density
//               rho         : gas density
//               vmax        : effective speed of light
//               dh          : cell size
//               fdir        : flux direction (0=x, 1=y, 2=z)
//               MicroPhy    : Microphysics object containing CR_sigma, CR_sigma_perp, CR_max_opacity,
//                             CR_taufact, CR_tau_asym_lim, and CR_vel_flx_flag
//
// Return      : v_diff along the specified direction (fdir)
//
// Reference   : Athena++ src/cr/integrators/cr_transport.cpp lines 77-170
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static real CR_ComputeVdiff( const real sigma_adv,
                             const real Bx, const real By, const real Bz,
                             const real Ec, const real rho,
                             const real vmax, const real dh, const int fdir,
                             const MicroPhy_t *MicroPhy )
{
   const real edd = (real)1.0 / (real)3.0;
   const real sigma_adv_perp = MicroPhy->CR_max_opacity;
   const real sigma_diff = MicroPhy->CR_sigma;
   const real sigma_diff_perp = MicroPhy->CR_sigma_perp;

// compute B-field angles for rotation
   real sint, cost, sinp, cosp;
   CR_ComputeBFieldAngles( Bx, By, Bz, sint, cost, sinp, cosp );

// total sigma: combine sigma_diff with the streaming opacity sigma_adv only if streaming is enabled
// In B-aligned frame: sigma_x is parallel to B, sigma_y and sigma_z are perpendicular
   real sigma_x = sigma_diff;
   real sigma_y = sigma_diff_perp;
   real sigma_z = sigma_diff_perp;
   if ( MicroPhy->CR_stream ) {
      sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv );
      sigma_y = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
      sigma_z = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
   }

// compute tau and diffv for each B-aligned direction
// x direction (parallel to B)
   real tau_x = MicroPhy->CR_taufact * sigma_x * dh;
   tau_x = tau_x * tau_x / ( (real)2.0 * edd );
   real diffv_x;
   if ( tau_x < MicroPhy->CR_tau_asym_lim )
      diffv_x = SQRT( (real)1.0 - (real)0.5 * tau_x );
   else
      diffv_x = SQRT( ( (real)1.0 - EXP(-tau_x) ) / tau_x );

// y direction (perpendicular to B)
   real tau_y = MicroPhy->CR_taufact * sigma_y * dh;
   tau_y = tau_y * tau_y / ( (real)2.0 * edd );
   real diffv_y;
   if ( tau_y < MicroPhy->CR_tau_asym_lim )
      diffv_y = SQRT( (real)1.0 - (real)0.5 * tau_y );
   else
      diffv_y = SQRT( ( (real)1.0 - EXP(-tau_y) ) / tau_y );

// z direction (perpendicular to B)
   real tau_z = MicroPhy->CR_taufact * sigma_z * dh;
   tau_z = tau_z * tau_z / ( (real)2.0 * edd );
   real diffv_z;
   if ( tau_z < MicroPhy->CR_tau_asym_lim )
      diffv_z = SQRT( (real)1.0 - (real)0.5 * tau_z );
   else
      diffv_z = SQRT( ( (real)1.0 - EXP(-tau_z) ) / tau_z );

// v_diff in B-aligned frame (component 1 = parallel to B, components 2/3 = perpendicular)
   real vdiff_1 = vmax * SQRT(edd) * diffv_x;
   real vdiff_2 = vmax * SQRT(edd) * diffv_y;
   real vdiff_3 = vmax * SQRT(edd) * diffv_z;

// rotate from B-aligned frame to lab frame
// --> vdiff_1/2/3 hold the lab-frame x/y/z components from here on
   InvRotateVec( sint, cost, sinp, cosp, vdiff_1, vdiff_2, vdiff_3 );

// take absolute value
   vdiff_1 = FABS( vdiff_1 );
   vdiff_2 = FABS( vdiff_2 );
   vdiff_3 = FABS( vdiff_3 );

// add CR sound speed for stability
   const real cr_sound = MicroPhy->CR_vel_flx_flag * SQRT( ((real)4.0/(real)9.0) * Ec / rho );
   vdiff_1 += cr_sound;
   vdiff_2 += cr_sound;
   vdiff_3 += cr_sound;

// return component along flux direction
   if ( fdir == 0 )      return vdiff_1;
   else if ( fdir == 1 ) return vdiff_2;
   else                  return vdiff_3;

} // FUNCTION : CR_ComputeVdiff



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeHLLEFlux
// Description : Compute HLLE flux for CR variables at a cell interface
//
// Note        : 1. Based on Athena++ cr_flux.cpp
//               2. fdir is the flux direction (0=x, 1=y, 2=z)
//
// Parameter   : Ec_L/R      : CR energy density (left/right)
//               Fc_L/R[3]   : CR flux (left/right)
//               v_L/R       : gas velocity along flux direction (left/right)
//               vdiff_L/R   : diffusion velocity (left/right)
//               vmax        : effective speed of light
//               fdir        : flux direction (0=CRF1, 1=CRF2, 2=CRF3)
//               flux_E      : output flux for CR energy
//               flux_F[3]   : output flux for CR flux components
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_ComputeHLLEFlux( const real Ec_L, const real Ec_R,
                                const real Fc_L[3], const real Fc_R[3],
                                const real v_L, const real v_R,
                                const real vdiff_L, const real vdiff_R,
                                const real vmax,
                                const int fdir,
                                real &flux_E, real flux_F[3] )
{
// Eddington factors (isotropic)
   const real edd   = (real)1.0 / (real)3.0;   // diagonal component
   const real edd12 = (real)0.0;               // off-diagonal components
   const real edd13 = (real)0.0;
   const real edd23 = (real)0.0;

// compute wave speeds
   const real meanadv  = (real)0.5 * ( v_L + v_R );
   const real meandiffv = (real)0.5 * ( vdiff_L + vdiff_R );

   real al = FMIN( meanadv - meandiffv, v_L - vdiff_L );
   real ar = FMAX( meanadv + meandiffv, v_R + vdiff_R );

// clamp to light speed
   const real sqrt_edd = SQRT( edd );
   ar = FMIN( ar,  vmax * sqrt_edd );
   al = FMAX( al, -vmax * sqrt_edd );

   const real bp = FMAX( ar, (real)0.0 );   // positive wave speed
   const real bm = FMIN( al, (real)0.0 );   // negative wave speed

// compute L/R fluxes minus wave transport: F_L - bm*U_L, F_R - bp*U_R
   real fl_e, fr_e;
   real fl_f1, fr_f1, fl_f2, fr_f2, fl_f3, fr_f3;

// energy equation flux: vmax * F_c - b * E_c
   fl_e = vmax * Fc_L[fdir] - bm * Ec_L;
   fr_e = vmax * Fc_R[fdir] - bp * Ec_R;

// flux equation fluxes depend on direction
   if ( fdir == 0 ) {       // x-direction (CRF1)
      fl_f1 = vmax * edd   * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd   * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd12 * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd12 * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd13 * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd13 * Ec_R - bp * Fc_R[2];
   }
   else if ( fdir == 1 ) {  // y-direction (CRF2)
      fl_f1 = vmax * edd12 * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd12 * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd   * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd   * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd23 * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd23 * Ec_R - bp * Fc_R[2];
   }
   else {                   // z-direction (CRF3)
      fl_f1 = vmax * edd13 * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd13 * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd23 * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd23 * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd   * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd   * Ec_R - bp * Fc_R[2];
   }

// HLLE flux formula
   real tmp = (real)0.0;
   if ( FABS(bm - bp) > TINY_NUMBER )
      tmp = (real)0.5 * (bp + bm) / (bp - bm);

   flux_E    = (real)0.5 * (fl_e  + fr_e ) + (fl_e  - fr_e ) * tmp;
   flux_F[0] = (real)0.5 * (fl_f1 + fr_f1) + (fl_f1 - fr_f1) * tmp;
   flux_F[1] = (real)0.5 * (fl_f2 + fr_f2) + (fl_f2 - fr_f2) * tmp;
   flux_F[2] = (real)0.5 * (fl_f3 + fr_f3) + (fl_f3 - fr_f3) * tmp;

} // FUNCTION : CR_ComputeHLLEFlux



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentFlux_HalfStep
//
// Description : Compute the half-step CR two-moment HLLE fluxes using donor-cell (first-order) states
//
// Note        : 1. Donor-cell reconstruction matches Athena++'s forced first-order reconstruction
//                  at VL2 stage 1
//               2. vdiff is evaluated at cell centers with the stored ADV_SIGMA
//               3. Overwrites the CR slots of g_Flux_Half[] computed by Hydro_RiemannPredict_Flux()
//                  and zeroes the ADV_* slots (sigma_adv/v_adv are work arrays recomputed each step,
//                  never advected --> matches Athena++, where they are not part of u_cr)
//
// Reference   : Athena++ src/cr/integrators/cr_transport.cpp, cr_flux.cpp
//
// Parameter   : g_Con_Var   : Array storing the input cell-centered conserved fluid variables
//               g_Flux_Half : Array with hydrodynamic fluxes for adding the cosmic-ray diffusive fluxes
//               g_FC_B      : Array storing the input face-centered B field
//               g_CC_B      : Array storing the input cell-centered B field
//               dh          : Cell size
//               MicroPhy    : Microphysics object
//
// Return      : g_Flux_Half[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentFlux_HalfStep( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                  real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],   //unuse
                            const real g_CC_B[][ CUBE(FLU_NXT) ],
                            const real dh, const MicroPhy_t *MicroPhy )
{


   const int  didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const int  flux_offset  = 1;  // skip the additional fluxes along the transverse directions for computing the CT electric field
   const real _dh          = (real)1.0 / dh;   //unuse

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1   //unuse
      const int TDir2 = (d+2)%3;    // transverse direction 2   //unuse

      int sizeB_i, sizeB_j, stride_fc_B;   //unuse (assigned in the switch below but only feed the unused idx_fc_B)
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i   = N_HF_FLUX-1;              size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = 0;                        j_offset = flux_offset;              k_offset    = flux_offset;
                  sizeB_i  = FLU_NXT_P1;               sizeB_j  = FLU_NXT;                  stride_fc_B = 1;
                  break;

         case 1 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-1;              size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = flux_offset;              j_offset = 0;                        k_offset    = flux_offset;
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT_P1;               stride_fc_B = FLU_NXT;
                  break;

         case 2 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-1;
                  i_offset = flux_offset;              j_offset = flux_offset;              k_offset    = 0;
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT;                  stride_fc_B = SQR(FLU_NXT);
                  break;
      } // switch ( d )

      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux   = idx % size_i           + i_offset;
         const int j_flux   = idx % size_ij / size_i + j_offset;
         const int k_flux   = idx / size_ij          + k_offset;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

//       conserved variable and cell-centered magnetic field index
         const int i_cvar   = i_flux;
         const int j_cvar   = j_flux;
         const int k_cvar   = k_flux;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

//       face-centered magnetic field index
         const int idx_fc_B = IDX321( i_cvar, j_cvar, k_cvar, sizeB_i, sizeB_j ) + stride_fc_B;   //unuse

//       get left and right cell indices
         const int idx_L = idx_cvar;
         const int idx_R = idx_cvar + didx_cvar[d];

//       1. get CR variables (left and right states) - using donor cell (1st order)
         const real Ec_L   = g_ConVar[CR_E ][idx_L];
         const real Ec_R   = g_ConVar[CR_E ][idx_R];
         const real Fc_L[3] = { g_ConVar[CR_F1][idx_L], g_ConVar[CR_F2][idx_L], g_ConVar[CR_F3][idx_L] };
         const real Fc_R[3] = { g_ConVar[CR_F1][idx_R], g_ConVar[CR_F2][idx_R], g_ConVar[CR_F3][idx_R] };

//       2. get density and velocity
         const real rho_L = g_ConVar[DENS][idx_L];
         const real rho_R = g_ConVar[DENS][idx_R];
         const real v_L   = g_ConVar[MOMX+d][idx_L] / rho_L;
         const real v_R   = g_ConVar[MOMX+d][idx_R] / rho_R;

//       3. get cell-centered magnetic field for left and right cells
         const real Bx_L = g_CC_B[0][idx_L];
         const real By_L = g_CC_B[1][idx_L];
         const real Bz_L = g_CC_B[2][idx_L];

         const real Bx_R = g_CC_B[0][idx_R];
         const real By_R = g_CC_B[1][idx_R];
         const real Bz_R = g_CC_B[2][idx_R];

//       4. READ sigma_adv from stored fields (computed at init or previous flux call)
//          Use cell-centered values directly (not averaged)
         const real sigma_adv_L = g_ConVar[ADV_SIGMA][idx_L];
         const real sigma_adv_R = g_ConVar[ADV_SIGMA][idx_R];

//       5. get CR parameters
         const real vmax = MicroPhy->CR_vmax;

//       6. compute v_diff at cell centers using the new helper function
//          vdiff_L is computed from LEFT cell's cell-centered data
//          vdiff_R is computed from RIGHT cell's cell-centered data
//          This matches Athena++ pattern: vdiff_l_(i) = v_diff[i-1], vdiff_r_(i) = v_diff[i]
         const real vdiff_L = CR_ComputeVdiff( sigma_adv_L, Bx_L, By_L, Bz_L, Ec_L, rho_L, vmax, dh, d, MicroPhy );
         const real vdiff_R = CR_ComputeVdiff( sigma_adv_R, Bx_R, By_R, Bz_R, Ec_R, rho_R, vmax, dh, d, MicroPhy );

//       7. compute HLLE flux
         real flux_E, flux_F[3];
         CR_ComputeHLLEFlux( Ec_L, Ec_R, Fc_L, Fc_R, v_L, v_R, vdiff_L, vdiff_R, vmax, d,
                             flux_E, flux_F );

//       8. add fluxes to output arrays
         g_Flux_Half[d][CR_E ][idx_flux] = flux_E;
         g_Flux_Half[d][CR_F1][idx_flux] = flux_F[0];
         g_Flux_Half[d][CR_F2][idx_flux] = flux_F[1];
         g_Flux_Half[d][CR_F3][idx_flux] = flux_F[2];

//       9. zero the ADV_* flux slots so the flux-divergence update passes these work arrays
//          through unchanged, instead of advecting them like physical passive scalars
//          --> the hydro Riemann solvers fill these slots with passive-advection fluxes,
//              but Athena++ never advects sigma_adv/v_adv (they are not part of u_cr);
//              this loop covers every face read by Hydro_RiemannPredict(), so dflux[ADV_*]=0
         g_Flux_Half[d][ADV_SIGMA][idx_flux] = (real)0.0;
         g_Flux_Half[d][ADV_VX   ][idx_flux] = (real)0.0;
         g_Flux_Half[d][ADV_VY   ][idx_flux] = (real)0.0;
         g_Flux_Half[d][ADV_VZ   ][idx_flux] = (real)0.0;

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)


} // FUNCTION : CR_TwoMomentFlux_HalfStep



//-----------------------------------------------------------------------------------------
// Function    : CR_TwoMomentFlux_FullStep
//
// Description : Compute full-step CR two-moment fluxes using reconstructed face-centered values
//
// Note        : 1. Uses g_FC_Var for reconstructed left/right states (Ec, Fc, rho, mom)
//               2. Uses g_PriVar_Half for cell-centered values needed by vdiff (B, sigma_adv)
//               3. Uses NFlux as the stride for flux array access
//               4. Overwrites the CR slots of g_FC_Flux[] computed by Hydro_ComputeFlux() and
//                  zeroes the ADV_* slots (sigma_adv/v_adv are work arrays recomputed each step,
//                  never advected --> matches Athena++, where they are not part of u_cr)
//
// Reference   : Athena++ cr_transport.cpp, cr_flux.cpp
//
// Parameter   : g_FC_Var      : Array storing the reconstructed face-centered conserved variables
//               g_PriVar_Half : Array storing the cell-centered half-step primitive variables
//                               (for B field and sigma_adv used by vdiff computation)
//               g_FC_Flux     : Array with hydrodynamic fluxes for adding the cosmic-ray diffusive fluxes
//               g_FC_B_Half   : Array storing the input face-centered, half-step magnetic field
//               NFlux         : Stride for accessing g_FC_Flux[]
//               NSkip_N       : Number of fluxes to skip at the beginning and end in the normal direction
//               NSkip_T       : Number of fluxes to skip at the beginning and end in the transverse direction
//               dh            : Cell size
//               MicroPhy      : Microphysics object
//
// Return      : g_FC_Flux[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentFlux_FullStep( const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                 const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                       real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                 const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],   //unuse
                                 const int NFlux, const int NSkip_N, const int NSkip_T,
                                 const real dh, const MicroPhy_t *MicroPhy )
{
   const int  didx_fc[3]   = { 1, N_FC_VAR, SQR(N_FC_VAR) };
   const int  didx_pvar[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const real _dh          = (real)1.0 / dh;   //unuse

// offset from g_FC_Var index to g_PriVar_Half index
// g_FC_Var has size N_FC_VAR, g_PriVar_Half has size N_HF_VAR
// For MHM_RP with MHD: N_FC_VAR = PS2+2 = 18, N_HF_VAR = FLU_NXT-2
   const int fc2pvar_offset = ( N_HF_VAR - N_FC_VAR ) / 2;

   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;        // left face index
      const int faceR = faceL + 1;  // right face index
      const int TDir1 = (d+1)%3;    // transverse direction 1   //unuse
      const int TDir2 = (d+2)%3;    // transverse direction 2   //unuse

      int sizeB_i, sizeB_j, stride_fc_B;   //unuse (assigned in the switch below but only feed the unused idx_fc_B)
      int idx_fc_s[3], idx_flux_e[3];

      switch ( d )
      {
         case 0 : idx_fc_s  [0] = NSkip_N;              idx_fc_s  [1] = NSkip_T;              idx_fc_s  [2] = NSkip_T;
                  idx_flux_e[0] = N_FC_VAR-1-2*NSkip_N; idx_flux_e[1] = N_FC_VAR-2*NSkip_T;   idx_flux_e[2] = N_FC_VAR-2*NSkip_T;
                  sizeB_i  = FLU_NXT_P1;                sizeB_j  = FLU_NXT;                   stride_fc_B = 1;
                  break;

         case 1 : idx_fc_s  [0] = NSkip_T;              idx_fc_s  [1] = NSkip_N;              idx_fc_s  [2] = NSkip_T;
                  idx_flux_e[0] = N_FC_VAR-2*NSkip_T;   idx_flux_e[1] = N_FC_VAR-1-2*NSkip_N; idx_flux_e[2] = N_FC_VAR-2*NSkip_T;
                  sizeB_i  = FLU_NXT;                   sizeB_j  = FLU_NXT_P1;                stride_fc_B = FLU_NXT;
                  break;

         case 2 : idx_fc_s  [0] = NSkip_T;              idx_fc_s  [1] = NSkip_T;              idx_fc_s  [2] = NSkip_N;
                  idx_flux_e[0] = N_FC_VAR-2*NSkip_T;   idx_flux_e[1] = N_FC_VAR-2*NSkip_T;   idx_flux_e[2] = N_FC_VAR-1-2*NSkip_N;
                  sizeB_i  = FLU_NXT;                   sizeB_j  = FLU_NXT;                   stride_fc_B = SQR(FLU_NXT);
                  break;
      } // switch ( d )

      const int size_ij = idx_flux_e[0]*idx_flux_e[1];

      CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
      {
//       flux index
         const int i_flux   = idx % idx_flux_e[0];
         const int j_flux   = idx % size_ij / idx_flux_e[0];
         const int k_flux   = idx / size_ij;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

//       face-centered variable index (for reconstructed values in g_FC_Var)
         const int i_fc     = i_flux + idx_fc_s[0];
         const int j_fc     = j_flux + idx_fc_s[1];
         const int k_fc     = k_flux + idx_fc_s[2];
         const int idx_fc   = IDX321( i_fc, j_fc, k_fc, N_FC_VAR, N_FC_VAR );

//       primitive variable index (for cell-centered values in g_PriVar_Half, used by vdiff)
         const int i_pvar   = i_fc + fc2pvar_offset;
         const int j_pvar   = j_fc + fc2pvar_offset;
         const int k_pvar   = k_fc + fc2pvar_offset;
         const int idx_pvar = IDX321( i_pvar, j_pvar, k_pvar, N_HF_VAR, N_HF_VAR );

//       face-centered B index
         const int idx_fc_B = IDX321( i_pvar, j_pvar, k_pvar, sizeB_i, sizeB_j ) + stride_fc_B;   //unuse

//       get left and right cell indices for cell-centered arrays (g_PriVar_Half)
         const int idx_pvar_L = idx_pvar;
         const int idx_pvar_R = idx_pvar + didx_pvar[d];

//       1. get CR variables from RECONSTRUCTED face-centered values
//          Left state:  right face of left cell  -> g_FC_Var[faceR][v][idx_fc]
//          Right state: left face of right cell  -> g_FC_Var[faceL][v][idx_fc + didx_fc[d]]
         const real Ec_L   = g_FC_Var[faceR][CR_E ][idx_fc];
         const real Ec_R   = g_FC_Var[faceL][CR_E ][idx_fc + didx_fc[d]];
         const real Fc_L[3] = { g_FC_Var[faceR][CR_F1][idx_fc], g_FC_Var[faceR][CR_F2][idx_fc], g_FC_Var[faceR][CR_F3][idx_fc] };
         const real Fc_R[3] = { g_FC_Var[faceL][CR_F1][idx_fc + didx_fc[d]], g_FC_Var[faceL][CR_F2][idx_fc + didx_fc[d]], g_FC_Var[faceL][CR_F3][idx_fc + didx_fc[d]] };

//       2. get density and velocity from RECONSTRUCTED face-centered values
//          g_FC_Var stores conserved variables (density and momentum), so compute velocity = mom/rho
         const real rho_L = g_FC_Var[faceR][DENS][idx_fc];
         const real rho_R = g_FC_Var[faceL][DENS][idx_fc + didx_fc[d]];
         const real v_L   = g_FC_Var[faceR][MOMX+d][idx_fc] / rho_L;
         const real v_R   = g_FC_Var[faceL][MOMX+d][idx_fc + didx_fc[d]] / rho_R;

//       3. get cell-centered magnetic field for left and right cells (for vdiff computation)
//          Cell-centered B is stored after MAG_OFFSET in primitive variables
//          Use donor-cell (cell-centered) values for vdiff, not reconstructed values
         const real Bx_L = g_PriVar_Half[MAG_OFFSET+0][idx_pvar_L];
         const real By_L = g_PriVar_Half[MAG_OFFSET+1][idx_pvar_L];
         const real Bz_L = g_PriVar_Half[MAG_OFFSET+2][idx_pvar_L];

         const real Bx_R = g_PriVar_Half[MAG_OFFSET+0][idx_pvar_R];
         const real By_R = g_PriVar_Half[MAG_OFFSET+1][idx_pvar_R];
         const real Bz_R = g_PriVar_Half[MAG_OFFSET+2][idx_pvar_R];

//       4. READ sigma_adv from stored fields (computed at init or previous flux call)
//          Use cell-centered values directly (not averaged) for vdiff computation
         const real sigma_adv_L = g_PriVar_Half[ADV_SIGMA][idx_pvar_L];
         const real sigma_adv_R = g_PriVar_Half[ADV_SIGMA][idx_pvar_R];

//       5. get CR parameters
         const real vmax = MicroPhy->CR_vmax;

//       6. compute v_diff at cell centers using the new helper function
//          vdiff_L is computed from LEFT cell's cell-centered data (donor-cell)
//          vdiff_R is computed from RIGHT cell's cell-centered data (donor-cell)
//          This matches Athena++ pattern: vdiff_l_(i) = v_diff[i-1], vdiff_r_(i) = v_diff[i]
//          Use cell-centered Ec and rho from g_PriVar_Half for vdiff (not reconstructed)
         const real Ec_cc_L  = g_PriVar_Half[CR_E][idx_pvar_L];
         const real Ec_cc_R  = g_PriVar_Half[CR_E][idx_pvar_R];
         const real rho_cc_L = g_PriVar_Half[DENS][idx_pvar_L];
         const real rho_cc_R = g_PriVar_Half[DENS][idx_pvar_R];
         const real vdiff_L = CR_ComputeVdiff( sigma_adv_L, Bx_L, By_L, Bz_L, Ec_cc_L, rho_cc_L, vmax, dh, d, MicroPhy );
         const real vdiff_R = CR_ComputeVdiff( sigma_adv_R, Bx_R, By_R, Bz_R, Ec_cc_R, rho_cc_R, vmax, dh, d, MicroPhy );

//       7. compute HLLE flux
         real flux_E, flux_F[3];
         CR_ComputeHLLEFlux( Ec_L, Ec_R, Fc_L, Fc_R, v_L, v_R, vdiff_L, vdiff_R, vmax, d,
                             flux_E, flux_F );

//       8. add fluxes to output arrays (replace, not add, since hydro doesn't handle CR)
         g_FC_Flux[d][CR_E ][idx_flux] = flux_E;
         g_FC_Flux[d][CR_F1][idx_flux] = flux_F[0];
         g_FC_Flux[d][CR_F2][idx_flux] = flux_F[1];
         g_FC_Flux[d][CR_F3][idx_flux] = flux_F[2];

//       9. zero the ADV_* flux slots so the flux-divergence update passes these work arrays
//          through unchanged, instead of advecting them like physical passive scalars
//          --> the hydro Riemann solvers fill these slots with passive-advection fluxes,
//              but Athena++ never advects sigma_adv/v_adv (they are not part of u_cr);
//              this loop covers every face read by Hydro_FullStepUpdate(), so dFlux[ADV_*]=0
         g_FC_Flux[d][ADV_SIGMA][idx_flux] = (real)0.0;
         g_FC_Flux[d][ADV_VX   ][idx_flux] = (real)0.0;
         g_FC_Flux[d][ADV_VY   ][idx_flux] = (real)0.0;
         g_FC_Flux[d][ADV_VZ   ][idx_flux] = (real)0.0;

      } // CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
   } // for (int d=0; d<3; d++)


} // FUNCTION : CR_TwoMomentFlux_FullStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_HalfStep
//
// Description : Apply implicit source term update for CR two-moment equations at half-step
//
// Note        : 1. The main flux divergence is already applied in the Hydro_RiemannPredict loop
//               2. This function applies the implicit source term solve at half-step with dt_source = 0.5*dt
//               3. Athena++ applies source terms at BOTH half-step and full-step for VL2 integrator
//               4. Gas state (rho, momentum) is read from the stage-updated OneCell[] (i.e. after the
//                  half-step flux divergence), matching Athena++'s use of the post-transport u in
//                  AddSourceTerms() at every stage
//               5. Back-reaction to gas momentum/energy is applied here as well (when CR_SOURCE is on),
//                  matching Athena++'s stage-agnostic AddSourceTerms(); the kicked half-step gas state
//                  then feeds the full-step reconstruction and the full-step source
//               6. The ec_source heating term uses grad(Pc) built from the CR flux divergence of the
//                  half-step fluxes (matching Athena++'s grad_pc_), with the gas velocity read from the
//                  t^n state (g_ConVar_In): Athena++ evaluates ec_source_ in CalculateFluxes() from the
//                  PRE-stage w, while the implicit solve uses the stage-updated u
//               7. The rotation angles use the t^n cell-centered B field (B_n[], recomputed by the
//                  caller from the t^n face-centered B), NOT the CT-half-step B in OneCell[MAG_OFFSET+*]:
//                  Athena++'s stage-1 source and ec_source use b_angle computed from the t^n bcc (set by
//                  the previous step's CRTC_OPACITY); the full-step source instead uses the half-step B
//                  (Athena++'s stage-2 bcc)
//
// Reference   : Athena++ cr_source.cpp, cr_transport.cpp, time_integrator.cpp
//
// Parameter   : OneCell     : Single-cell fluid array (already updated with flux divergence)
//               g_ConVar_In : Array storing the input cell-centered conserved variables
//               g_Flux_Half : Array storing the input face-centered fluxes
//               B_n         : t^n cell-centered B field [3] (see Note 7)
//               idx_in      : Index of accessing g_ConVar_In[]
//               didx_in     : Index increment of g_ConVar_In[]
//               idx_flux    : Index of accessing g_Flux_Half[]
//               didx_flux   : Index increment of g_Flux_Half[]
//               dt          : Full time step (source uses 0.5*dt for half-step)
//               dh          : Cell size
//               EoS         : EoS object
//               MicroPhy    : Microphysics object
//
// Return      : OneCell[] (modified CR and gas fields)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_HalfStep( real OneCell[NCOMP_TOTAL_PLUS_MAG],
                            const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                            const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real B_n[],
                            const int idx_in, const int didx_in[3],   //unuse: didx_in
                            const int idx_flux, const int didx_flux[3],
                            const real dt, const real dh, const EoS_t *EoS , const MicroPhy_t *MicroPhy )   //unuse: EoS
{
// The flux divergence update for CR_E, CR_F1, CR_F2, CR_F3 is already done
// in the main Hydro_RiemannPredict loop above where out_con is updated.

// CR parameters
   const real vmax   = MicroPhy->CR_vmax;
   const real invlim = (real)1.0 / vmax;

// Half-step uses dt_source = 0.5 * dt (matching Athena++ VL2 stage 1 with beta=0.5)
   const real dt_source = (real)0.5 * dt;

// 1. Get current CR state (already has flux divergence applied)
//    Floor CR_E before using it, matching Athena++'s post-transport floor (cr_transport.cpp:
//    CRE is floored to TINY_NUMBER after every transport step): the half-step flux divergence
//    can leave CR_E<0 in the strong-CR regime, and the generic passive-scalar floor in
//    Hydro_RiemannPredict() runs only AFTER this function; a negative ec would otherwise give
//    a negative Pc and be preserved by the "new_ec = ec_old" reset below
   real ec  = FMAX( OneCell[CR_E], TINY_NUMBER );
   real fc1 = OneCell[CR_F1];
   real fc2 = OneCell[CR_F2];
   real fc3 = OneCell[CR_F3];

// 2. Get gas density and velocity from the STAGE-UPDATED conserved variables (OneCell already has
//    the half-step flux divergence applied), matching Athena++'s post-transport u in AddSourceTerms()
//    --> guard rho with TINY_NUMBER since the density floor is applied only AFTER this function;
//        note Athena++ instead applies its EoS density floor here (cr_source.cpp), so floored
//        cells will not match
   const real rho = FMAX( OneCell[DENS], TINY_NUMBER );
   real v1 = OneCell[MOMX] / rho;
   real v2 = OneCell[MOMY] / rho;
   real v3 = OneCell[MOMZ] / rho;

// 3. Get B field: t^n cell-centered B passed in by the caller, matching Athena++'s stage-1 b_angle
//    computed from the t^n bcc — NOT the CT-half-step-updated OneCell[MAG_OFFSET+*] (see Note 7)
//    (MHD is compile-enforced for CR_STREAMING, see Aux_Check_Parameter.cpp)
   const real Bx = B_n[MAGX];
   const real By = B_n[MAGY];
   const real Bz = B_n[MAGZ];

// 4. Compute B-field angles for rotation
   real sint, cost, sinp, cosp;
   CR_ComputeBFieldAngles( Bx, By, Bz, sint, cost, sinp, cosp );

// 5. READ sigma_adv and v_adv from stored fields (updated by flux function)
//    This follows the Python/Athena++ pattern where add_source() READS sigma_adv/v_adv
   const real sigma_adv_para = g_ConVar_In[ADV_SIGMA][idx_in];
   const real v_adv_x = g_ConVar_In[ADV_VX][idx_in];
   const real v_adv_y = g_ConVar_In[ADV_VY][idx_in];
   const real v_adv_z = g_ConVar_In[ADV_VZ][idx_in];
   const real sigma_adv_perp = MicroPhy->CR_max_opacity;
   const bool CR_source      = MicroPhy->CR_source;      // flag to enable back-reaction to gas
   const bool CR_stream      = MicroPhy->CR_stream;      // flag to enable streaming
   const bool CR_Ec_source   = MicroPhy->CR_Ec_source;   // flag to include the CR energy source term

// Total velocity = gas velocity + streaming velocity (streaming added only if enabled)
   real vtot1 = v1;
   real vtot2 = v2;
   real vtot3 = v3;
   if ( CR_stream ) {
      vtot1 += v_adv_x;
      vtot2 += v_adv_y;
      vtot3 += v_adv_z;
   }

// 6. Save original CR energy and flux for the floor check and gas back-reaction
   const real ec_old  = ec;
   const real fc1_old = fc1;
   const real fc2_old = fc2;
   const real fc3_old = fc3;

// 7. Rotate all vectors to B-aligned frame
   real fr1 = fc1, fr2 = fc2, fr3 = fc3;

#  ifdef MHD
   RotateVec( sint, cost, sinp, cosp, v1, v2, v3 );
   RotateVec( sint, cost, sinp, cosp, fr1, fr2, fr3 );
   RotateVec( sint, cost, sinp, cosp, vtot1, vtot2, vtot3 );

// Perpendicular components have no streaming velocity contribution in B-frame
   vtot2 = (real)0.0;
   vtot3 = (real)0.0;
#  endif

// 8. Compute effective sigma: combine sigma_diff with sigma_adv only if streaming is enabled
   const real sigma_diff = MicroPhy->CR_sigma;
   const real sigma_diff_perp = MicroPhy->CR_sigma_perp;

   real sigma_x = sigma_diff;
   real sigma_y = sigma_diff_perp;
   real sigma_z = sigma_diff_perp;
   if ( CR_stream ) {
      sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_para );
      sigma_y = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
      sigma_z = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
   }

// 9. Build implicit matrix and solve
//     Source terms:
//     dEc/dt = -vtot · sigma · (Fc - v*Ec*(4/3)/vmax)
//     dFc/dt = -vmax · sigma · (Fc - v*Ec*(4/3)/vmax)

   const real rhs1 = ec;
   const real rhs2 = fr1;
   const real rhs3 = fr2;
   const real rhs4 = fr3;

// Coefficients for Ec equation (row 1)
   const real coef_11 = (real)1.0 - dt_source * sigma_x * vtot1 * v1 * invlim * ((real)4.0/(real)3.0)
                                  - dt_source * sigma_y * vtot2 * v2 * invlim * ((real)4.0/(real)3.0)
                                  - dt_source * sigma_z * vtot3 * v3 * invlim * ((real)4.0/(real)3.0);
   const real coef_12 = dt_source * sigma_x * vtot1;
   const real coef_13 = dt_source * sigma_y * vtot2;
   const real coef_14 = dt_source * sigma_z * vtot3;

// Coefficients for Fr1 equation (row 2)
   const real coef_21 = -dt_source * v1 * sigma_x * ((real)4.0/(real)3.0);
   const real coef_22 = (real)1.0 + dt_source * vmax * sigma_x;

// Coefficients for Fr2 equation (row 3)
   const real coef_31 = -dt_source * v2 * sigma_y * ((real)4.0/(real)3.0);
   const real coef_33 = (real)1.0 + dt_source * vmax * sigma_y;

// Coefficients for Fr3 equation (row 4)
   const real coef_41 = -dt_source * v3 * sigma_z * ((real)4.0/(real)3.0);
   const real coef_44 = (real)1.0 + dt_source * vmax * sigma_z;

// Solve by substitution (since flux equations are decoupled from each other)
   const real e_coef = coef_11 - coef_12 * coef_21 / coef_22 
                               - coef_13 * coef_31 / coef_33
                               - coef_14 * coef_41 / coef_44;

   real new_ec = rhs1 - coef_12 * rhs2 / coef_22
                      - coef_13 * rhs3 / coef_33
                      - coef_14 * rhs4 / coef_44;
   new_ec /= e_coef;

// CR_Ec_source=0 drops the CR energy source term so that the conservative system
// dEc/dt = -div(Fc) (Eq. 19 of Jiang & Oh 2018) is solved; the flux equations keep
// their full source terms
   if ( !CR_Ec_source )   new_ec = rhs1;

// Back-substitute to get new flux
   real newfr1 = ( rhs2 - coef_21 * new_ec ) / coef_22;
   real newfr2 = ( rhs3 - coef_31 * new_ec ) / coef_33;
   real newfr3 = ( rhs4 - coef_41 * new_ec ) / coef_44;

// 10. Rotate back to lab frame
#  ifdef MHD
   InvRotateVec( sint, cost, sinp, cosp, newfr1, newfr2, newfr3 );
#  endif

// 11. Compute perpendicular heating term (ec_source)
//     This is the work done by perpendicular gas flow: v_perp · grad(Pc)
//     grad(Pc) is built from the CR flux divergence of the half-step fluxes (same discretization
//     as CR_UpdateStreaming()), matching Athena++'s grad_pc_ used for ec_source_ (cr_transport.cpp):
//     each component n is the FULL divergence of the flux vector of the CR_F(n+1) equation,
//     summing all three flux directions (9 terms total)
#  ifdef MHD
   const real _dh = (real)1.0 / dh;
   const int  CRF_v[3] = { CR_F1, CR_F2, CR_F3 };  // CR_F* field indices DESCEND, hence the explicit list
   real grad_pc[3];
   for (int n=0; n<3; n++) {
      grad_pc[n] = (real)0.0;
      for (int d=0; d<3; d++)
         grad_pc[n] += g_Flux_Half[d][ CRF_v[n] ][idx_flux] - g_Flux_Half[d][ CRF_v[n] ][ idx_flux - didx_flux[d] ];
      grad_pc[n] *= _dh * invlim;
   }

   real dpcdx_B = grad_pc[0], dpcdy_B = grad_pc[1], dpcdz_B = grad_pc[2];
   RotateVec( sint, cost, sinp, cosp, dpcdx_B, dpcdy_B, dpcdz_B );

// Perpendicular velocity (in B-frame, gas velocity only)
// --> deliberately read from the t^n state: Athena++ evaluates ec_source_ in CalculateFluxes()
//     from the PRE-stage w, while the implicit solve above uses the stage-updated state
   const real rho_n = g_ConVar_In[DENS][idx_in];
   real v1_B = g_ConVar_In[MOMX][idx_in] / rho_n;
   real v2_B = g_ConVar_In[MOMY][idx_in] / rho_n;
   real v3_B = g_ConVar_In[MOMZ][idx_in] / rho_n;
   RotateVec( sint, cost, sinp, cosp, v1_B, v2_B, v3_B );

   const real ec_source = v2_B * dpcdy_B + v3_B * dpcdz_B;
   if ( CR_Ec_source )   new_ec += dt_source * ec_source;
#  endif

// 12. Floor CR energy
   if ( new_ec < TINY_NUMBER )
      new_ec = ec_old;

// 13. Apply back-reaction to gas momentum and energy, matching Athena++'s stage-agnostic
//     AddSourceTerms(); the kicked half-step gas state feeds the full-step reconstruction
//     (the effective source dt is already contained in newfr*/new_ec via dt_source)
   if ( CR_source ) {
//    momentum change: delta_p = -(new_Fc - old_Fc) / vmax
      OneCell[MOMX] += -( newfr1 - fc1_old ) * invlim;
      OneCell[MOMY] += -( newfr2 - fc2_old ) * invlim;
      OneCell[MOMZ] += -( newfr3 - fc3_old ) * invlim;

//    energy change: delta_E = -(new_Ec - old_Ec)
//    (CR energy lost goes to gas thermal energy)
      real new_eg = OneCell[ENGY] - ( new_ec - ec_old );
      if ( new_eg < (real)0.0 )
         new_eg = OneCell[ENGY];
      OneCell[ENGY] = new_eg;
   }

// 14. Update CR fields
   OneCell[CR_E ] = new_ec;
   OneCell[CR_F1] = newfr1;
   OneCell[CR_F2] = newfr2;
   OneCell[CR_F3] = newfr3;

} // FUNCTION : CR_TwoMomentSource_HalfStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_FullStep
//
// Description : Apply implicit source term update for CR two-moment equations
//
// Note        : 1. Solves the coupled source terms for Ec and Fc implicitly
//               2. Includes rotation to B-aligned frame for anisotropic transport
//               3. Applies back-reaction to gas momentum and energy if CR_source is enabled
//               4. Gas state (rho, momentum) is read from the stage-updated g_Output[] (i.e. after the
//                  full-step flux divergence in Hydro_FullStepUpdate()), matching Athena++'s use of the
//                  post-transport u in AddSourceTerms() at every stage; B field and sigma_adv/v_adv are
//                  still read from the half-step g_PriVar_Half[], matching Athena++'s stage-2 bcc
//               5. The ec_source heating term uses grad(Pc) built from the CR flux divergence of the
//                  full-step fluxes (matching Athena++'s grad_pc_), with the gas velocity read from the
//                  half-step primitives (g_PriVar_Half): Athena++ evaluates ec_source_ in
//                  CalculateFluxes() from the PRE-stage w, while the implicit solve uses the
//                  stage-updated u
//               6. Based on Athena++ cr_source.cpp
//
// Reference   : Athena++ cr_source.cpp, cr_transport.cpp
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered primitive variables
//               g_Output      : Array to store the updated fluid data (already has flux divergence applied)
//               g_Flux        : Array storing the input face-centered fluxes
//               g_FC_Var      : Array storing the input face-centered conserved variables
//               dt            : Time interval to advance solution
//               dh            : Cell size
//               EoS           : EoS object
//               MicroPhy      : Microphysics object
//
// Return      : g_Output[] with updated CR and gas fields
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_Output[][ CUBE(PS2) ],
                                const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],   //unuse
                                const real dt, const real dh, const EoS_t *EoS, const MicroPhy_t *MicroPhy )   //unuse: EoS
{
   const int  didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };

// offset between output array and primitive array
   const int  pvar_offset  = ( N_HF_VAR - PS2 ) / 2;

// CR parameters
   const real vmax   = MicroPhy->CR_vmax;
   const real invlim = (real)1.0 / vmax;
   const bool CR_source = MicroPhy->CR_source;         // flag to enable back-reaction to gas
   const bool CR_stream = MicroPhy->CR_stream;         // flag to enable streaming
   const bool CR_Ec_source = MicroPhy->CR_Ec_source;   // flag to include the CR energy source term

   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
//    output array indices
      const int i_out   = idx_out % PS2;
      const int j_out   = idx_out % SQR(PS2) / PS2;
      const int k_out   = idx_out / SQR(PS2);

//    primitive variable indices (with offset)
      const int i_pvar  = i_out + pvar_offset;
      const int j_pvar  = j_out + pvar_offset;
      const int k_pvar  = k_out + pvar_offset;
      const int idx_pvar = IDX321( i_pvar, j_pvar, k_pvar, N_HF_VAR, N_HF_VAR );

//    1. get current CR state from g_Output (which already has flux divergence applied)
//       CR_E is already non-negative here in practice: it is registered with FLOOR_YES, so
//       Hydro_FullStepUpdate() floors it right after the flux-divergence update (matching
//       Athena++'s post-transport floor in cr_transport.cpp); the FMAX below is kept as a
//       cheap defensive floor so this function stays correct on its own
      real ec  = FMAX( g_Output[CR_E ][idx_out], TINY_NUMBER );
      real fc1 = g_Output[CR_F1][idx_out];
      real fc2 = g_Output[CR_F2][idx_out];
      real fc3 = g_Output[CR_F3][idx_out];

//    2. get gas density and velocity from the STAGE-UPDATED conserved variables (g_Output already
//       has the full-step flux divergence applied), matching Athena++'s post-transport u in
//       AddSourceTerms()
//       --> guard rho with TINY_NUMBER since Hydro_FullStepUpdate() deliberately does not floor
//           the density (deferred to Flu_Close()); note Athena++ instead applies its EoS density
//           floor here (cr_source.cpp), so floored cells will not match
      const real rho = FMAX( g_Output[DENS][idx_out], TINY_NUMBER );
      real v1 = g_Output[MOMX][idx_out] / rho;
      real v2 = g_Output[MOMY][idx_out] / rho;
      real v3 = g_Output[MOMZ][idx_out] / rho;

//    3. get cell-centered B field from half-step primitive variables
//       (MHD is compile-enforced for CR_STREAMING, see Aux_Check_Parameter.cpp)
      const real Bx = g_PriVar_Half[MAG_OFFSET+MAGX][idx_pvar];
      const real By = g_PriVar_Half[MAG_OFFSET+MAGY][idx_pvar];
      const real Bz = g_PriVar_Half[MAG_OFFSET+MAGZ][idx_pvar];

//    4. compute B-field angles for rotation
      real sint, cost, sinp, cosp;
      CR_ComputeBFieldAngles( Bx, By, Bz, sint, cost, sinp, cosp );

//    5. READ sigma_adv and v_adv from stored fields (updated by flux function)
//       This follows the Python/Athena++ pattern where add_source() READS sigma_adv/v_adv
      const real sigma_adv_para = g_PriVar_Half[ADV_SIGMA][idx_pvar];
      const real v_adv_x = g_PriVar_Half[ADV_VX][idx_pvar];
      const real v_adv_y = g_PriVar_Half[ADV_VY][idx_pvar];
      const real v_adv_z = g_PriVar_Half[ADV_VZ][idx_pvar];
      const real sigma_adv_perp = MicroPhy->CR_max_opacity;

//    total velocity = gas velocity + streaming velocity (streaming added only if enabled)
      real vtot1 = v1;
      real vtot2 = v2;
      real vtot3 = v3;
      if ( CR_stream ) {
         vtot1 += v_adv_x;
         vtot2 += v_adv_y;
         vtot3 += v_adv_z;
      }

//    6. save original CR flux for back-reaction calculation
      const real fc1_old = fc1;
      const real fc2_old = fc2;
      const real fc3_old = fc3;
      const real ec_old  = ec;

//    7. rotate all vectors to B-aligned frame
      real fr1 = fc1, fr2 = fc2, fr3 = fc3;

#     ifdef MHD
      RotateVec( sint, cost, sinp, cosp, v1, v2, v3 );
      RotateVec( sint, cost, sinp, cosp, fr1, fr2, fr3 );
      RotateVec( sint, cost, sinp, cosp, vtot1, vtot2, vtot3 );

//    perpendicular components have no streaming velocity contribution in B-frame
      vtot2 = (real)0.0;
      vtot3 = (real)0.0;
#     endif

//    8. compute effective sigma: combine sigma_diff with sigma_adv only if streaming is enabled
      const real sigma_diff = MicroPhy->CR_sigma;
      const real sigma_diff_perp = MicroPhy->CR_sigma_perp;

      real sigma_x = sigma_diff;
      real sigma_y = sigma_diff_perp;
      real sigma_z = sigma_diff_perp;
      if ( CR_stream ) {
         sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_para );
         sigma_y = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
         sigma_z = (real)1.0 / ( (real)1.0/sigma_diff_perp + (real)1.0/sigma_adv_perp );
      }

//    9. build implicit matrix and solve
//        Source terms:
//        dEc/dt = -vtot · sigma · (Fc - v*Ec*(4/3)/vmax)
//        dFc/dt = -vmax · sigma · (Fc - v*Ec*(4/3)/vmax)
//
//        Matrix form: A * [Ec_new, Fr1_new, Fr2_new, Fr3_new]^T = [rhs1, rhs2, rhs3, rhs4]^T

      const real rhs1 = ec;
      const real rhs2 = fr1;
      const real rhs3 = fr2;
      const real rhs4 = fr3;

//    coefficients for Ec equation (row 1)
      const real coef_11 = (real)1.0 - dt * sigma_x * vtot1 * v1 * invlim * ((real)4.0/(real)3.0)
                                     - dt * sigma_y * vtot2 * v2 * invlim * ((real)4.0/(real)3.0)
                                     - dt * sigma_z * vtot3 * v3 * invlim * ((real)4.0/(real)3.0);
      const real coef_12 = dt * sigma_x * vtot1;
      const real coef_13 = dt * sigma_y * vtot2;
      const real coef_14 = dt * sigma_z * vtot3;

//    coefficients for Fr1 equation (row 2)
      const real coef_21 = -dt * v1 * sigma_x * ((real)4.0/(real)3.0);
      const real coef_22 = (real)1.0 + dt * vmax * sigma_x;

//    coefficients for Fr2 equation (row 3)
      const real coef_31 = -dt * v2 * sigma_y * ((real)4.0/(real)3.0);
      const real coef_33 = (real)1.0 + dt * vmax * sigma_y;

//    coefficients for Fr3 equation (row 4)
      const real coef_41 = -dt * v3 * sigma_z * ((real)4.0/(real)3.0);
      const real coef_44 = (real)1.0 + dt * vmax * sigma_z;

//    solve by substitution (since flux equations are decoupled from each other)
//    Fr_new = (rhs - coef_*1 * Ec_new) / coef_**
//    Substitute into Ec equation and solve for Ec_new

      const real e_coef = coef_11 - coef_12 * coef_21 / coef_22 
                                  - coef_13 * coef_31 / coef_33
                                  - coef_14 * coef_41 / coef_44;

      real new_ec = rhs1 - coef_12 * rhs2 / coef_22
                        - coef_13 * rhs3 / coef_33
                        - coef_14 * rhs4 / coef_44;
      new_ec /= e_coef;

//    CR_Ec_source=0 drops the CR energy source term so that the conservative system
//    dEc/dt = -div(Fc) (Eq. 19 of Jiang & Oh 2018) is solved; the flux equations keep
//    their full source terms
      if ( !CR_Ec_source )   new_ec = rhs1;

//    back-substitute to get new flux
      real newfr1 = ( rhs2 - coef_21 * new_ec ) / coef_22;
      real newfr2 = ( rhs3 - coef_31 * new_ec ) / coef_33;
      real newfr3 = ( rhs4 - coef_41 * new_ec ) / coef_44;

//    10. rotate back to lab frame
#     ifdef MHD
      InvRotateVec( sint, cost, sinp, cosp, newfr1, newfr2, newfr3 );
#     endif

//    11. compute perpendicular heating term (ec_source)
//        This is the work done by perpendicular gas flow: v_perp · grad(Pc)
//        In B-aligned frame: v_perp = (0, v2, v3), grad_pc_perp = (0, dPc/dy', dPc/dz')
//        grad(Pc) is built from the CR flux divergence of the full-step fluxes (same
//        discretization as CR_UpdateStreaming()), matching Athena++'s grad_pc_ used for
//        ec_source_ (cr_transport.cpp): each component n is the FULL divergence of the flux
//        vector of the CR_F(n+1) equation, summing all three flux directions (9 terms total)
#     ifdef MHD
      const real _dh = (real)1.0 / dh;

//    flux index: one flux ring is skipped along each transverse direction for the CT electric
//    field, hence the +1 offset (same mapping as Hydro_FullStepUpdate() with MHD, which is
//    compile-enforced for CR_STREAMING)
      const int idx_flux = IDX321( i_out+1, j_out+1, k_out+1, N_FL_FLUX, N_FL_FLUX );

      const int CRF_v[3] = { CR_F1, CR_F2, CR_F3 };  // CR_F* field indices DESCEND, hence the explicit list
      real grad_pc[3];
      for (int n=0; n<3; n++) {
         grad_pc[n] = (real)0.0;
         for (int d=0; d<3; d++)
            grad_pc[n] += g_Flux[d][ CRF_v[n] ][idx_flux] - g_Flux[d][ CRF_v[n] ][ idx_flux - didx_flux[d] ];
         grad_pc[n] *= _dh * invlim;
      }

      real dpcdx_B = grad_pc[0], dpcdy_B = grad_pc[1], dpcdz_B = grad_pc[2];
      RotateVec( sint, cost, sinp, cosp, dpcdx_B, dpcdy_B, dpcdz_B );

//    perpendicular velocity (in B-frame, gas velocity only)
//    Note: primitive array layout is [DENS=0, VelX=1, VelY=2, VelZ=3, ...]
//    --> deliberately read from the half-step primitives: Athena++ evaluates ec_source_ in
//        CalculateFluxes() from the PRE-stage w, while the implicit solve above uses the
//        stage-updated state
      real v1_B = g_PriVar_Half[1][idx_pvar];
      real v2_B = g_PriVar_Half[2][idx_pvar];
      real v3_B = g_PriVar_Half[3][idx_pvar];
      RotateVec( sint, cost, sinp, cosp, v1_B, v2_B, v3_B );

      const real ec_source = v2_B * dpcdy_B + v3_B * dpcdz_B;
      if ( CR_Ec_source )   new_ec += dt * ec_source;
#     endif

//    12. floor CR energy
      if ( new_ec < TINY_NUMBER )
         new_ec = ec_old;

//    13. apply back-reaction to gas momentum and energy
      if ( CR_source ) {
//       momentum change: delta_p = -(new_Fc - old_Fc) / vmax
         g_Output[MOMX][idx_out] += -( newfr1 - fc1_old ) * invlim;
         g_Output[MOMY][idx_out] += -( newfr2 - fc2_old ) * invlim;
         g_Output[MOMZ][idx_out] += -( newfr3 - fc3_old ) * invlim;

//       energy change: delta_E = -(new_Ec - old_Ec)
//       (CR energy lost goes to gas thermal energy)
         real new_eg = g_Output[ENGY][idx_out] - ( new_ec - ec_old );
         if ( new_eg < (real)0.0 )
            new_eg = g_Output[ENGY][idx_out];
         g_Output[ENGY][idx_out] = new_eg;
      }

//    14. update CR fields
      g_Output[CR_E ][idx_out] = new_ec;
      g_Output[CR_F1][idx_out] = newfr1;
      g_Output[CR_F2][idx_out] = newfr2;
      g_Output[CR_F3][idx_out] = newfr3;

   } // CGPU_LOOP( idx_out, CUBE(PS2) )

} // FUNCTION : CR_TwoMomentSource_FullStep

#endif // #ifdef CR_STREAMING
