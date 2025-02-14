#include "CUAPI.h"
#include "CUFLU.h"

#ifdef GPU



extern real (*d_Flu_Array_F_In )[FLU_NIN ][ CUBE(FLU_NXT) ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ CUBE(PS2) ];
extern real (*d_Flux_Array)[9][NFLUX_TOTAL][ SQR(PS2) ];
#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ CUBE(USG_NXT_F) ];
extern double (*d_Corner_Array_F)[3];
#endif
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_F_Out)[ CUBE(PS2) ];
#endif
#ifdef MHD
extern real (*d_Mag_Array_F_In )[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_Mag_Array_F_Out)[NCOMP_MAG][ PS2P1*SQR(PS2)          ];
extern real (*d_Ele_Array      )[9][NCOMP_ELE][ PS2P1*PS2 ];
extern real (*d_Mag_Array_T)[NCOMP_MAG][ PS1P1*SQR(PS1) ];
extern real (*d_Mag_Array_S_In)[NCOMP_MAG][ SRC_NXT_P1*SQR(SRC_NXT) ];
#endif
extern real *d_dt_Array_T;
extern real (*d_Flu_Array_T)[FLU_NIN_T][ CUBE(PS1) ];
extern real (*d_Flu_Array_S_In )[FLU_NIN_S ][ CUBE(SRC_NXT) ];
extern real (*d_Flu_Array_S_Out)[FLU_NOUT_S][ CUBE(PS1)     ];
extern double (*d_Corner_Array_S)[3];
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ];
extern real (*d_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ];
extern real (*d_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ];
extern real (*d_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ];
#ifdef MHD
extern real (*d_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ];
#endif
#endif // FLU_SCHEME

#if ( MODEL == ELBDM )
extern bool (*d_IsCompletelyRefined);
#endif

#if ( ELBDM_SCHEME == ELBDM_HYBRID )
extern bool (*d_HasWaveCounterpart)[ CUBE(HYB_NXT) ];
#endif

#if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
extern gramfe_matmul_float (*d_Flu_TimeEvo)[ 2*FLU_NXT ];
#endif

#if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#endif

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree_Fluid
// Description :  Free the GPU and CPU memory previously allocated by CUAPI_MemAllocate_Fluid()
//
// Parameter   :  GPU_NStream : Number of CUDA streams for the asynchronous memory copy
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemFree_Fluid( const int GPU_NStream )
{

// free the device memory
   if ( d_Flu_Array_F_In      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_F_In      )  );  d_Flu_Array_F_In      = NULL; }
   if ( d_Flu_Array_F_Out     != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_F_Out     )  );  d_Flu_Array_F_Out     = NULL; }
   if ( d_Flux_Array          != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flux_Array          )  );  d_Flux_Array          = NULL; }
#  ifdef UNSPLIT_GRAVITY
   if ( d_Pot_Array_USG_F     != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Pot_Array_USG_F     )  );  d_Pot_Array_USG_F     = NULL; }
   if ( d_Corner_Array_F      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Corner_Array_F      )  );  d_Corner_Array_F      = NULL; }
#  endif
#  ifdef DUAL_ENERGY
   if ( d_DE_Array_F_Out      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_DE_Array_F_Out      )  );  d_DE_Array_F_Out      = NULL; }
#  endif
#  ifdef MHD
   if ( d_Mag_Array_F_In      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Mag_Array_F_In      )  );  d_Mag_Array_F_In      = NULL; }
   if ( d_Mag_Array_F_Out     != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Mag_Array_F_Out     )  );  d_Mag_Array_F_Out     = NULL; }
   if ( d_Ele_Array           != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Ele_Array           )  );  d_Ele_Array           = NULL; }
   if ( d_Mag_Array_T         != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Mag_Array_T         )  );  d_Mag_Array_T         = NULL; }
   if ( d_Mag_Array_S_In      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Mag_Array_S_In      )  );  d_Mag_Array_S_In      = NULL; }
#  endif
   if ( d_dt_Array_T          != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_dt_Array_T          )  );  d_dt_Array_T          = NULL; }
   if ( d_Flu_Array_T         != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_T         )  );  d_Flu_Array_T         = NULL; }
   if ( d_Flu_Array_S_In      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_S_In      )  );  d_Flu_Array_S_In      = NULL; }
   if ( d_Flu_Array_S_Out     != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_S_Out     )  );  d_Flu_Array_S_Out     = NULL; }
   if ( d_Corner_Array_S      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Corner_Array_S      )  );  d_Corner_Array_S      = NULL; }
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   if ( d_PriVar              != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_PriVar              )  );  d_PriVar              = NULL; }
   if ( d_Slope_PPM           != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Slope_PPM           )  );  d_Slope_PPM           = NULL; }
   if ( d_FC_Var              != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_FC_Var              )  );  d_FC_Var              = NULL; }
   if ( d_FC_Flux             != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_FC_Flux             )  );  d_FC_Flux             = NULL; }
#  ifdef MHD
   if ( d_FC_Mag_Half         != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_FC_Mag_Half         )  );  d_FC_Mag_Half         = NULL; }
   if ( d_EC_Ele              != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_EC_Ele              )  );  d_EC_Ele              = NULL; }
#  endif
#  endif // FLU_SCHEME

#  if ( MODEL == ELBDM )
   if ( d_IsCompletelyRefined != NULL ) {  CUDA_CHECK_ERROR (  cudaFree( d_IsCompletelyRefined)  );  d_IsCompletelyRefined = NULL; }
#  endif

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( d_HasWaveCounterpart  != NULL ) {  CUDA_CHECK_ERROR (  cudaFree( d_HasWaveCounterpart )  );  d_HasWaveCounterpart  = NULL; }
#  endif

#  if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   if ( d_Flu_TimeEvo         != NULL ) {  CUDA_CHECK_ERROR (   cudaFree ( d_Flu_TimeEvo      )  );  d_Flu_TimeEvo         = NULL; }
#  endif

#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#    warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif


// free the host memory allocated by CUDA
   for (int t=0; t<2; t++)
   {
      if ( h_Flu_Array_F_In     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_F_In     [t] )  );  h_Flu_Array_F_In     [t] = NULL; }
      if ( h_Flu_Array_F_Out    [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_F_Out    [t] )  );  h_Flu_Array_F_Out    [t] = NULL; }
      if ( h_Flux_Array         [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flux_Array         [t] )  );  h_Flux_Array         [t] = NULL; }
#     ifdef UNSPLIT_GRAVITY
      if ( h_Pot_Array_USG_F    [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Pot_Array_USG_F    [t] )  );  h_Pot_Array_USG_F    [t] = NULL; }
      if ( h_Corner_Array_F     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Corner_Array_F     [t] )  );  h_Corner_Array_F     [t] = NULL; }
#     endif
#     ifdef DUAL_ENERGY
      if ( h_DE_Array_F_Out     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_DE_Array_F_Out     [t] )  );  h_DE_Array_F_Out     [t] = NULL; }
#     endif
#     ifdef MHD
      if ( h_Mag_Array_F_In     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Mag_Array_F_In     [t] )  );  h_Mag_Array_F_In     [t] = NULL; }
      if ( h_Mag_Array_F_Out    [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Mag_Array_F_Out    [t] )  );  h_Mag_Array_F_Out    [t] = NULL; }
      if ( h_Ele_Array          [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Ele_Array          [t] )  );  h_Ele_Array          [t] = NULL; }
      if ( h_Mag_Array_T        [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Mag_Array_T        [t] )  );  h_Mag_Array_T        [t] = NULL; }
      if ( h_Mag_Array_S_In     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Mag_Array_S_In     [t] )  );  h_Mag_Array_S_In     [t] = NULL; }
#     endif
      if ( h_dt_Array_T         [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_dt_Array_T         [t] )  );  h_dt_Array_T         [t] = NULL; }
      if ( h_Flu_Array_T        [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_T        [t] )  );  h_Flu_Array_T        [t] = NULL; }
      if ( h_Flu_Array_S_In     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_S_In     [t] )  );  h_Flu_Array_S_In     [t] = NULL; }
      if ( h_Flu_Array_S_Out    [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_S_Out    [t] )  );  h_Flu_Array_S_Out    [t] = NULL; }
      if ( h_Corner_Array_S     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Corner_Array_S     [t] )  );  h_Corner_Array_S     [t] = NULL; }

#     if ( MODEL == ELBDM )
      if ( h_IsCompletelyRefined[t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_IsCompletelyRefined[t] )  );  h_IsCompletelyRefined[t] = NULL; }
#     endif

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( h_HasWaveCounterpart [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_HasWaveCounterpart [t] )  ); h_HasWaveCounterpart  [t] = NULL; }
#     endif
   } // for (int t=0; t<2; t++)

#  if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   if ( h_GramFE_TimeEvo != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost ( h_GramFE_TimeEvo )  );  h_GramFE_TimeEvo = NULL; }
#  endif

// destroy streams
   if ( Stream != NULL )
   {
      for (int s=0; s<GPU_NStream; s++)   CUDA_CHECK_ERROR(  cudaStreamDestroy( Stream[s] )  );

      delete [] Stream;
      Stream = NULL;
   }

} // FUNCTION : CUAPI_MemFree_Fluid



#endif // #ifdef GPU
