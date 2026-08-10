#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "Global.h"
#include "CUDA_CheckError.h"
#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"

extern real *d_Turb_AccTable[2];

#endif // #ifdef __CUDACC__


// local function prototypes
#ifndef __CUDACC__

void Src_SetAuxArray_Turbulence( double [], int [] );
void Src_SetCPUFunc_Turbulence( SrcFunc_t & );
#ifdef GPU
void Src_SetGPUFunc_Turbulence( SrcFunc_t & );
#endif
void Src_SetConstMemory_Turbulence( const double AuxArray_Flt[], const int AuxArray_Int[],
                                    double *&DevPtr_Flt, int *&DevPtr_Int );
void Src_PassData2GPU_Turbulence();

#endif



/********************************************************
1. Turbulence source term
   --> Enabled by the runtime option "SRC_TURBULENCE"

2. This file is shared by both CPU and GPU

   CUSRC_Src_Turbulence.cu -> CPU_Src_Turbulence.cpp

3. Four steps are required to implement a source term

   I.   Set auxiliary arrays
   II.  Implement the source-term function
   III. [Optional] Add the work to be done every time
        before calling the major source-term function
   IV.  Set initialization functions

4. The source-term function must be thread-safe and
   not use any global variable
********************************************************/



// =======================
// I. Set auxiliary arrays
// =======================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetAuxArray_Turbulence
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by Src_Init_Turbulence()
//                2. AuxArray_Flt/Int[] have the size of SRC_NAUX_DLEP defined in Macro.h (default = 5)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_SetAuxArray_Turbulence( double AuxArray_Flt[], int AuxArray_Int[] )
{

// TBF

} // FUNCTION : Src_SetAuxArray_Turbulence
#endif // #ifndef __CUDACC__



// ======================================
// II. Implement the source-term function
// ======================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Turbulence
// Description :  Major source-term function
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. See Src_SetAuxArray_Turbulence() for the values stored in AuxArray_Flt/Int[]
//                3. Shared by both CPU and GPU
//
// Parameter   :  fluid             : Fluid array storing both the input and updated values
//                                    --> Including both active and passive variables
//                B                 : Cell-centered magnetic field
//                SrcTerms          : Structure storing all source-term variables
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                x/y/z             : Target physical coordinates
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                EoS               : EoS object
//                AuxArray_*        : Auxiliary arrays (see the Note above)
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void Src_Turbulence( real fluid[], const real B[],
                            const SrcTerms_t *SrcTerms, const real dt, const real dh,
                            const double x, const double y, const double z,
                            const double TimeNew, const double TimeOld,
                            const real MinDens, const real MinPres, const real MinEint, const long PassiveFloor,
                            const EoS_t *EoS, const double AuxArray_Flt[], const int AuxArray_Int[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif

// TBF

} // FUNCTION : Src_Turbulence



// ==================================================
// III. [Optional] Add the work to be done every time
//      before calling the major source-term function
// ==================================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_WorkBeforeMajorFunc_Turbulence
// Description :  Specify work to be done every time before calling the major source-term function
//
// Note        :  1. Invoked by Src_WorkBeforeMajorFunc()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  lv               : Target refinement level
//                TimeNew          : Target physical time to reach
//                TimeOld          : Physical time before update
//                                   --> The major source-term function will update the system from TimeOld to TimeNew
//                dt               : Time interval to advance solution
//                                   --> Physical coordinates : TimeNew - TimeOld == dt
//                                       Comoving coordinates : TimeNew - TimeOld == delta(scale factor) != dt
//                AuxArray_Flt/Int : Auxiliary arrays
//                                   --> Can be used and/or modified here
//                                   --> Must call Src_SetConstMemory_Turbulence() after modification
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_WorkBeforeMajorFunc_Turbulence( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                         double AuxArray_Flt[], int AuxArray_Int[] )
{

   // check update

} // FUNCTION : Src_WorkBeforeMajorFunc_Turbulence
#endif



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_PassData2GPU_Turbulence
// Description :  Transfer data to GPU
//
// Note        :  1. Invoked by Src_WorkBeforeMajorFunc_Turbulence()
//                2. Use synchronous transfer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_PassData2GPU_Turbulence()
{

   const long Size_Data   = sizeof(real)*    ;

// use synchronous transfer
   CUDA_CHECK_ERROR(  cudaMemcpy( d_Turb_AccTable,   d_Turb_AccTable,   Size_Data,   cudaMemcpyHostToDevice )  );

} // FUNCTION : Src_PassData2GPU_Turbulence
#endif // #ifdef __CUDACC__



// ================================
// IV. Set initialization functions
// ================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE SrcFunc_t SrcFunc_Ptr = Src_Turbulence;

//-----------------------------------------------------------------------------------------
// Function    :  Src_SetCPU/GPUFunc_Turbulence
// Description :  Return the function pointer of the CPU/GPU source-term function
//
// Note        :  1. Invoked by Src_Init_Turbulence()
//                2. Call-by-reference
//
// Parameter   :  SrcFunc_CPU/GPUPtr : CPU/GPU function pointer to be set
//
// Return      :  SrcFunc_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void Src_SetGPUFunc_Turbulence( SrcFunc_t &SrcFunc_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &SrcFunc_GPUPtr, SrcFunc_Ptr, sizeof(SrcFunc_t) )  );
}

#else

void Src_SetCPUFunc_Turbulence( SrcFunc_t &SrcFunc_CPUPtr )
{
   SrcFunc_CPUPtr = SrcFunc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetConstMemory_Turbulence
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Src_Init_Turbulence() and, if necessary, Src_WorkBeforeMajorFunc_Turbulence()
//                3. SRC_NAUX_DLEP is defined in Macro.h
//
// Parameter   :  AuxArray_Flt/Int : Auxiliary arrays to be copied to the constant memory
//                DevPtr_Flt/Int   : Pointers to store the addresses of constant memory arrays
//
// Return      :  c_Src_Dlep_AuxArray_Flt[], c_Src_Dlep_AuxArray_Int[], DevPtr_Flt, DevPtr_Int
//---------------------------------------------------------------------------------------------------
void Src_SetConstMemory_Turbulence( const double AuxArray_Flt[], const int AuxArray_Int[],
                                    double *&DevPtr_Flt, int *&DevPtr_Int )
{

// copy data to constant memory
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_Dlep_AuxArray_Flt, AuxArray_Flt, SRC_NAUX_DLEP*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_Dlep_AuxArray_Int, AuxArray_Int, SRC_NAUX_DLEP*sizeof(int   ) )  );

// obtain the constant-memory pointers
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Flt, c_Src_Dlep_AuxArray_Flt )  );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Int, c_Src_Dlep_AuxArray_Int )  );

} // FUNCTION : Src_SetConstMemory_Turbulence
#endif // #ifdef __CUDACC__



#ifndef __CUDACC__

//-----------------------------------------------------------------------------------------
// Function    :  Src_Init_Turbulence
// Description :  Initialize the turbulence source term
//
// Note        :  1. Set auxiliary arrays by invoking Src_SetAuxArray_*()
//                   --> Copy to the GPU constant memory and store the associated addresses
//                2. Set the source-term function by invoking Src_SetCPU/GPUFunc_*()
//                3. Invoked by Src_Init()
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_Init_Turbulence()
{

// set the auxiliary arrays
   Src_SetAuxArray_Turbulence( Src_Turb_AuxArray_Flt, Src_Turb_AuxArray_Int );

// copy the auxiliary arrays to the GPU constant memory and store the associated addresses
#  ifdef GPU
   Src_SetConstMemory_Turbulence( Src_Turb_AuxArray_Flt, Src_Turb_AuxArray_Int,
                                  SrcTerms.Turb_AuxArrayDevPtr_Flt, SrcTerms.Turb_AuxArrayDevPtr_Int );
#  else
   SrcTerms.Turb_AuxArrayDevPtr_Flt = Src_Turb_AuxArray_Flt;
   SrcTerms.Turb_AuxArrayDevPtr_Int = Src_Turb_AuxArray_Int;
#  endif

// set the major source-term function
   Src_SetCPUFunc_Turbulence( SrcTerms.Turb_CPUPtr );

#  ifdef GPU
   Src_SetGPUFunc_Turbulence( SrcTerms.Turb_GPUPtr );
   SrcTerms.Turb_FuncPtr = SrcTerms.Turb_GPUPtr;
#  else
   SrcTerms.Turb_FuncPtr = SrcTerms.Turb_CPUPtr;
#  endif

} // FUNCTION : Src_Init_Turbulence



//-----------------------------------------------------------------------------------------
// Function    :  Src_End_Turbulence
// Description :  Release the resources used by the Turbulence source term
//
// Note        :  1. Invoked by Src_End()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_End_Turbulence()
{

// TBF

} // FUNCTION : Src_End_Turbulence

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
