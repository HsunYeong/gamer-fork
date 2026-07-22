#ifndef __TURBULENCE_H__
#define __TURBULENCE_H__



#include "Macro.h"

#ifndef TURBULENCE
#  error : ERROR : TURBULENCE is not defined !!
#endif



//-------------------------------------------------------------------------------------------------------
// Structure   :  Turbulence_t
// Description :  Data structure of turbulence
//
// Data Member :  NMode            : Number of non-zero k modes
//                RSeed            : Current random seed value
//                Time             : Time to update turbulence pattern
//                Tdecay           : Turbulence correlation time
//                dt               : dt to update turbulence pattern
//                OUvar            : Ornstein-Uhlenbeck process variance
//                Kmode            : Non-zero k modes
//                Sin              : Array to store pre-computed sin modes
//                Cos              : Array to store pre-computed cos modes
//                Amplitude        : Amplitude of each k mode
//                OUphase          : Random phases updated by Ornstein-Uhlenbeck process
//
// Method      :  Turbulence_t        : Constructor
//               ~Turbulence_t        : Destructor
//-------------------------------------------------------------------------------------------------------
struct Turbulence_t
{

// data members
// ===================================================================================
   int     RSeed;
   long    NMode;
   double  Time;
   double  Tdecay;
   double  dt;
   double  OUvar;

   double *Kmode[3];
   double *Sin  [3];
   double *Cos  [3];
   double *OUphase;
   double *Amplitude;

   //===================================================================================
   // Constructor :  Particle_t
   // Description :  Constructor of the structure "Particle_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  None
   //===================================================================================
   Turbulence_t()
   {
      NMode    =  0;
      RSeed    =  1;
      Time     =  0.0;
      Tdecay   =  1.0;
      dt       =  1.0;
      OUvar    =  1.0;

      for (int i = 0; i < 3; i ++)
      {
         Kmode[i] = NULL;
         Sin  [i] = NULL;
         Cos  [i] = NULL;
      }
      Amplitude = NULL;
      OUphase   = NULL;

   } // METHOD : Turbulence_t



   //===================================================================================
   // Destructor  :  ~Turbulence_t
   // Description :  Destructor of the structure "Turbulence_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~Turbulence_t()
   {
      for (int i = 0; i < 3; i ++)
      {
         if ( Kmode[i] != NULL ) delete [] Kmode[i];
         if ( Sin  [i] != NULL ) delete [] Sin  [i];
         if ( Cos  [i] != NULL ) delete [] Cos  [i];
      }
      if ( Amplitude != NULL ) delete [] Amplitude;
      if ( OUphase   != NULL ) delete [] OUphase;


   } // METHOD : ~Turbulence_t


}; // struct Turbulence_t



#endif // #ifndef __TURBULENCE_H__
