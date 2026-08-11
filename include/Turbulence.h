#ifndef __TURBULENCE_H__
#define __TURBULENCE_H__



#include "Macro.h"



//-------------------------------------------------------------------------------------------------------
// Structure   :  Turbulence_t
// Description :  Data structure of turbulence
//
// Data Member :  NMode            : Number of non-zero k modes
//                RSeed            : Current random seed value
//                Tdecay           : Turbulence correlation time
//                dt               : dt to update turbulence pattern
//                OUvar            : Ornstein-Uhlenbeck process variance
//                TimeLast         : Time for last turbulence field
//                TimeNext         : Time for next turbulence field
//                IdxLast          : Index for last turbulence field
//                IdxNext          : Index for next turbulence field
//                Amplitude        : Amplitude of each k mode
//                Kmode            : Non-zero k modes
//                Sin              : Array to store pre-computed sin modes
//                Cos              : Array to store pre-computed cos modes
//                OUphase          : Random phases updated by Ornstein-Uhlenbeck process, two entries for TimeLast and TimeNext
//
// Method      :  Turbulence_t     : Constructor
//               ~Turbulence_t     : Destructor
//-------------------------------------------------------------------------------------------------------
struct Turbulence_t
{

// data members
// ===================================================================================
   int     RSeed;
   long    NMode;
   double  Tdecay;
   double  dt;
   double  OUvar;
   double  TimeLast;
   double  TimeNext;
   int     IdxLast;
   int     IdxNext;

   double *Amplitude;
   double *Kmode[3];
   double *Sin  [3];
   double *Cos  [3];
   double *OUphase[2];

   //===================================================================================
   // Constructor :  Turbulence_t
   // Description :  Constructor of the structure "Turbulence_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  None
   //===================================================================================
   Turbulence_t()
   {
      NMode    =  0;
      RSeed    =  1;
      Tdecay   =  1.0;
      dt       =  1.0;
      OUvar    =  1.0;
      TimeLast = -1.0;
      TimeNext = -1.0;
      IdxLast  =  0;
      IdxNext  =  1;

      Amplitude  = NULL;

      for (int i = 0; i < 3; i ++)
      {
         Kmode[i] = NULL;
         Sin  [i] = NULL;
         Cos  [i] = NULL;
      }

      for (int i = 0; i < 2; i ++)
         OUphase[i] = NULL;

   } // METHOD : Turbulence_t



   //===================================================================================
   // Destructor  :  ~Turbulence_t
   // Description :  Destructor of the structure "Turbulence_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~Turbulence_t()
   {
      if ( Amplitude != NULL ) delete [] Amplitude;

      for (int i = 0; i < 3; i ++)
      {
         if ( Kmode[i] != NULL ) delete [] Kmode[i];
         if ( Sin  [i] != NULL ) delete [] Sin  [i];
         if ( Cos  [i] != NULL ) delete [] Cos  [i];
      }

      for (int i = 0; i < 2; i ++)
         if ( OUphase[i] != NULL ) delete [] OUphase[i];

   } // METHOD : ~Turbulence_t

}; // struct Turbulence_t



#endif // #ifndef __TURBULENCE_H__
