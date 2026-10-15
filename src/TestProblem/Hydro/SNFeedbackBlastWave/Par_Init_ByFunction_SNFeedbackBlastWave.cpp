#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

extern double SNBlast_ParMass;
extern double SNBlast_ParCenter[3];
extern double SNBlast_ParVelocity[3];
extern double SNBlast_ParMetalMassFrac;
#ifdef TRACER
extern double SNBlast_Use_Tracers;
extern int    SNBlast_NPar[3];
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_SNFeedbackBlastWave
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank   : Number of particles to be set by this MPI rank
//                NPar_AllRank    : Total Number of particles in all MPI ranks
//                ParMass         : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z     : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z     : Particle velocity array with the size of NPar_ThisRank
//                ParTime         : Particle time     array with the size of NPar_ThisRank
//                ParType         : Particle type     array with the size of NPar_ThisRank
//                AllAttributeFlt : Pointer array for all particle floating-point attributes
//                                  --> Dimension = [PAR_NATT_FLT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                      to access the data
//                AllAttributeInt : Pointer array for all particle integer attributes
//                                  --> Dimension = [PAR_NATT_INT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttributeFlt, AllAttributeInt
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_SNFeedbackBlastWave( const long NPar_ThisRank, const long NPar_AllRank,
                                              real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                              real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                              long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                              long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   long NPar_All = 1;
   if ( SNBlast_Use_Tracers ) NPar_All += SNBlast_NPar[0]*SNBlast_NPar[1]*SNBlast_NPar[2];

   if ( NPar_All != NPar_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles found [%ld] != expect [%ld] !!\n",
                 NPar_All, NPar_AllRank );

// define the particle attribute arrays
   real_par *ParFltData_AllRank[PAR_NATT_FLT_TOTAL];
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   ParFltData_AllRank[v] = NULL;
   long_par *ParIntData_AllRank[PAR_NATT_INT_TOTAL];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   ParIntData_AllRank[v] = NULL;

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 ) {

//    allocate memory for particle attribute arrays
      ParFltData_AllRank[PAR_MASS] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSZ] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELZ] = new real_par [NPar_AllRank];
      ParIntData_AllRank[PAR_TYPE] = new long_par [NPar_AllRank];
      ParIntData_AllRank[PAR_PUID] = new long_par [NPar_AllRank];
#     ifdef FEEDBACK
      ParFltData_AllRank[Idx_ParSNIITime ] = new real_par [NPar_AllRank];
#     endif
      ParFltData_AllRank[Idx_ParMetalFrac] = new real_par [NPar_AllRank];

      long p = 0;

      ParFltData_AllRank[PAR_MASS][p] = (real_par)SNBlast_ParMass;
      ParFltData_AllRank[PAR_POSX][p] = (real_par)SNBlast_ParCenter[0];
      ParFltData_AllRank[PAR_POSY][p] = (real_par)SNBlast_ParCenter[1];
      ParFltData_AllRank[PAR_POSZ][p] = (real_par)SNBlast_ParCenter[2];
      ParFltData_AllRank[PAR_VELX][p] = (real_par)SNBlast_ParVelocity[0];
      ParFltData_AllRank[PAR_VELY][p] = (real_par)SNBlast_ParVelocity[1];
      ParFltData_AllRank[PAR_VELZ][p] = (real_par)SNBlast_ParVelocity[2];
      ParIntData_AllRank[PAR_TYPE][p] = PTYPE_STAR;
      ParIntData_AllRank[PAR_PUID][p] = PUID_TBA;
#     ifdef FEEDBACK
      ParFltData_AllRank[Idx_ParSNIITime ][p] = (real_par)Time[0] + (real_par)FB_RESOLVED_SNEII_DELAY_TIME;
#     endif
      ParFltData_AllRank[Idx_ParMetalFrac][p] = (real_par)SNBlast_ParMetalMassFrac;

      p++;

#     ifdef TRACER
      if ( SNBlast_Use_Tracers ) {

         const double delta_p[3] = { amr->BoxSize[0]/(SNBlast_NPar[0]+1),
                                     amr->BoxSize[1]/(SNBlast_NPar[1]+1),
                                     amr->BoxSize[2]/(SNBlast_NPar[2]+1) };

         for (long kk=0; kk<SNBlast_NPar[2]; kk++)
         for (long jj=0; jj<SNBlast_NPar[1]; jj++)
         for (long ii=0; ii<SNBlast_NPar[0]; ii++)
         {

//          tracer particles have no mass
            ParFltData_AllRank[PAR_MASS][p] = (real_par)0.0;

//          assign positions
            ParFltData_AllRank[PAR_POSX][p] = real_par( (ii+1)*delta_p[0] );
            ParFltData_AllRank[PAR_POSY][p] = real_par( (jj+1)*delta_p[1] );
            ParFltData_AllRank[PAR_POSZ][p] = real_par( (kk+1)*delta_p[2] );

//          set velocities to zero (these will be updated from the grid later)
            ParFltData_AllRank[PAR_VELX][p] = (real_par)0.0;
            ParFltData_AllRank[PAR_VELY][p] = (real_par)0.0;
            ParFltData_AllRank[PAR_VELZ][p] = (real_par)0.0;

//          set the particle type to be tracer
            ParIntData_AllRank[PAR_TYPE][p] = PTYPE_TRACER;
            ParIntData_AllRank[PAR_PUID][p] = PUID_TBA;
#           ifdef FEEDBACK
            ParFltData_AllRank[Idx_ParSNIITime ][p] = (real_par)0.0;
#           endif
            ParFltData_AllRank[Idx_ParMetalFrac][p] = (real_par)0.0;


            p++;
         }
      } // if ( ParTest_Use_Tracers )
#     endif
   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL|
#                           ifdef FEEDBACK
                            BIDX(Idx_ParSNIITime)|
#                           endif
                            BIDX(Idx_ParMetalFrac), _PAR_TYPE|_PAR_PUID,
                            ParFltData_AllRank, ParIntData_AllRank, AllAttributeFlt, AllAttributeInt );

// synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)
      ParTime[p] = (real_par)Time[0];

// free resource
   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   delete [] ParFltData_AllRank[v];
      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   delete [] ParIntData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_SNFeedbackBlastWave



#endif // #ifdef MASSIVE_PARTICLES
