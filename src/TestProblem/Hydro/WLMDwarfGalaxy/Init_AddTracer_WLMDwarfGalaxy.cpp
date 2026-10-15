#include "GAMER.h"

#if ( defined TRACER && defined SUPPORT_GSL )
#include <gsl/gsl_qrng.h>

extern bool   WLMDwarfGalaxy_AddTracer;
extern int    WLMDwarfGalaxy_AddTracerMinLv;
extern double WLMDwarfGalaxy_TracerMass;


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_AddTracer_WLMDwarfGalaxy()
// Description :  Add a new disk from an existing snapshot
//
// Note        :  Must enable OPT__RESTART_RESET and AddParWhenRestart
//-------------------------------------------------------------------------------------------------------
void Init_AddTracer_WLMDwarfGalaxy()
{

   if ( !WLMDwarfGalaxy_AddTracer )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, 3);

   for (int lv=WLMDwarfGalaxy_AddTracerMinLv; lv<=MAX_LEVEL; lv++)
   {
      const double dh             = amr->dh[lv];
      const real   dv             = CUBE( dh );
      const int    FluSg          = amr->FluSg[lv];

      const int    MaxNewParPerPatch = 8*CUBE(PS1);
      real_par   (*NewParAttFlt)[PAR_NATT_FLT_TOTAL] = new real_par [MaxNewParPerPatch][PAR_NATT_FLT_TOTAL];
      long_par   (*NewParAttInt)[PAR_NATT_INT_TOTAL] = new long_par [MaxNewParPerPatch][PAR_NATT_INT_TOTAL];
      long        *NewParID                          = new long     [MaxNewParPerPatch];


#     pragma omp for schedule( static )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         double x0   = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         double y0   = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         double z0   = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;
         int NNewPar = 0;

         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            if ( NNewPar >= MaxNewParPerPatch )
               Aux_Error( ERROR_INFO, "NNewPar (%d) >= MaxNewParPerPatch (%d) !!\n", NNewPar, MaxNewParPerPatch );

            double GasMass = amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i] * dv;

            double x = x0 + i*dh;
            double y = y0 + j*dh;
            double z = z0 + k*dh;

//          determine number of tracers to be added
            int NTracer = (lv < MAX_LEVEL)? 1 : MAX( 1, std::round( GasMass/WLMDwarfGalaxy_TracerMass ) );
            NTracer = MIN( 8, NTracer );

            if ( NTracer == 8 )
            {
               for (int kc=0; kc<2; kc++)
               for (int jc=0; jc<2; jc++)
               for (int ic=0; ic<2; ic++)
               {
                  NewParAttFlt[NNewPar][PAR_MASS] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_POSX] = x - 0.25*dh + 0.5*dh*ic;
                  NewParAttFlt[NNewPar][PAR_POSY] = y - 0.25*dh + 0.5*dh*jc;
                  NewParAttFlt[NNewPar][PAR_POSZ] = z - 0.25*dh + 0.5*dh*kc;
                  NewParAttFlt[NNewPar][PAR_VELX] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_VELY] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_VELZ] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_TIME] = (real_par)Time[0];
                  NewParAttInt[NNewPar][PAR_TYPE] = PTYPE_TRACER;
                  NewParAttInt[NNewPar][PAR_PUID] = PUID_TBA;
                  NewParAttInt[NNewPar][PAR_FLAG] = PFLAG_NO;
#                 ifdef STAR_FORMATION
                  NewParAttFlt[NNewPar][Idx_ParCreTime  ] = (real_par)0.0;
#                 endif
#                 ifdef FEEDBACK
                  NewParAttFlt[NNewPar][Idx_ParSNIITime ] = (real_par)0.0;
#                 endif
                  NewParAttFlt[NNewPar][Idx_ParMetalFrac] = (real_par)0.0;

                  NNewPar ++;
               }
            }
//          add one particle at center
            else if ( NTracer == 1 )
            {
               NewParAttFlt[NNewPar][PAR_MASS] = (real_par)0.0;
               NewParAttFlt[NNewPar][PAR_POSX] = x;
               NewParAttFlt[NNewPar][PAR_POSY] = y;
               NewParAttFlt[NNewPar][PAR_POSZ] = z;
               NewParAttFlt[NNewPar][PAR_VELX] = (real_par)0.0;
               NewParAttFlt[NNewPar][PAR_VELY] = (real_par)0.0;
               NewParAttFlt[NNewPar][PAR_VELZ] = (real_par)0.0;
               NewParAttFlt[NNewPar][PAR_TIME] = (real_par)Time[0];
               NewParAttInt[NNewPar][PAR_TYPE] = PTYPE_TRACER;
               NewParAttInt[NNewPar][PAR_PUID] = PUID_TBA;
               NewParAttInt[NNewPar][PAR_FLAG] = PFLAG_NO;
#              ifdef STAR_FORMATION
               NewParAttFlt[NNewPar][Idx_ParCreTime  ] = (real_par)0.0;
#              endif
#              ifdef FEEDBACK
               NewParAttFlt[NNewPar][Idx_ParSNIITime ] = (real_par)0.0;
#              endif
               NewParAttFlt[NNewPar][Idx_ParMetalFrac] = (real_par)0.0;

               NNewPar ++;
            }
            else
            {
               for (int nt=0; nt<NTracer; nt++)
               {
                  double pos[3];
                  gsl_qrng_get (q, pos);

                  NewParAttFlt[NNewPar][PAR_MASS] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_POSX] = x - 0.5*dh + pos[0]*dh;
                  NewParAttFlt[NNewPar][PAR_POSY] = y - 0.5*dh + pos[1]*dh;
                  NewParAttFlt[NNewPar][PAR_POSZ] = z - 0.5*dh + pos[2]*dh;
                  NewParAttFlt[NNewPar][PAR_VELX] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_VELY] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_VELZ] = (real_par)0.0;
                  NewParAttFlt[NNewPar][PAR_TIME] = (real_par)Time[0];
                  NewParAttInt[NNewPar][PAR_TYPE] = PTYPE_TRACER;
                  NewParAttInt[NNewPar][PAR_PUID] = PUID_TBA;
                  NewParAttInt[NNewPar][PAR_FLAG] = PFLAG_NO;
#                 ifdef STAR_FORMATION
                  NewParAttFlt[NNewPar][Idx_ParCreTime  ] = (real_par)0.0;
#                 endif
#                 ifdef FEEDBACK
                  NewParAttFlt[NNewPar][Idx_ParSNIITime ] = (real_par)0.0;
#                 endif
                  NewParAttFlt[NNewPar][Idx_ParMetalFrac] = (real_par)0.0;

                  NNewPar ++;
               }
            }

         } // for k,j,i

#        pragma omp critical
         {
//          add particles to the particle repository
            for (int p=0; p<NNewPar; p++)
               NewParID[p] = amr->Par->AddOneParticle( NewParAttFlt[p], NewParAttInt[p] );

//          add particles to the patch
            const long_par *PType = amr->Par->Type;

            amr->patch[0][lv][PID]->AddParticle( NNewPar, NewParID, &amr->Par->NPar_Lv[lv], PType );

         } // pragma omp critical

      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

//    free memory
      delete [] NewParAttFlt;
      delete [] NewParAttInt;
      delete [] NewParID;

      MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

   } // for lv

   gsl_qrng_free (q);

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_AddTracer_WLMDwarfGalaxy
#endif // #if ( defined TRACER && defined SUPPORT_GSL )


