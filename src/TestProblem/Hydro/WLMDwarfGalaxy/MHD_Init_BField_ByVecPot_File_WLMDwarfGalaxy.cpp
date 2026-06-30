#include "GAMER.h"

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static int nAx, nAy, nAz;
static double Axmin, Aymin, Azmin;
static double Adx, Ady, Adz;
static double *Axcoord, *Aycoord, *Azcoord;

#ifdef SUPPORT_HDF5
static void VecPot_ReadField( hid_t mag_file_id, const int ibegin, const int jbegin,
                              const int kbegin, double Ax[], double Ay[], double Az[] );
#endif

//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Init_BField_ByVecPot_File
// Description :  Use the input uniform-mesh array stored in the file "B_Filename"
//                to assign a magnetic vector potential to all real patches on level
//                "B_lv" and take the curl of this potential to compute the magnetic field
//
// Note        :
//
// Parameter   :  B_lv         : Target AMR level
//
// Return      :  amr->patch->magnetic
//-------------------------------------------------------------------------------------------------------
void MHD_Init_BField_ByVecPot_File_WLMDwarfGalaxy( const int B_lv )
{

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  else

#  ifndef SUPPORT_HDF5
   Aux_Error( ERROR_INFO, "SUPPORT_HDF5 must be set to load a vector potential from a file !!\n" );
#  endif
   char B_Filename[2*MAX_STRING];
   sprintf( B_Filename, "B_IC_lv%02d", B_lv );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading the magnetic field from the input file %s ...\n", B_Filename );

   const double dh = amr->dh[B_lv];

   double *Axf, *Ayf, *Azf;

   if ( !Aux_CheckFileExist(B_Filename) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", B_Filename );

// Open the magnetic field file and determine the dimensionality of the vector
// potential grid

#  ifdef SUPPORT_HDF5

   herr_t status;
   hid_t dataset, dataspace;
   hsize_t dims[3], maxdims[3];
   int ndim;

   hid_t mag_file_id = H5Fopen(B_Filename, H5F_ACC_RDONLY, H5P_DEFAULT);



   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_z", H5P_DEFAULT);

   dataspace = H5Dget_space(dataset);

   ndim = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

   if ( ndim != 3 ) Aux_Error( ERROR_INFO, "Incorrect dimensionality of vector potential ndim=%d !!\n", ndim );

   H5Sclose(dataspace);
   H5Dclose(dataset);

// NOTE: Magnetic vector potential arrays are stored in column-major order,
// i.e., Ax[nAx][nAy][nAz]
// A_z is cell-cented along z and edge-centered along x and y
// store the cell-centered nA here
   nAx = dims[0] - 1;
   nAy = dims[1] - 1;
   nAz = dims[2];

#  endif

// Read the edge-centered coordinate information from the vector potential grid
   Axcoord = new double [ nAx+1 ];
   Aycoord = new double [ nAy+1 ];
   Azcoord = new double [ nAz+1 ];

#  ifdef SUPPORT_HDF5

   dataset = H5Dopen(mag_file_id, "x", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Axcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load x-coordinate !!\n" );

   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "y", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Aycoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load y-coordinate !!\n" );
   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "z", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Azcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load z-coordinate !!\n" );
   H5Dclose(dataset);

#  endif

// Cell spacing and left edge of vector potential grid
// dx, dy, dz should be equal to dh
   Adx = Axcoord[1]-Axcoord[0];
   Ady = Aycoord[1]-Aycoord[0];
   Adz = Azcoord[1]-Azcoord[0];

   Axmin = Axcoord[0];
   Aymin = Aycoord[0];
   Azmin = Azcoord[0];

   double Axmax = Axcoord[nAx];
   double Aymax = Aycoord[nAy];
   double Azmax = Azcoord[nAz];

   if ( fabs(Adx-dh)/dh > 1e-10 || fabs(Ady-dh)/dh > 1e-10 || fabs(Adz-dh)/dh > 1e-10 )
       Aux_Error(ERROR_INFO, "   vector potential table spacing = [%13.7e, %13.7e, %13.7e] does not match AMR (dh = %13.7e)!!\n",
                                 Adx, Ady, Adz, dh);

   if (MPI_Rank == 0)
   Aux_Message(stderr, "      vector potential table %s:\n"
                       "         nA   = [ %13d, %13d, %13d ]\n"
                       "         Adh  = [ %13.7e, %13.7e, %13.7e ]\n"
                       "         Amin = [ %13.7e, %13.7e, %13.7e ]\n"
                       "         Amax = [ %13.7e, %13.7e, %13.7e ]\n",
                              B_Filename, Adx, Ady, Adz, nAx, nAy, nAz, Axmin, Aymin, Azmin, Axmax, Aymax, Azmax );

   double *Ax = new double [ CUBE(PS1+1) ];
   double *Ay = new double [ CUBE(PS1+1) ];
   double *Az = new double [ CUBE(PS1+1) ];

   for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++) {

      double EdgeL[3];
      double EdgeR[3];

      for (int i=0; i<3; i++) {
         EdgeL[i] = amr->patch[0][B_lv][PID]->EdgeL[i];
         EdgeR[i] = amr->patch[0][B_lv][PID]->EdgeR[i];
      }

// min(EdgeL[0]) should = Axmin, but add 0.1*dh to avoid rounding error
   if ( EdgeL[0] < Axmin-0.1*dh || EdgeR[0] > Axmax+0.1*dh ||
        EdgeL[1] < Aymin-0.1*dh || EdgeR[1] > Aymax+0.1*dh ||
        EdgeL[2] < Azmin-0.1*dh || EdgeR[2] > Azmax+0.1*dh )
      Aux_Error( ERROR_INFO, "Input patch lies outside the vector potential table!!\n"
                             "   MPI_Rank=%d, lv=%d PID=%d,\n"
                             "   EdgeL = [ %13.7e, %13.7e, %13.7e ], EdgeR = [ %13.7e, %13.7e, %13.7e ]\n"
                             "   Amin  = [ %13.7e, %13.7e, %13.7e ], Amax  = [ %13.7e, %13.7e, %13.7e ]\n",
                              MPI_Rank, B_lv, PID, EdgeL[0], EdgeL[1], EdgeL[2], EdgeR[0], EdgeR[1], EdgeR[2],
                              Axmin, Aymin, Azmin, Axmax, Aymax, Azmax );

      int ibegin = round((EdgeL[0]-Axmin)/dh);
      int jbegin = round((EdgeL[1]-Aymin)/dh);
      int kbegin = round((EdgeL[2]-Azmin)/dh);

//    Allocate for the data on the vector potential grid local to this patch and
//    read it from the file

      Axf = new double [PS1*(PS1+1)*(PS1+1)];
      Ayf = new double [(PS1+1)*PS1*(PS1+1)];
      Azf = new double [(PS1+1)*(PS1+1)*PS1];

#     ifdef SUPPORT_HDF5
      VecPot_ReadField( mag_file_id, ibegin, jbegin, kbegin,
                        Axf, Ayf, Azf );
#     endif

//    Loop over the indices in this patch and assign the vector potential
      for (int k=0; k<PS1+1; k++) {  const double z0 = EdgeL[2] + k*dh;
      for (int j=0; j<PS1+1; j++) {  const double y0 = EdgeL[1] + j*dh;
      for (int i=0; i<PS1+1; i++) {  const double x0 = EdgeL[0] + i*dh;

//       #define IDX321( i, j, k, Ni, Nj )   (  ( (k)*(Nj) + (j) )*(Ni) + (i)  )
         int idx = IDX321( i, j, k, PS1+1, PS1+1 );

         Ax[idx] = 0.0;
         Ay[idx] = 0.0;
         Az[idx] = 0.0;

         if ( i != PS1 ) {
            //const double x = x0 + (ii+0.5)*dh;
            const int idx2 = i*(PS1+1)*(PS1+1) + j*(PS1+1) + k;
            Ax[idx] = Axf[idx2];
         }
         if ( j != PS1 ) {
            //const double y = y0 + (jj+0.5)*dh;
            const int idx2 = i*(PS1+1)*PS1 + j*(PS1+1) + k;
            Ay[idx] = Ayf[idx2];
         }
         if ( k != PS1 ) {
            //const double z = z0 + (kk+0.5)*dh;
            const int idx2 = i*PS1*(PS1+1) + j*PS1 + k;
            Az[idx] = Azf[idx2];
         }
      }}}

//    Calculate Bx from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1+1; i++) {
         int idx  = IDX321   ( i, j,   k,   PS1+1, PS1+1 );
         int idxj = IDX321   ( i, j+1, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i, j,   k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BX( i, j,   k,   PS1,   PS1   );
         real Bx = ( Az[idxj] - Az[idx] - Ay[idxk] + Ay[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[0][idxB] = Bx;
      }}}

//    Calculate By from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1+1; j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j, k,   PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i,   j, k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BY( i,   j, k,   PS1,   PS1   );
         real By = ( Ax[idxk] - Ax[idx] - Az[idxi] + Az[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[1][idxB] = By;
      }}}

//    Calculate Bz from vector potential
      for (int k=0; k<PS1+1; k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j,   k, PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j,   k, PS1+1, PS1+1 );
         int idxj = IDX321   ( i,   j+1, k, PS1+1, PS1+1 );
         int idxB = IDX321_BZ( i,   j,   k, PS1,   PS1   );
         real Bz = ( Ay[idxi] - Ay[idx] - Ax[idxj] + Ax[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[2][idxB] = Bz;
      }}}

      delete [] Axf;
      delete [] Ayf;
      delete [] Azf;

   } // for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++)

// Close the magnetic field file

#  ifdef SUPPORT_HDF5
   H5Fclose(mag_file_id);
#  endif

   delete [] Ax;
   delete [] Ay;
   delete [] Az;

   delete [] Axcoord;
   delete [] Aycoord;
   delete [] Azcoord;

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   Loading the magnetic field from the input file ... done\n" );
#  endif // ifndef MHD ... else

} // FUNCTION : MHD_Init_BField_ByVecPot_File_WLMDwarfGalaxy

#ifdef SUPPORT_HDF5

void VecPot_ReadField( hid_t mag_file_id, const int ibegin, const int jbegin,
                       const int kbegin, double Ax[], double Ay[], double Az[] )
{
   hid_t dataset, dataspace, memspace, dxfer_template;

   herr_t status;

   hsize_t start[3], stride[3], dimsx[3], dimsy[3], dimsz[3];

   int rank, ierr;

   rank = 3;

   start[0] = ibegin;
   start[1] = jbegin;
   start[2] = kbegin;

   stride[0] = 1;
   stride[1] = 1;
   stride[2] = 1;

   dimsx[0] = PS1;
   dimsx[1] = PS1+1;
   dimsx[2] = PS1+1;

   dimsy[0] = PS1+1;
   dimsy[1] = PS1;
   dimsy[2] = PS1+1;

   dimsz[0] = PS1+1;
   dimsz[1] = PS1+1;
   dimsz[2] = PS1;

// Read Ax
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_x", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, dimsx, NULL);
   memspace = H5Screate_simple(rank, dimsx, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, Ax);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_x !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Ay
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_y", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, dimsy, NULL);
   memspace = H5Screate_simple(rank, dimsy, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, Ay);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_y !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Az
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_z", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, dimsz, NULL);
   memspace = H5Screate_simple(rank, dimsz, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, Az);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_z !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   return;

} // FUNCTION : VecPot_ReadField

#endif // #ifdef SUPPORT_HDF5
