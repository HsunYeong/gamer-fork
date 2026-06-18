import h5py
import numpy as np
import yt

# Vector potential parameters
r0     = 30.0    # in kpc
rho0   = 1.0     # in code_density
B0     = 1.0e-9  # in gauss
Rho_bg = 1.0e-36 # background gas to add if the data on the grid is zero, in g/cm^3

# This part should be consistent with construct_ic.py!!
# ---------------------------------------------------------------------------------------------------------------------
# GIZMO data information
GIZMO_Filename                  = 'snap_000.hdf5' # filename of the input GIZMO data
GIZMO_UnitLength_in_cm          = 3.085678e+21    # unit of time in GIZMO, in cm/h, where h = HubbleParam in GIZMO
GIZMO_UnitMass_in_g             = 1.989e+43       # unit of mass in GIZMO, in g/h, where h = HubbleParam in GIZMO
GIZMO_UnitVelocity_in_cm_per_s  = 100000          # unit of velocity in GIZMO, in cm/s

unit_base = {
              'UnitLength_in_cm'        : GIZMO_UnitLength_in_cm,
              'UnitMass_in_g'           : GIZMO_UnitMass_in_g,
              'UnitVelocity_in_cm_per_s': GIZMO_UnitVelocity_in_cm_per_s,
            }

# Set grid information
Target_N_x_base = 256   # number of Lv0 cells along x-direction of UM_IC; actual number may be larger to satisfy the proper nesting of AMR
N_sub           = 1     # number of cells for sub-sampling
# Width of region in each refined level     [   Lv0,  Lv1,  Lv2,  Lv3,  Lv4,  Lv5,  Lv6, Lv7, Lv8 ]
Target_Refined_Width_x  = np.array([450.0, 16.0, 16.0, 16.0, 16.0, 14.0, 11.2, 6.4, 2.8])
Target_Refined_Width_y  = np.array([450.0, 16.0, 16.0, 16.0, 16.0, 14.0, 11.2, 6.4, 2.8])
Target_Refined_Width_z  = np.array([450.0,  8.0,  8.0,  8.0,  8.0,  7.0,  5.6, 3.2, 1.4])
PatchSize               = 8
factorMassNormalization = 7.08799838 / 9.40795185 # factor for rescaling the density to ensure total mass consistent
# ---------------------------------------------------------------------------------------------------------------------

# Check
assert Target_N_x_base % (4*PatchSize) == 0,                              'Target_N_x_base should be multiples of 4*PatchSize'
assert len(Target_Refined_Width_x) > 0,                                   'Length of Target_Refined_Width_x should be larger than 0'
assert len(Target_Refined_Width_y) == len(Target_Refined_Width_x),        'Given length of Target_Refined_Width_y is inconsistent with Target_Refined_Width_x'
assert len(Target_Refined_Width_z) == len(Target_Refined_Width_x),        'Given length of Target_Refined_Width_z is inconsistent with Target_Refined_Width_x'
assert np.all(Target_Refined_Width_x[:-1] >= Target_Refined_Width_x[1:]), 'Target_Refined_Width_x should be monotonically decreasing'
assert np.all(Target_Refined_Width_y[:-1] >= Target_Refined_Width_y[1:]), 'Target_Refined_Width_y should be monotonically decreasing'
assert np.all(Target_Refined_Width_z[:-1] >= Target_Refined_Width_z[1:]), 'Target_Refined_Width_z should be monotonically decreasing'


# Set other parameters
NLEVEL = len(Target_Refined_Width_x)
dh     = np.array( [ Target_Refined_Width_x[0]/(Target_N_x_base*2**lv) for lv in range( 0, NLEVEL, 1 ) ] )
N_x    = np.array( [ 4*PatchSize*np.ceil( Target_Refined_Width_x[lv]/( 4*PatchSize*dh[lv] ) ) for lv in range( 0, NLEVEL, 1 ) ], dtype=np.uint32 )
N_y    = np.array( [ 4*PatchSize*np.ceil( Target_Refined_Width_y[lv]/( 4*PatchSize*dh[lv] ) ) for lv in range( 0, NLEVEL, 1 ) ], dtype=np.uint32 )
N_z    = np.array( [ 4*PatchSize*np.ceil( Target_Refined_Width_z[lv]/( 4*PatchSize*dh[lv] ) ) for lv in range( 0, NLEVEL, 1 ) ], dtype=np.uint32 )
UNIT_P = GIZMO_UnitMass_in_g * GIZMO_UnitVelocity_in_cm_per_s**2.0 / GIZMO_UnitLength_in_cm**3.0 # energy density unit
UNIT_B = ( UNIT_P*4.0*np.pi )**0.5                                                               # magnetic field unit

for lv in range( NLEVEL-2, -1, -1 ):
    N_x[lv] = max( N_x[lv], N_x[lv+1]//2 + 2*PatchSize + 2*PatchSize*int( (N_x[lv+1]//2)%(4*PatchSize) == 0 ) ) # Make sure the boundaries of range lv and lv+1 are separated by at least 1 patch of lv
    N_y[lv] = max( N_y[lv], N_y[lv+1]//2 + 2*PatchSize + 2*PatchSize*int( (N_y[lv+1]//2)%(4*PatchSize) == 0 ) )
    N_z[lv] = max( N_z[lv], N_z[lv+1]//2 + 2*PatchSize + 2*PatchSize*int( (N_z[lv+1]//2)%(4*PatchSize) == 0 ) )

width_x = N_x*dh
width_y = N_y*dh
width_z = N_z*dh

EdgeL = -0.5*np.array([width_x, width_y, width_z])
EdgeR =  0.5*np.array([width_x, width_y, width_z])

np.set_printoptions( formatter={'float':'{: 16.8e}'.format, 'int':'{: 16d}'.format}, linewidth=1000 )
print( '' )
print( '------------------------------------------------------------------------------------------------' )
print( 'B_IC information' )
print( '------------------------------------------------------------------------------------------------' )
print( f'{PatchSize         = }' )
print( f'{NLEVEL            = }' )
print( f'{N_x               = }' )
print( f'{N_y               = }' )
print( f'{N_z               = }' )
print( '' )
print( f'{width_x           = }' )
print( f'{width_y           = }' )
print( f'{width_z           = }' )
print( '' )
print( f'{dh                = }' )
print( '' )
print( f'{EdgeL[0]             = }' )
print( f'{EdgeR[0]             = }' )
print( '' )
print( f'{EdgeL[1]             = }' )
print( f'{EdgeR[1]             = }' )
print( '' )
print( f'{EdgeL[2]             = }' )
print( f'{EdgeR[2]             = }' )
print( '------------------------------------------------------------------------------------------------' )
print( '' )

bbox = [ [EdgeL[0,0] - 0.5*dh[0], EdgeR[0,0] + 0.5*dh[0]],
         [EdgeL[1,0] - 0.5*dh[0], EdgeR[1,0] + 0.5*dh[0]],
         [EdgeL[2,0] - 0.5*dh[0], EdgeR[2,0] + 0.5*dh[0]] ]

ds = yt.load( GIZMO_Filename, unit_base=unit_base, bounding_box=bbox )

fine_data = {}

for lv in range(NLEVEL-1, -1, -1):

   print(f'Processing level {lv}')
   x = np.linspace( EdgeL[0,lv], EdgeR[0,lv], N_x[lv] + 1 )
   y = np.linspace( EdgeL[1,lv], EdgeR[1,lv], N_y[lv] + 1 )
   z = np.linspace( EdgeL[2,lv], EdgeR[2,lv], N_z[lv] + 1 )

#  Find the grid cell midpoints
   x_cell = 0.5*(x[1:] + x[:-1])
   y_cell = 0.5*(y[1:] + y[:-1])
   z_cell = 0.5*(z[1:] + z[:-1])

   f = h5py.File(f'B_IC_lv{lv:02d}', "w")

#  Write coordinate arrays
   f.create_dataset("x", data=x+0.5*width_x[0])
   f.create_dataset("y", data=y+0.5*width_y[0])
   f.create_dataset("z", data=z+0.5*width_z[0])

#  loop through 3 directions
   for i in range(3):
      dh_arr    = np.array([dh[lv], dh[lv], dh[lv]])
      dh_arr[i] = 0
      dims      = np.array([N_x[lv], N_y[lv], N_z[lv]]) + 1
      dims[i]  -= 1
      ag = ds.arbitrary_grid( EdgeL[:,lv] - 0.5*dh_arr, EdgeR[:,lv] + 0.5*dh_arr, dims=dims )
      ag_density = ag['PartType0', 'density']
      dens = ( 1.0/N_sub**3 )*factorMassNormalization*ag_density
      dens[ dens <= 0 ] += ds.quan( Rho_bg, 'g/cm**3')
      dens = dens.in_units('code_density').d

#     Use the 1-D coordinate arrays to consruct 3D coordinate arrays
#     that we will use to compute an analytic vector potential
      xx, yy, zz = np.meshgrid( x_cell, y, z, sparse=False, indexing='ij' ) if i == 0 else \
                   np.meshgrid( x, y_cell, z, sparse=False, indexing='ij' ) if i == 1 else \
                   np.meshgrid( x, y, z_cell, sparse=False, indexing='ij' )

#     Assign vector potential
      r   =  np.sqrt(xx**2 + yy**2 + zz**2)
      amp =  np.where( r < r0, 1.0, np.exp( -(r-r0)/r0) )
      A   =  yy*B0*amp*(dens/rho0)**(2.0/3.0) / UNIT_B if i == 0 else \
            -xx*B0*amp*(dens/rho0)**(2.0/3.0) / UNIT_B if i == 1 else \
             np.zeros_like(zz)                / UNIT_B

#     Restrict data
      if lv < NLEVEL-1:
#        Get fine grid data
         A_f = fine_data['A%s'%(chr(ord('x') + i))]
         assert A_f.shape[ i     ] % 2 == 0
         assert A_f.shape[(i+1)%3] % 2 == 1
         assert A_f.shape[(i+2)%3] % 2 == 1

#        Take average for every two data
         if i == 0:
            A_r = 0.5*( A_f[0::2,:,:] + A_f[1::2,:,:] )
            A_r = A_r[:,0::2,0::2]
         elif i == 1:
            A_r = 0.5*( A_f[:,0::2,:] + A_f[:,1::2,:] )
            A_r = A_r[0::2,:,0::2]
         else:
            A_r = 0.5*( A_f[:,:,0::2] + A_f[:,:,1::2] )
            A_r = A_r[0::2,0::2,:]

#        Overwrite the corresponding coarse cells
         nx, ny, nz = A_r.shape
         i0 = int(round((EdgeL[0,lv+1]-EdgeL[0,lv])/dh[lv]))
         j0 = int(round((EdgeL[1,lv+1]-EdgeL[1,lv])/dh[lv]))
         k0 = int(round((EdgeL[2,lv+1]-EdgeL[2,lv])/dh[lv]))

         A[i0:i0+nx, j0:j0+ny, k0:k0+nz] = A_r
      f.create_dataset("magnetic_vector_potential_%s"%(chr(ord('x') + i)), data=A)
      fine_data['A%s'%(chr(ord('x') + i))] = A.copy()

#  Close the file
   f.flush()
   f.close()


