import yt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import gc



######################################################################################################
# README
######################################################################################################
# - To convert the GIZMO snapshot into the GAMER initial conditions file format (PAR_IC and UM_IC)
#    - PAR_IC for the PartType1 (halo) and PartType2 (disk) particles
#       - PAR_INIT          =  3   (ByFile)
#       - PAR_IC_FORMAT     =  1   ([attribute][id])
#       - PAR_IC_MASS       = -1.0 (off)
#       - PAR_IC_TYPE       = -1   (off)
#    - UM_IC for the PartType0 (gas) particles
#       - OPT__INIT         =  3   (ByFile)
#       - OPT__UM_IC_NVAR   = -1   (Hydro=5+passive)
#       - OPT__UM_IC_FORMAT =  1   (vzyx)
# - The units are the same as the input GIZMO snapshot
# - Assume the region center of each level is the same = the domain center of the GIZMO snapshot
#    - The region width of each level is the only degree of freedom
# - Number of cells of each level in each direction must be a multiples of 4*PatchSize
#    - The left-side and the rigth-side are symmetric
#    - Must have integer multiples of patch group (2*PatchSize) in each side
######################################################################################################

######################################################################################################
# Table of contents
######################################################################################################
# 0.0 Set parameters
# 0.1 Set input GIZMO data information
# 0.2 Set output PAR_IC information
# 0.3 Set output UM_IC information
# 0.4 Check
# 0.5 Set other parameters
#
# 1.0 Load Data
# 1.1 Print GIZMO data information
#
# 2.0 Particle information
# 2.1 Plot the projections for particles in GIZMO data
#
# 3.0 Prepare the particle data for PartType1 and PartType2
# 3.1 Plot the PartType1 and PartType2 particles
# 3.2 Construct PAR_IC
# 3.3 Output PAR_IC information
#
# 4.0 Prepare AMR grid data for the PartType0 particles
# 4.1 Create Arbitrary Grid
# 4.2 Save the arbitrary grid object as dataset
# 4.3 Convert to conserved variables in GAMER
# 4.4 Add the background to the density and energy density
# 4.5 Plot the results of grid data
# 4.6 Construct UM_IC
# 4.7 Create Input__UM_IC_RefineRegion
# 4.8 Output UM_IC information
#
# 5.0 Output GAMER Input__Parameter information
######################################################################################################



# 0.0 Set parameters
VERBOSE                         = True            # whether to print detailed information


# 0.1 Set input GIZMO data information
GIZMO_Filename                  = 'snap_000.hdf5' # filename of the input GIZMO data
GIZMO_UnitLength_in_cm          = 3.085678e+21    # unit of time in GIZMO, in cm/h, where h = HubbleParam in GIZMO
GIZMO_UnitMass_in_g             = 1.989e+43       # unit of mass in GIZMO, in g/h, where h = HubbleParam in GIZMO
GIZMO_UnitVelocity_in_cm_per_s  = 100000          # unit of velocity in GIZMO, in cm/s

# 0.1.1 Check the virial 200 quantities in MakeDiskGalaxy
MakeDiskGalaxy_R200             = 45.0            # R200, in kpc


# 0.2 Set output PAR_IC information
PAR_IC_FLOAT8                   = False           # whether the PAR_IC floating-pointer number is in double precision
PAR_IC_INT8                     = False           # whether the PAR_IC integer is in double precision
OutputParCreTime                = True            # whether to output particle creation time as attribute, needed when GAMER has star formation
OutputParSNIITime               = True            # whether to output particle SNII explosion time as attribute, needed when GAMER has SNe feedback
OutputParMetalFrac              = True            # whether to output particle metal fraction as attribute, needed when GAMER has metal



# 0.3 Set output UM_IC information
Target_N_x_base                 = 256             # number of Lv0 cells along x-direction of UM_IC; actual number may be larger to satisfy the proper nesting of AMR
# Width of region in each refined level     [   Lv0,  Lv1,  Lv2,  Lv3,  Lv4,  Lv5,  Lv6, Lv7, Lv8 ]
Target_Refined_Width_x          = np.array( [ 450.0, 16.0, 16.0, 16.0, 16.0, 14.0, 11.2, 6.4, 2.8 ] ) # in GIZMO code_length
Target_Refined_Width_y          = np.array( [ 450.0, 16.0, 16.0, 16.0, 16.0, 14.0, 11.2, 6.4, 2.8 ] ) # actual width of the results may be larger to have integer multiples of patche groups
Target_Refined_Width_z          = np.array( [ 450.0,  8.0,  8.0,  8.0,  8.0,  7.0,  5.6, 3.2, 1.4 ] )
PatchSize                       = 8               # PATCH_SIZE in GAMER
UM_IC_FLOAT8                    = False           # whether the UM_IC is in double precision
Rho_bg                          = 1.0e-36         # background gas to add if the data on the grid is zero, in g/cm^3
OutputMetal                     = True            # whether to output the metal density as a passive field
MetalMassFrac                   = 0.002           # metal mass fraction
OutputDust                      = True            # whether to output the dust density as a passive field
DustMassFrac                    = 0.001           # dust mass fraction
N_sub                           = 1               # number of cells for sub-sampling
factorMassNormalization         = 7.08799838 / 9.40795185 # factor for rescaling the density to ensure total mass consistent


# 0.4 Check
assert Target_N_x_base % (4*PatchSize) == 0,                              'Target_N_x_base should be multiples of 4*PatchSize'
assert len(Target_Refined_Width_x) > 0,                                   'Length of Target_Refined_Width_x should be larger than 0'
assert len(Target_Refined_Width_y) == len(Target_Refined_Width_x),        'Given length of Target_Refined_Width_y is inconsistent with Target_Refined_Width_x'
assert len(Target_Refined_Width_z) == len(Target_Refined_Width_x),        'Given length of Target_Refined_Width_z is inconsistent with Target_Refined_Width_x'
assert np.all(Target_Refined_Width_x[:-1] >= Target_Refined_Width_x[1:]), 'Target_Refined_Width_x should be monotonically decreasing'
assert np.all(Target_Refined_Width_y[:-1] >= Target_Refined_Width_y[1:]), 'Target_Refined_Width_y should be monotonically decreasing'
assert np.all(Target_Refined_Width_z[:-1] >= Target_Refined_Width_z[1:]), 'Target_Refined_Width_z should be monotonically decreasing'


# 0.5 Set other parameters
dtype_float_PAR_IC  = np.double if PAR_IC_FLOAT8 else np.single
dtype_int_PAR_IC    = np.int64  if PAR_IC_INT8   else np.int32
dtype_UM_IC         = np.double if UM_IC_FLOAT8  else np.single

UM_IC_NLEVEL        = len(Target_Refined_Width_x)         # number of levels of UM_IC
UM_IC_dh            = np.array( [ Target_Refined_Width_x[0]/(Target_N_x_base*2**lv) for lv in range( 0, UM_IC_NLEVEL, 1 ) ] )                                               # cell size of the UM_IC for each level
UM_IC_N_x           = np.array( [ 4*PatchSize*np.ceil( Target_Refined_Width_x[lv]/( 4*PatchSize*UM_IC_dh[lv] ) ) for lv in range( 0, UM_IC_NLEVEL, 1 ) ], dtype=np.uint32 ) # number of cells on each level in the x direction
UM_IC_N_y           = np.array( [ 4*PatchSize*np.ceil( Target_Refined_Width_y[lv]/( 4*PatchSize*UM_IC_dh[lv] ) ) for lv in range( 0, UM_IC_NLEVEL, 1 ) ], dtype=np.uint32 ) # ...                                  y direction
UM_IC_N_z           = np.array( [ 4*PatchSize*np.ceil( Target_Refined_Width_z[lv]/( 4*PatchSize*UM_IC_dh[lv] ) ) for lv in range( 0, UM_IC_NLEVEL, 1 ) ], dtype=np.uint32 ) # ...                                  z direction

for lv in range( UM_IC_NLEVEL-2, -1, -1 ):
    UM_IC_N_x[lv] = max( UM_IC_N_x[lv],  UM_IC_N_x[lv+1]//2 + 2*PatchSize + 2*PatchSize*int( (UM_IC_N_x[lv+1]//2)%(4*PatchSize) == 0 ) ) # Make sure the boundaries of range lv and lv+1 are separated by at least 1 patch of lv
    UM_IC_N_y[lv] = max( UM_IC_N_y[lv],  UM_IC_N_y[lv+1]//2 + 2*PatchSize + 2*PatchSize*int( (UM_IC_N_y[lv+1]//2)%(4*PatchSize) == 0 ) )
    UM_IC_N_z[lv] = max( UM_IC_N_z[lv],  UM_IC_N_z[lv+1]//2 + 2*PatchSize + 2*PatchSize*int( (UM_IC_N_z[lv+1]//2)%(4*PatchSize) == 0 ) )

UM_IC_Width_x       = UM_IC_N_x*UM_IC_dh   # width on each level in the x direction
UM_IC_Width_y       = UM_IC_N_y*UM_IC_dh   # ...                        y direction
UM_IC_Width_z       = UM_IC_N_z*UM_IC_dh   # ...                        z direction

UM_IC_x0            = np.zeros( UM_IC_NLEVEL )                   # left edge of the refinement region for each level in the x direction
UM_IC_y0            = np.zeros( UM_IC_NLEVEL )                   # left ...                                                 y direction
UM_IC_z0            = np.zeros( UM_IC_NLEVEL )                   # left ...                                                 z direction
UM_IC_x1            = np.zeros( UM_IC_NLEVEL )                   # right ...                                                x direction
UM_IC_y1            = np.zeros( UM_IC_NLEVEL )                   # right ...                                                y direction
UM_IC_z1            = np.zeros( UM_IC_NLEVEL )                   # right ...                                                z direction
UM_IC_dLv           = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # Lv
UM_IC_NP_Skip_xL    = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # number of patches on the parent level to be skipped in the x direction from the left edge of the parent refinement region
UM_IC_NP_Skip_xR    = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # ...                                                        x ...                right ...
UM_IC_NP_Skip_yL    = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # ...                                                        y ...                left  ...
UM_IC_NP_Skip_yR    = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # ...                                                        y ...                right ...
UM_IC_NP_Skip_zL    = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # ...                                                        z ...                left ...
UM_IC_NP_Skip_zR    = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # ...                                                        z ...                right ...
UM_IC_NP_x          = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # number of patches on each level in the x direction
UM_IC_NP_y          = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # ...                                    y direction
UM_IC_NP_z          = np.zeros( UM_IC_NLEVEL, dtype=np.uint32 )  # ...                                    z direction

for lv in range( 0, UM_IC_NLEVEL, 1 ):
    UM_IC_dLv       [lv] = lv
    UM_IC_x0        [lv] = 0.0                      if lv == 0 else 0.5*UM_IC_Width_x[0] - 0.5*UM_IC_N_x[lv]*UM_IC_dh[lv]
    UM_IC_y0        [lv] = 0.0                      if lv == 0 else 0.5*UM_IC_Width_y[0] - 0.5*UM_IC_N_y[lv]*UM_IC_dh[lv]
    UM_IC_z0        [lv] = 0.0                      if lv == 0 else 0.5*UM_IC_Width_z[0] - 0.5*UM_IC_N_z[lv]*UM_IC_dh[lv]
    UM_IC_x1        [lv] = UM_IC_Width_x[0]         if lv == 0 else 0.5*UM_IC_Width_x[0] + 0.5*UM_IC_N_x[lv]*UM_IC_dh[lv]
    UM_IC_y1        [lv] = UM_IC_Width_y[0]         if lv == 0 else 0.5*UM_IC_Width_y[0] + 0.5*UM_IC_N_y[lv]*UM_IC_dh[lv]
    UM_IC_z1        [lv] = UM_IC_Width_z[0]         if lv == 0 else 0.5*UM_IC_Width_z[0] + 0.5*UM_IC_N_z[lv]*UM_IC_dh[lv]
    UM_IC_NP_Skip_xL[lv] = 0                        if lv == 0 else int( ( UM_IC_x0[lv  ] - UM_IC_x0[lv-1] )/( PatchSize * UM_IC_dh[lv-1] ) )
    UM_IC_NP_Skip_xR[lv] = 0                        if lv == 0 else int( ( UM_IC_x1[lv-1] - UM_IC_x1[lv  ] )/( PatchSize * UM_IC_dh[lv-1] ) )
    UM_IC_NP_Skip_yL[lv] = 0                        if lv == 0 else int( ( UM_IC_y0[lv  ] - UM_IC_y0[lv-1] )/( PatchSize * UM_IC_dh[lv-1] ) )
    UM_IC_NP_Skip_yR[lv] = 0                        if lv == 0 else int( ( UM_IC_y1[lv-1] - UM_IC_y1[lv  ] )/( PatchSize * UM_IC_dh[lv-1] ) )
    UM_IC_NP_Skip_zL[lv] = 0                        if lv == 0 else int( ( UM_IC_z0[lv  ] - UM_IC_z0[lv-1] )/( PatchSize * UM_IC_dh[lv-1] ) )
    UM_IC_NP_Skip_zR[lv] = 0                        if lv == 0 else int( ( UM_IC_z1[lv-1] - UM_IC_z1[lv  ] )/( PatchSize * UM_IC_dh[lv-1] ) )
    UM_IC_NP_x      [lv] = UM_IC_N_x[0]/PatchSize   if lv == 0 else 2*(UM_IC_NP_x[lv-1]-UM_IC_NP_Skip_xL[lv]-UM_IC_NP_Skip_xR[lv])
    UM_IC_NP_y      [lv] = UM_IC_N_y[0]/PatchSize   if lv == 0 else 2*(UM_IC_NP_y[lv-1]-UM_IC_NP_Skip_yL[lv]-UM_IC_NP_Skip_yR[lv])
    UM_IC_NP_z      [lv] = UM_IC_N_z[0]/PatchSize   if lv == 0 else 2*(UM_IC_NP_z[lv-1]-UM_IC_NP_Skip_zL[lv]-UM_IC_NP_Skip_zR[lv])



# 1.0 Load Data
# Reference: https://yt-project.org/docs/dev/examining/loading_data.html#units-and-bounding-boxes
bbox = [ [-0.5*UM_IC_Width_x[0], 0.5*UM_IC_Width_x[0]],
         [-0.5*UM_IC_Width_y[0], 0.5*UM_IC_Width_y[0]],
         [-0.5*UM_IC_Width_z[0], 0.5*UM_IC_Width_z[0]] ]

unit_base = {
              'UnitLength_in_cm'        : GIZMO_UnitLength_in_cm,
              'UnitMass_in_g'           : GIZMO_UnitMass_in_g,
              'UnitVelocity_in_cm_per_s': GIZMO_UnitVelocity_in_cm_per_s,
            }

ds = yt.load( GIZMO_Filename, unit_base=unit_base, bounding_box=bbox )


# 1.1 Print GIZMO data information
if VERBOSE:
    print( '' )
    print( '------------------------------------------------------------------------------------------------' )
    print( 'GIZMO data information' )
    print( '------------------------------------------------------------------------------------------------' )
    print( f'{GIZMO_Filename                  = }' )
    print( f'{GIZMO_UnitLength_in_cm          = }' )
    print( f'{GIZMO_UnitMass_in_g             = }' )
    print( f'{GIZMO_UnitVelocity_in_cm_per_s  = }' )
    print( f'GIZMO_HubbleParam (h)           = {ds.parameters["HubbleParam"]}' )

    print( '' )
    print( 'parameters:' )
    for i in ds.parameters:
        print( '   ', i, '  =', ds.parameters[i] )

    print( f'   {ds.sph_smoothing_style   = }' )
    print( f'   {ds.use_sph_normalization = }' )
    print( f'   {ds.default_kernel_name   = }' )

    print( '' )
    print( 'particle_type_counts:' )
    for i in ds.particle_type_counts:
        print( '  ', i, '  =', ds.particle_type_counts[i] )

    print( '' )
    print( 'units:' )
    print( f'   Note h = {ds.parameters["HubbleParam"] = } = {ds.hubble_constant = }' )
    print( f'   {ds.quan( 1.0, "code_length"   ).in_units("cm/h")   = }' )
    print(  '                                                      = ', ds.quan( 1.0, "code_length"   ).in_units("cm")      )
    print(  '                                                      = ', ds.quan( 1.0, "code_length"   ).in_units("pc")      )
    print(  '                                                      = ', ds.quan( 1.0, "code_length"   ).in_units("kpc")     )
    print( '' )
    print( f'   {ds.quan( 1.0, "code_velocity" ).in_units("cm/s")   = }' )
    print(  '                                                      = ', ds.quan( 1.0, "code_velocity" ).in_units("km/s")    )
    print(  '                                                      = ', ds.quan( 1.0, "code_velocity" ).in_units("pc/Myr")  )
    print(  '                                                      = ', ds.quan( 1.0, "code_velocity" ).in_units("kpc/Myr") )
    print( '' )
    print( f'   {ds.quan( 1.0, "code_mass"     ).in_units("g/h")    = }' )
    print(  '                                                      = ', ds.quan( 1.0, "code_mass"     ).in_units("g")       )
    print(  '                                                      = ', ds.quan( 1.0, "code_mass"     ).in_units("Msun/h")  )
    print(  '                                                      = ', ds.quan( 1.0, "code_mass"     ).in_units("Msun")    )
    print( '' )
    print( f'   {ds.quan( 1.0, "code_time"     ).in_units("s/h")    = }' )
    print(  '                                                      = ', ds.quan( 1.0, "code_time"     ).in_units("s")       )
    print(  '                                                      = ', ds.quan( 1.0, "code_time"     ).in_units("Myr")     )
    print(  '                                                      = ', ds.quan( 1.0, "code_time"     ).in_units("Gyr")     )

    print( '' )
    print( 'field_list:' )
    for i in ds.field_list+ds.derived_field_list:
        print( '  ', i )

    print( '------------------------------------------------------------------------------------------------' )
    print( '\n' )



# 2.0 Particle information
ad = ds.all_data()
if VERBOSE:
    print( '' )
    print( '------------------------------------------------------------------------------------------------' )
    print( 'Particle (in the bounding_box) information' )
    print( '------------------------------------------------------------------------------------------------' )
    print( 'PartType1 particles:' )
    print( '   particle_mass       :', stats.describe( ad['PartType1', 'particle_mass'],       axis=None ) )
    print( '   particle_position_x :', stats.describe( ad['PartType1', 'particle_position_x'], axis=None ) )
    print( '   particle_position_y :', stats.describe( ad['PartType1', 'particle_position_y'], axis=None ) )
    print( '   particle_position_z :', stats.describe( ad['PartType1', 'particle_position_z'], axis=None ) )
    print( '   particle_velocity_x :', stats.describe( ad['PartType1', 'particle_velocity_x'], axis=None ) )
    print( '   particle_velocity_y :', stats.describe( ad['PartType1', 'particle_velocity_y'], axis=None ) )
    print( '   particle_velocity_z :', stats.describe( ad['PartType1', 'particle_velocity_z'], axis=None ) )
    print( '' )
    print( 'PartType2 particles:' )
    print( '   particle_mass       :', stats.describe( ad['PartType2', 'particle_mass'],       axis=None ) )
    print( '   particle_position_x :', stats.describe( ad['PartType2', 'particle_position_x'], axis=None ) )
    print( '   particle_position_y :', stats.describe( ad['PartType2', 'particle_position_y'], axis=None ) )
    print( '   particle_position_z :', stats.describe( ad['PartType2', 'particle_position_z'], axis=None ) )
    print( '   particle_velocity_x :', stats.describe( ad['PartType2', 'particle_velocity_x'], axis=None ) )
    print( '   particle_velocity_y :', stats.describe( ad['PartType2', 'particle_velocity_y'], axis=None ) )
    print( '   particle_velocity_z :', stats.describe( ad['PartType2', 'particle_velocity_z'], axis=None ) )
    print( '' )
    print( 'PartType0 particles:' )
    print( '   particle_mass       :', stats.describe( ad['PartType0', 'particle_mass'],       axis=None ) )
    print( '   particle_position_x :', stats.describe( ad['PartType0', 'particle_position_x'], axis=None ) )
    print( '   particle_position_y :', stats.describe( ad['PartType0', 'particle_position_y'], axis=None ) )
    print( '   particle_position_z :', stats.describe( ad['PartType0', 'particle_position_z'], axis=None ) )
    print( '   particle_velocity_x :', stats.describe( ad['PartType0', 'particle_velocity_x'], axis=None ) )
    print( '   particle_velocity_y :', stats.describe( ad['PartType0', 'particle_velocity_y'], axis=None ) )
    print( '   particle_velocity_z :', stats.describe( ad['PartType0', 'particle_velocity_z'], axis=None ) )
    print( '   smoothing_length    :', stats.describe( ad['PartType0', 'smoothing_length'],    axis=None ) )
    print( '   Density             :', stats.describe( ad['PartType0', 'Density'],             axis=None ) )
    print( '   InternalEnergy      :', stats.describe( ad['PartType0', 'InternalEnergy'],      axis=None ) )
    print( '' )
    print( 'Total particle mass:' )
    print( '   PartType1: {: >16.8e}'.format( ad.quantities.total_quantity( ('PartType1', 'particle_mass') ).in_units('Msun') ) )
    print( '   PartType2: {: >16.8e}'.format( ad.quantities.total_quantity( ('PartType2', 'particle_mass') ).in_units('Msun') ) )
    print( '   PartType0: {: >16.8e}'.format( ad.quantities.total_quantity( ('PartType0', 'particle_mass') ).in_units('Msun') ) )
    print( '   All : {: >16.8e}'.format( ad.quantities.total_quantity( ('all',  'particle_mass') ).in_units('Msun') ) )

    if MakeDiskGalaxy_R200 > 0:
        R200 = ds.quan( MakeDiskGalaxy_R200, 'kpc' )
        M200 = ds.sphere( ds.domain_center, R200 ).quantities.total_quantity( ('all', 'particle_mass') )
        print( '' )
        print( f'Total particle mass in {MakeDiskGalaxy_R200 = } kpc' )
        print( '   PartType1: {: >16.8e}'.format( ds.sphere( ds.domain_center, R200 ).quantities.total_quantity( ('PartType1', 'particle_mass') ).in_units('Msun') ) )
        print( '   PartType2: {: >16.8e}'.format( ds.sphere( ds.domain_center, R200 ).quantities.total_quantity( ('PartType2', 'particle_mass') ).in_units('Msun') ) )
        print( '   PartType0: {: >16.8e}'.format( ds.sphere( ds.domain_center, R200 ).quantities.total_quantity( ('PartType0', 'particle_mass') ).in_units('Msun') ) )
        print( '' )
        print( '   M200= {: >16.8e}'.format( M200.in_units('Msun') ) )
        print( '   V200= {: >16.8e}'.format( np.sqrt( ds.units.newtons_constant*M200/R200 ).in_units('km/s') ) )

    print( '------------------------------------------------------------------------------------------------' )
    print( '\n' )


# 2.1 Plot the projections for particles in GIZMO data
for lv in range( 0, UM_IC_NLEVEL, 1 ):

    Width_plot = ds.quan( np.max( [ UM_IC_Width_x[lv], UM_IC_Width_y[lv], UM_IC_Width_z[lv] ] ), 'code_length' ) # width to plot the projections

    for yz in [ 'particle_position_y', 'particle_position_z' ]:
        # PartType1
        p = yt.ParticlePlot( ds,
                            ('PartType1', 'particle_position_x'),
                            ('PartType1', yz ),
                            ('PartType1', 'particle_mass'), width=Width_plot )
        p.set_unit( ('PartType1', 'particle_mass'), 'code_mass' )
        p.set_cmap( ('PartType1', 'particle_mass'), 'inferno' )
        p.save( 'fig_Projection_Lv%02d_PartType1_particle_position_x-%s.png'%(lv,yz) )

        for dlv in range( lv, UM_IC_NLEVEL, 1 ):
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
        p.save( 'fig_Projection_Lv%02d_PartType1_particle_position_x-%s_withAMRboundary.png'%(lv,yz) )

        # PartType2
        p = yt.ParticlePlot( ds,
                            ('PartType2', 'particle_position_x'),
                            ('PartType2', yz ),
                            ('PartType2', 'particle_mass'), width=Width_plot )
        p.set_unit( ('PartType2', 'particle_mass'), 'code_mass' )
        p.set_cmap( ('PartType2', 'particle_mass'), 'arbre' )
        p.save( 'fig_Projection_Lv%02d_PartType2_particle_position_x-%s.png'%(lv,yz) )

        for dlv in range( lv, UM_IC_NLEVEL, 1 ):
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
        p.save( 'fig_Projection_Lv%02d_PartType2_particle_position_x-%s_withAMRboundary.png'%(lv,yz) )

        # PartType0
        p = yt.ParticlePlot( ds,
                            ('PartType0', 'particle_position_x'),
                            ('PartType0', yz ),
                            ('PartType0', 'particle_mass'), width=Width_plot )
        p.set_unit( ('PartType0', 'particle_mass'), 'code_mass' )
        p.set_cmap( ('PartType0', 'particle_mass'), 'viridis' )
        p.save( 'fig_Projection_Lv%02d_PartType0_particle_position_x-%s.png'%(lv,yz) )

        for dlv in range( lv, UM_IC_NLEVEL, 1 ):
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
        p.save( 'fig_Projection_Lv%02d_PartType0_particle_position_x-%s_withAMRboundary.png'%(lv,yz) )

    for yz in [ 'y', 'z' ]:
        # PartType0 density
        p = yt.ProjectionPlot( ds, yz, ('gas', 'density'), width=Width_plot )
        p.set_unit( ('gas', 'density'), 'code_density*code_length' )
        p.set_cmap( ('gas', 'density'), 'viridis' )
        p.set_zlim( ('gas', 'density'), 1.0e-8, 1.0e-2 )
        p.save( 'fig_Projection_Lv%02d_%s_gas_density.png'%(lv,yz) )

        for dlv in range( lv, UM_IC_NLEVEL, 1 ):
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y1[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x1[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
            p.annotate_line( [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z1[dlv]-0.5*UM_IC_Width_z[0]], [UM_IC_x0[dlv]-0.5*UM_IC_Width_x[0], UM_IC_y0[dlv]-0.5*UM_IC_Width_y[0], UM_IC_z0[dlv]-0.5*UM_IC_Width_z[0]], coord_system='data', color='red' )
        p.save( 'fig_Projection_Lv%02d_%s_gas_density_withAMRboundary.png'%(lv,yz) )

        # PartType0 momentum density x
        p = yt.ProjectionPlot( ds, yz, ('gas', 'momentum_density_x'), width=Width_plot )
        p.set_unit( ('gas', 'momentum_density_x'), 'code_velocity*code_density*code_length' )
        p.set_cmap( ('gas', 'momentum_density_x'), 'viridis' )
        p.save( 'fig_Projection_Lv%02d_%s_gas_momentum_density_x.png'%(lv,yz) )

        # PartType0 momentum density y
        p = yt.ProjectionPlot( ds, yz, ('gas', 'momentum_density_y'), width=Width_plot )
        p.set_unit( ('gas', 'momentum_density_y'), 'code_velocity*code_density*code_length' )
        p.set_cmap( ('gas', 'momentum_density_y'), 'viridis' )
        p.save( 'fig_Projection_Lv%02d_%s_gas_momentum_density_y.png'%(lv,yz) )

        # PartType0 momentum density z
        p = yt.ProjectionPlot( ds, yz, ('gas', 'momentum_density_z'), width=Width_plot )
        p.set_unit( ('gas', 'momentum_density_z'), 'code_velocity*code_density*code_length' )
        p.set_cmap( ('gas', 'momentum_density_z'), 'viridis' )
        p.save( 'fig_Projection_Lv%02d_%s_gas_momentum_density_z.png'%(lv,yz) )

        # PartType0 temperature
#        p = yt.ProjectionPlot( ds, yz, ('gas', 'kT'), width=Width_plot )
#        p.set_unit( ('gas', 'kT'), 'code_velocity**2*code_mass*code_length' )
#        p.set_cmap( ('gas', 'kT'), 'viridis' )
#        p.set_zlim( ('gas', 'kT'), 1.0e-69, 1.0e-64 )
#        p.save( 'fig_Projection_Lv%02d_%s_gas_kT.png'%(lv,yz) )

        # PartType0 kinetic energy density
        p = yt.ProjectionPlot( ds, yz, ('gas', 'kinetic_energy_density'), width=Width_plot )
        p.set_unit( ('gas', 'kinetic_energy_density'), 'code_velocity**2*code_density*code_length' )
        p.set_cmap( ('gas', 'kinetic_energy_density'), 'viridis' )
        p.set_zlim( ('gas', 'kinetic_energy_density'), 1.0e-6, 1.0e0 )
        p.save( 'fig_Projection_Lv%02d_%s_gas_kinetic_energy_density.png'%(lv,yz) )



# 3.0 Prepare the particle data for PartType1 and PartType2
Array_PartType1_ParMass =   ad['PartType1', 'particle_mass'].in_units('code_mass').d
Array_PartType1_ParPosX = ( ad['PartType1', 'particle_position_x'] - ds.domain_center[0] + 0.5*ds.domain_width[0] ).in_units('code_length').d  # Shift the coordinates to be in [0, BoxSize]
Array_PartType1_ParPosY = ( ad['PartType1', 'particle_position_y'] - ds.domain_center[1] + 0.5*ds.domain_width[1] ).in_units('code_length').d
Array_PartType1_ParPosZ = ( ad['PartType1', 'particle_position_z'] - ds.domain_center[2] + 0.5*ds.domain_width[2] ).in_units('code_length').d
Array_PartType1_ParVelX =   ad['PartType1', 'particle_velocity_x'].in_units('code_velocity').d
Array_PartType1_ParVelY =   ad['PartType1', 'particle_velocity_y'].in_units('code_velocity').d
Array_PartType1_ParVelZ =   ad['PartType1', 'particle_velocity_z'].in_units('code_velocity').d
Array_PartType1_ParType =   ad['PartType1', 'particle_ones'].d * 2   # ParType = 2 in GAMER

Array_PartType2_ParMass =   ad['PartType2', 'particle_mass'].in_units('code_mass').d
Array_PartType2_ParPosX = ( ad['PartType2', 'particle_position_x'] - ds.domain_center[0] + 0.5*ds.domain_width[0] ).in_units('code_length').d  # Shift the coordinates to be in [0, BoxSize]
Array_PartType2_ParPosY = ( ad['PartType2', 'particle_position_y'] - ds.domain_center[1] + 0.5*ds.domain_width[1] ).in_units('code_length').d
Array_PartType2_ParPosZ = ( ad['PartType2', 'particle_position_z'] - ds.domain_center[2] + 0.5*ds.domain_width[2] ).in_units('code_length').d
Array_PartType2_ParVelX =   ad['PartType2', 'particle_velocity_x'].in_units('code_velocity').d
Array_PartType2_ParVelY =   ad['PartType2', 'particle_velocity_y'].in_units('code_velocity').d
Array_PartType2_ParVelZ =   ad['PartType2', 'particle_velocity_z'].in_units('code_velocity').d
Array_PartType2_ParType =   ad['PartType2', 'particle_ones'].d * 3   # ParType = 3 in GAMER


# 3.1 Plot the PartType1 and PartType2 particles
fig = plt.figure( 1, (8,3) )
fig.suptitle( 'ParPos' )
ax1 = fig.add_subplot(121)
ax1.scatter( Array_PartType1_ParPosX, Array_PartType1_ParPosY, c='r', s=0.5, label='PartType1', alpha=0.3, edgecolors='none' )
ax1.scatter( Array_PartType2_ParPosX, Array_PartType2_ParPosY, c='b', s=0.5, label='PartType2', alpha=0.3, edgecolors='none' )
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

ax2 = fig.add_subplot(122)
ax2.scatter( Array_PartType1_ParPosZ, Array_PartType1_ParPosX, c='r', s=0.5, label='PartType1', alpha=0.3, edgecolors='none' )
ax2.scatter( Array_PartType2_ParPosZ, Array_PartType2_ParPosX, c='b', s=0.5, label='PartType2', alpha=0.3, edgecolors='none' )
ax2.set_xlabel('Z')
ax2.set_ylabel('X')
ax2.legend()
plt.tight_layout()
fig.savefig( 'fig_Array_ParPos.png', dpi=300 )
plt.close()


# 3.2 Construct PAR_IC
Array_ParMass = np.concatenate( (Array_PartType1_ParMass, Array_PartType2_ParMass), axis=None ).astype( dtype=dtype_float_PAR_IC )
Array_ParPosX = np.concatenate( (Array_PartType1_ParPosX, Array_PartType2_ParPosX), axis=None ).astype( dtype=dtype_float_PAR_IC )
Array_ParPosY = np.concatenate( (Array_PartType1_ParPosY, Array_PartType2_ParPosY), axis=None ).astype( dtype=dtype_float_PAR_IC )
Array_ParPosZ = np.concatenate( (Array_PartType1_ParPosZ, Array_PartType2_ParPosZ), axis=None ).astype( dtype=dtype_float_PAR_IC )
Array_ParVelX = np.concatenate( (Array_PartType1_ParVelX, Array_PartType2_ParVelX), axis=None ).astype( dtype=dtype_float_PAR_IC )
Array_ParVelY = np.concatenate( (Array_PartType1_ParVelY, Array_PartType2_ParVelY), axis=None ).astype( dtype=dtype_float_PAR_IC )
Array_ParVelZ = np.concatenate( (Array_PartType1_ParVelZ, Array_PartType2_ParVelZ), axis=None ).astype( dtype=dtype_float_PAR_IC )
Array_ParType = np.concatenate( (Array_PartType1_ParType, Array_PartType2_ParType), axis=None ).astype( dtype=dtype_int_PAR_IC   )

PAR_NPAR      = len(Array_ParMass)

with open( 'PAR_IC', 'wb' ) as f:
    f.write( Array_ParMass.tobytes() )
    f.write( Array_ParPosX.tobytes() )
    f.write( Array_ParPosY.tobytes() )
    f.write( Array_ParPosZ.tobytes() )
    f.write( Array_ParVelX.tobytes() )
    f.write( Array_ParVelY.tobytes() )
    f.write( Array_ParVelZ.tobytes() )
    if OutputParCreTime:
        f.write( (   -1.0*np.ones_like( Array_ParMass ) ).astype( dtype=dtype_float_PAR_IC ).tobytes() )
    if OutputParSNIITime:
        f.write( ( np.inf*np.ones_like( Array_ParMass ) ).astype( dtype=dtype_float_PAR_IC ).tobytes() )
    if OutputParMetalFrac:
        f.write( (   -1.0*np.ones_like( Array_ParMass ) ).astype( dtype=dtype_float_PAR_IC ).tobytes() )
    f.write( Array_ParType.tobytes() )
    f.close()


# 3.3 Output PAR_IC information
print( '' )
print( '------------------------------------------------------------------------------------------------' )
print( 'PAR_IC information' )
print( '------------------------------------------------------------------------------------------------' )
print( f'{PAR_IC_FLOAT8      = }' )
print( f'{PAR_IC_INT8        = }' )
print( f'{PAR_NPAR           = }' )
print( '' )
print( f'{OutputParCreTime   = }' )
print( f'{OutputParSNIITime  = }' )
print( f'{OutputParMetalFrac = }' )
print( '' )
print( 'Array_ParMass      :', stats.describe( Array_ParMass,      axis=None ) )
print( 'Array_ParPosX      :', stats.describe( Array_ParPosX,      axis=None ) )
print( 'Array_ParPosY      :', stats.describe( Array_ParPosY,      axis=None ) )
print( 'Array_ParPosZ      :', stats.describe( Array_ParPosZ,      axis=None ) )
print( 'Array_ParVelX      :', stats.describe( Array_ParVelX,      axis=None ) )
print( 'Array_ParVelY      :', stats.describe( Array_ParVelY,      axis=None ) )
print( 'Array_ParVelZ      :', stats.describe( Array_ParVelZ,      axis=None ) )
print( 'Array_ParType      :', stats.describe( Array_ParType,      axis=None ) )

if VERBOSE:
    print( '' )
    print( 'Array_PartType1_ParMass :', stats.describe( Array_PartType1_ParMass, axis=None ) )
    print( 'Array_PartType1_ParPosX :', stats.describe( Array_PartType1_ParPosX, axis=None ) )
    print( 'Array_PartType1_ParPosY :', stats.describe( Array_PartType1_ParPosY, axis=None ) )
    print( 'Array_PartType1_ParPosZ :', stats.describe( Array_PartType1_ParPosZ, axis=None ) )
    print( 'Array_PartType1_ParVelX :', stats.describe( Array_PartType1_ParVelX, axis=None ) )
    print( 'Array_PartType1_ParVelY :', stats.describe( Array_PartType1_ParVelY, axis=None ) )
    print( 'Array_PartType1_ParVelZ :', stats.describe( Array_PartType1_ParVelZ, axis=None ) )
    print( 'Array_PartType1_ParType :', stats.describe( Array_PartType1_ParType, axis=None ) )
    print( '' )
    print( 'Array_PartType2_ParMass :', stats.describe( Array_PartType2_ParMass, axis=None ) )
    print( 'Array_PartType2_ParPosX :', stats.describe( Array_PartType2_ParPosX, axis=None ) )
    print( 'Array_PartType2_ParPosY :', stats.describe( Array_PartType2_ParPosY, axis=None ) )
    print( 'Array_PartType2_ParPosZ :', stats.describe( Array_PartType2_ParPosZ, axis=None ) )
    print( 'Array_PartType2_ParVelX :', stats.describe( Array_PartType2_ParVelX, axis=None ) )
    print( 'Array_PartType2_ParVelY :', stats.describe( Array_PartType2_ParVelY, axis=None ) )
    print( 'Array_PartType2_ParVelZ :', stats.describe( Array_PartType2_ParVelZ, axis=None ) )
    print( 'Array_PartType2_ParType :', stats.describe( Array_PartType2_ParType, axis=None ) )

print( '------------------------------------------------------------------------------------------------' )
print( '\n' )


del(Array_ParMass)
del(Array_ParPosX)
del(Array_ParPosY)
del(Array_ParPosZ)
del(Array_ParVelX)
del(Array_ParVelY)
del(Array_ParVelZ)
del(Array_ParType)
del(Array_PartType1_ParMass)
del(Array_PartType1_ParPosX)
del(Array_PartType1_ParPosY)
del(Array_PartType1_ParPosZ)
del(Array_PartType1_ParVelX)
del(Array_PartType1_ParVelY)
del(Array_PartType1_ParVelZ)
del(Array_PartType1_ParType)
del(Array_PartType2_ParMass)
del(Array_PartType2_ParPosX)
del(Array_PartType2_ParPosY)
del(Array_PartType2_ParPosZ)
del(Array_PartType2_ParVelX)
del(Array_PartType2_ParVelY)
del(Array_PartType2_ParVelZ)
del(Array_PartType2_ParType)
gc.collect()



# 4.0 Prepare AMR grid data for the PartType0 particles
if VERBOSE:
    print( '' )
    print( '------------------------------------------------------------------------------------------------' )
    print( 'Gas grid data information' )
    print( '------------------------------------------------------------------------------------------------' )

for lv in range( 0, UM_IC_NLEVEL, 1 ):
    print( '---Lv = %02d-----------------'%lv )
    print( '' )


#   4.1 Create Arbitrary Grid
#   Reference: https://yt-project.org/doc/analyzing/objects.html#arbitrary-grids-objects
    ag = ds.arbitrary_grid( ds.domain_center-0.5*ds.arr( [ UM_IC_Width_x[lv], UM_IC_Width_y[lv], UM_IC_Width_z[lv] ], 'code_length' ),
                            ds.domain_center+0.5*ds.arr( [ UM_IC_Width_x[lv], UM_IC_Width_y[lv], UM_IC_Width_z[lv] ], 'code_length' ),
                            dims=[ UM_IC_N_x[lv]*N_sub, UM_IC_N_y[lv]*N_sub, UM_IC_N_z[lv]*N_sub ] )

    ag_density                 = ag['PartType0', 'density']
    ag_velocity_x              = ag['PartType0', 'velocity_x']
    ag_velocity_y              = ag['PartType0', 'velocity_y']
    ag_velocity_z              = ag['PartType0', 'velocity_z']
    ag_specific_thermal_energy = ag['PartType0', 'specific_thermal_energy']

    if VERBOSE:
        print( 'Number of sub-sampling cells = ', N_sub )
        print( 'Arbitrary Grid, Lv%02d'%lv )
        print( '   density                 :', stats.describe( ag_density                , axis=None ) )
        print( '   velocity_x              :', stats.describe( ag_velocity_x             , axis=None ) )
        print( '   velocity_y              :', stats.describe( ag_velocity_y             , axis=None ) )
        print( '   velocity_z              :', stats.describe( ag_velocity_z             , axis=None ) )
        print( '   specific_thermal_energy :', stats.describe( ag_specific_thermal_energy, axis=None ) )
        # print( f'{ ag['gas', 'momentum_density_x']      = }' ) # didn't work


#   4.2 Save the arbitrary grid object as dataset
    #fn      = ag.save_as_dataset( 'Arbitrary_Grid_Lv%02d.h5'%lv,
    #                               fields=[ ('PartType0', 'density'),
    #                                        ('PartType0', 'velocity_x'),
    #                                        ('PartType0', 'velocity_y'),
    #                                        ('PartType0', 'velocity_z'),
    #                                        ('PartType0', 'specific_thermal_energy'),
    #                                      ] )
    #ds_grid = yt.load( fn )


#   4.3 Convert to conserved variables in GAMER
    Mtot_gas_ag  = ag_density.in_units('code_density').sum() * ds.quan( (UM_IC_dh[lv]/N_sub)**3, 'code_length**3' )
    Mtot_gas_par = ds.box(
                           ds.domain_center-0.5*ds.arr( [ UM_IC_Width_x[lv], UM_IC_Width_y[lv], UM_IC_Width_z[lv] ], 'code_length' ),
                           ds.domain_center+0.5*ds.arr( [ UM_IC_Width_x[lv], UM_IC_Width_y[lv], UM_IC_Width_z[lv] ], 'code_length' )
                         ).quantities.total_quantity( ('PartType0',  'particle_mass') ).in_units('code_mass')

    if VERBOSE:
        print( '' )
        print( 'Normalization for mass conservation' )
        print( '   Total mass of the enclosed SPH gas particles  = {: >16.8e}'.format( Mtot_gas_par                 ) )
        print( '   Total mass of the interpolated Arbitrary Grid = {: >16.8e}'.format( Mtot_gas_ag                  ) )
        print( '   particle to grid total mass ratio             = {: >16.8e}'.format( (Mtot_gas_par/Mtot_gas_ag).d ) )
        print( '   Normalization factor for interpolated density = {: >16.8e}'.format( factorMassNormalization      ) )

    subArray_DENS = ( factorMassNormalization*(     ag_density.in_units('code_density').d ) )
    subArray_MOMX = ( factorMassNormalization*(     ag_density.in_units('code_density').d * ag_velocity_x.in_units('code_velocity').d ) )
    subArray_MOMY = ( factorMassNormalization*(     ag_density.in_units('code_density').d * ag_velocity_y.in_units('code_velocity').d ) )
    subArray_MOMZ = ( factorMassNormalization*(     ag_density.in_units('code_density').d * ag_velocity_z.in_units('code_velocity').d ) )
    subArray_EINT = ( factorMassNormalization*(     ag_density.in_units('code_density').d * ag_specific_thermal_energy.in_units('code_velocity**2').d ) )
    subArray_EKIN = ( factorMassNormalization*( 0.5*ag_density.in_units('code_density').d *(ag_velocity_x**2 +
                                                                                            ag_velocity_y**2 +
                                                                                            ag_velocity_z**2 ).in_units('code_velocity**2').d ) )
    subArray_ENGY = ( subArray_EINT + subArray_EKIN )

    Array_DENS = ( ( 1.0/N_sub**3 )*subArray_DENS.reshape( UM_IC_N_x[lv], N_sub, UM_IC_N_y[lv], N_sub, UM_IC_N_z[lv], N_sub ).sum(axis=1).sum(axis=2).sum(axis=3) ).astype( dtype=dtype_UM_IC )
    Array_MOMX = ( ( 1.0/N_sub**3 )*subArray_MOMX.reshape( UM_IC_N_x[lv], N_sub, UM_IC_N_y[lv], N_sub, UM_IC_N_z[lv], N_sub ).sum(axis=1).sum(axis=2).sum(axis=3) ).astype( dtype=dtype_UM_IC )
    Array_MOMY = ( ( 1.0/N_sub**3 )*subArray_MOMY.reshape( UM_IC_N_x[lv], N_sub, UM_IC_N_y[lv], N_sub, UM_IC_N_z[lv], N_sub ).sum(axis=1).sum(axis=2).sum(axis=3) ).astype( dtype=dtype_UM_IC )
    Array_MOMZ = ( ( 1.0/N_sub**3 )*subArray_MOMZ.reshape( UM_IC_N_x[lv], N_sub, UM_IC_N_y[lv], N_sub, UM_IC_N_z[lv], N_sub ).sum(axis=1).sum(axis=2).sum(axis=3) ).astype( dtype=dtype_UM_IC )
    Array_EINT = ( ( 1.0/N_sub**3 )*subArray_EINT.reshape( UM_IC_N_x[lv], N_sub, UM_IC_N_y[lv], N_sub, UM_IC_N_z[lv], N_sub ).sum(axis=1).sum(axis=2).sum(axis=3) ).astype( dtype=dtype_UM_IC )
    Array_EKIN = ( ( 1.0/N_sub**3 )*subArray_EKIN.reshape( UM_IC_N_x[lv], N_sub, UM_IC_N_y[lv], N_sub, UM_IC_N_z[lv], N_sub ).sum(axis=1).sum(axis=2).sum(axis=3) ).astype( dtype=dtype_UM_IC )
    Array_ENGY = ( ( 1.0/N_sub**3 )*subArray_ENGY.reshape( UM_IC_N_x[lv], N_sub, UM_IC_N_y[lv], N_sub, UM_IC_N_z[lv], N_sub ).sum(axis=1).sum(axis=2).sum(axis=3) ).astype( dtype=dtype_UM_IC )

    if VERBOSE:
        print( '' )
        print( 'Conserved Variables, before adding background, Lv%02d'%lv )
        print( '   Array_DENS              : ', stats.describe( Array_DENS, axis=None ) )
        print( '   Array_EINT              : ', stats.describe( Array_EINT, axis=None ) )
        print( '   Array_EKIN              : ', stats.describe( Array_EKIN, axis=None ) )
        print( '   Array_ENGY              : ', stats.describe( Array_ENGY, axis=None ) )


#   4.4 Add the background to the density and energy density
    Dens_bg = ds.quan( Rho_bg, 'g/cm**3' )
    Eint_bg = max(ad['PartType0', 'InternalEnergy'])* Dens_bg

    Array_DENS[ Array_DENS <= 0.0 ] = Dens_bg.in_units('code_density').d.astype( dtype=dtype_UM_IC )
    Array_EINT[ Array_EINT <= 0.0 ] = Eint_bg.in_units('code_density*code_velocity**2').d.astype( dtype=dtype_UM_IC )
    Array_ENGY = Array_EINT + Array_EKIN

    if VERBOSE:
        print( '' )
        print( 'Conserved Variables, after adding background, Lv%02d'%lv )
        print( '   Array_DENS              : ', stats.describe( Array_DENS, axis=None ) )
        print( '   Array_MOMX              : ', stats.describe( Array_MOMX, axis=None ) )
        print( '   Array_MOMY              : ', stats.describe( Array_MOMY, axis=None ) )
        print( '   Array_MOMZ              : ', stats.describe( Array_MOMZ, axis=None ) )
        print( '   Array_ENGY              : ', stats.describe( Array_ENGY, axis=None ) )
        print( '' )
        print( f'{OutputMetal       = }' )
        if OutputMetal:
            print( f'{MetalMassFrac     = }' )
        print( f'{OutputDust        = }' )
        if OutputDust:
            print( f'{DustMassFrac      = }' )
        print( '\n' )


#   4.5 Plot the results of grid data
    for name, Array in zip( (     'DENS',     'MOMX',     'MOMY',     'MOMZ',     'EINT',     'EKIN',     'ENGY'),
                            (Array_DENS, Array_MOMX, Array_MOMY, Array_MOMZ, Array_EINT, Array_EKIN, Array_ENGY ) ):
        fig = plt.figure( 1, (8,3) )
        fig.suptitle( '%s, Lv%02d'%(name,lv) )

        scale = 'symlog' if name[0] == 'M' else 'log'

        ax1 = fig.add_subplot(121)
        pos = ax1.imshow( np.sum( Array, axis=2 ).T, norm=scale, origin='lower', interpolation='none' )
        fig.colorbar( pos, ax=ax1 )

        ax2 = fig.add_subplot(122)
        pos = ax2.imshow( np.sum( Array, axis=1 ), norm=scale, origin='lower', interpolation='none' )
        fig.colorbar( pos, ax=ax2 )

        plt.tight_layout()
        fig.savefig( 'fig_Array_%s_Lv%02d.png'%(name,lv), dpi=300 )
        plt.close()


#   4.6 Construct UM_IC
    mode = 'wb' if lv == 0 else 'ab'
    with open( 'UM_IC', mode ) as f:
        f.write( np.ascontiguousarray( np.swapaxes( Array_DENS, 0, 2 ) ).astype( dtype=dtype_UM_IC ).tobytes() )
        f.write( np.ascontiguousarray( np.swapaxes( Array_MOMX, 0, 2 ) ).astype( dtype=dtype_UM_IC ).tobytes() )
        f.write( np.ascontiguousarray( np.swapaxes( Array_MOMY, 0, 2 ) ).astype( dtype=dtype_UM_IC ).tobytes() )
        f.write( np.ascontiguousarray( np.swapaxes( Array_MOMZ, 0, 2 ) ).astype( dtype=dtype_UM_IC ).tobytes() )
        f.write( np.ascontiguousarray( np.swapaxes( Array_ENGY, 0, 2 ) ).astype( dtype=dtype_UM_IC ).tobytes() )
        if OutputMetal:
            f.write( ( np.ascontiguousarray( np.swapaxes( Array_DENS, 0, 2 ) )*MetalMassFrac ).astype( dtype=dtype_UM_IC ).tobytes() )
        if OutputDust:
            f.write( ( np.ascontiguousarray( np.swapaxes( Array_DENS, 0, 2 ) )*DustMassFrac ).astype( dtype=dtype_UM_IC ).tobytes() )
        f.close()

#   free memory
    del(subArray_DENS)
    del(subArray_MOMX)
    del(subArray_MOMY)
    del(subArray_MOMZ)
    del(subArray_EINT)
    del(subArray_EKIN)
    del(subArray_ENGY)
    del(Array_DENS)
    del(Array_MOMX)
    del(Array_MOMY)
    del(Array_MOMZ)
    del(Array_EINT)
    del(Array_EKIN)
    del(Array_ENGY)
    del(ag_density)
    del(ag_velocity_x)
    del(ag_velocity_y)
    del(ag_velocity_z)
    del(ag_specific_thermal_energy)
    ag.clear_data()
    gc.collect()

print( '------------------------------------------------------------------------------------------------' )
print( '\n' )


# 4.7 Create Input__UM_IC_RefineRegion
Input__UM_IC_RefineRegion_header =  ' For loading an AMR initial condition from a file --> OPT__INIT=3 && OPT__UM_IC_NLEVEL>1\n'
Input__UM_IC_RefineRegion_header += '\n'
Input__UM_IC_RefineRegion_header += ' dLv       : target AMR level = OPT__UM_IC_LEVEL + dLv\n'
Input__UM_IC_RefineRegion_header += ' NP_Skip_xL: number of patches on the parent level to be skipped in the x direction\n'
Input__UM_IC_RefineRegion_header += '             from the left edge of the parent refinement region\n'
Input__UM_IC_RefineRegion_header += ' ==================================================================================\n'
Input__UM_IC_RefineRegion_header += '       dLv  NP_Skip_xL  NP_Skip_xR  NP_Skip_yL  NP_Skip_yR  NP_Skip_zL  NP_Skip_zR'

np.savetxt( 'Input__UM_IC_RefineRegion',
            np.column_stack( (UM_IC_dLv[1:],
                              UM_IC_NP_Skip_xL[1:],
                              UM_IC_NP_Skip_xR[1:],
                              UM_IC_NP_Skip_yL[1:],
                              UM_IC_NP_Skip_yR[1:],
                              UM_IC_NP_Skip_zL[1:],
                              UM_IC_NP_Skip_zR[1:]) ),
            fmt='%11d', header=Input__UM_IC_RefineRegion_header, comments='#' )


# 4.8 Output UM_IC information
np.set_printoptions( formatter={'float':'{: 16.8e}'.format, 'int':'{: 16d}'.format}, linewidth=1000 )
print( '' )
print( '------------------------------------------------------------------------------------------------' )
print( 'UM_IC information' )
print( '------------------------------------------------------------------------------------------------' )
print( f'{PatchSize         = }' )
print( f'{UM_IC_FLOAT8      = }' )
print( '' )
print( f'{UM_IC_NLEVEL      = }' )
print( '' )
print( f'{UM_IC_dLv         = }' )
print( '' )
print( f'{UM_IC_N_x         = }' )
print( f'{UM_IC_N_y         = }' )
print( f'{UM_IC_N_z         = }' )
print( '' )
print( f'{UM_IC_Width_x     = }' )
print( f'{UM_IC_Width_y     = }' )
print( f'{UM_IC_Width_z     = }' )
print( '' )
print( f'{UM_IC_dh          = }' )
print( '' )
print( f'{UM_IC_NP_Skip_xL  = }' )
print( f'{UM_IC_NP_Skip_xR  = }' )
print( f'{UM_IC_x0          = }' )
print( f'{UM_IC_x1          = }' )
print( '' )
print( f'{UM_IC_NP_Skip_yL  = }' )
print( f'{UM_IC_NP_Skip_yR  = }' )
print( f'{UM_IC_y0          = }' )
print( f'{UM_IC_y1          = }' )
print( '' )
print( f'{UM_IC_NP_Skip_zL  = }' )
print( f'{UM_IC_NP_Skip_zR  = }' )
print( f'{UM_IC_z0          = }' )
print( f'{UM_IC_z1          = }' )
print( '' )
print( f'{UM_IC_NP_x        = }' )
print( f'{UM_IC_NP_y        = }' )
print( f'{UM_IC_NP_z        = }' )
print( '------------------------------------------------------------------------------------------------' )
print( '' )



# 5.0 Output GAMER Input__Parameter information
print( '' )
print( '------------------------------------------------------------------------------------------------' )
print( 'GAMER Input__Parameter information' )
print( '------------------------------------------------------------------------------------------------' )
print( 'BOX_SIZE                     % 21.14e'%(  np.max( [ UM_IC_Width_x[0], UM_IC_Width_y[0], UM_IC_Width_z[0] ] ) ) )
print( 'NX0_TOT_X                    % 2d'%(     UM_IC_N_x[0]                                                        ) )
print( 'NX0_TOT_Y                    % 2d'%(     UM_IC_N_y[0]                                                        ) )
print( 'NX0_TOT_Z                    % 2d'%(     UM_IC_N_z[0]                                                        ) )
print( 'OPT__UNIT                    % 2d'%(     1                                                                   ) )
print( 'UNIT_L                       % 21.14e'%( ds.quan( 1.0, 'code_length'   ).in_units('cm').d                    ) )
print( 'UNIT_M                       % 21.14e'%( ds.quan( 1.0, 'code_mass'     ).in_units('g').d                     ) )
print( 'UNIT_V                       % 21.14e'%( ds.quan( 1.0, 'code_velocity' ).in_units('cm/s').d                  ) )
print( 'PAR_NPAR                     % 2d'%(     PAR_NPAR                                                            ) )
print( 'PAR_INIT                     % 2d'%(     3                                                                   ) )
print( 'PAR_IC_FORMAT                % 2d'%(     1                                                                   ) )
print( 'PAR_IC_FLOAT8                % 2d'%(     PAR_IC_FLOAT8                                                       ) )
print( 'PAR_IC_INT8                  % 2d'%(     PAR_IC_INT8                                                         ) )
print( 'PAR_IC_MASS                  % 2d'%(    -1.0                                                                 ) )
print( 'PAR_IC_TYPE                  % 2d'%(    -1                                                                   ) )
print( 'MAX_LEVEL                    % 2d'%(     UM_IC_NLEVEL-1                                                      ) )
print( 'OPT__INIT                    % 2d'%(     3                                                                   ) )
print( 'OPT__UM_IC_LEVEL             % 2d'%(     0                                                                   ) )
print( 'OPT__UM_IC_NLEVEL            % 2d'%(     UM_IC_NLEVEL                                                        ) )
print( 'OPT__UM_IC_NVAR              % 2d'%(    -1                                                                   ) )
print( 'OPT__UM_IC_FORMAT            % 2d'%(     1                                                                   ) )
print( 'OPT__UM_IC_FLOAT8            % 2d'%(     UM_IC_FLOAT8                                                        ) )
print( '------------------------------------------------------------------------------------------------' )
print( '' )
