import argparse
from matplotlib import patheffects
import sys
import yt
import matplotlib.pyplot as plt
import WLMDwarfGalaxy_load_datasets
import WLMDwarfGalaxy_derived_fields

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-c', action='store', required=False, type=str, dest='code',
                     help='simulation code [%(default)s]', default='GAMER' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start    = args.idx_start
idx_end      = args.idx_end
didx         = args.didx
prefix       = args.prefix
code         = args.code

colormap     = {
                 'density'                   :'viridis',
                 'dust'                      :'viridis',
                 'dust2gas'                  :'viridis',
                 'T'                         :'magma',
                 'kinetic_energy_density'    :'plasma',
                 'magnetic_energy_density'   :'plasma',
                 'cosmic_ray_energy_density' :'plasma',
                 'plasma_beta'               :'plasma',
                 'particle_density_on_grid'  :'algae',
                 'velocity_magnitude'        :'RdPu',
                 'particle'                  :'algae',
                 'resolution_size'           :'cividis_r',
               }
field_unit   = {
                 'density'                   :'Msun/pc**3',
                 'dust'                      :'Msun/pc**3',
                 'dust2gas'                  :'dimensionless',
                 'T'                         :'K',
                 'kinetic_energy_density'    :'Msun/pc**3*km**2/s**2',
                 'magnetic_energy_density'   :'G**2',
                 'cosmic_ray_energy_density' :'Msun/pc**3*km**2/s**2',
                 'plasma_beta'               :'dimensionless',
                 'particle_density_on_grid'  :'Msun/pc**3',
                 'velocity_magnitude'        :'km/s',
                 'resolution_size'           :'pc',
               }
zlim         = {
                 'density_s'                  :(1.0e-7, 1.0e+0),
                 'density_p'                  :(3.0e-2, 3.0e+1),
                 'dust_s'                     :(1.0e-10,1.0e-3),
                 'dust_p'                     :(3.0e-10,1.0e-3),
                 'dust2gas_s'                 :(5.0e-4, 2.0e-3),
                 'dust2gas_p'                 :(5.0e-4, 2.0e-3),
                 'T_s'                        :(1.0e+1, 1.0e+8),
                 'T_p'                        :(3.0e+0, 3.0e+6),
                 'kinetic_energy_density_s'   :(1.0e-3, 1.0e+2),
                 'kinetic_energy_density_p'   :(1.0e+2, 1.0e+6),
                 'magnetic_energy_density_s'  :(1.0e-28,1.0e-18),
                 'magnetic_energy_density_p'  :(1.0e-24,1.0e-14),
                 'cosmic_ray_energy_density_s':(1.0e-4, 1.0e+1),
                 'cosmic_ray_energy_density_p':(1.0e+1, 1.0e+5),
                 'plasma_beta_s'              :(1.0e+6, 1.0e+16),
                 'plasma_beta_p'              :(1.0e+6, 1.0e+16),
                 'particle_density_on_grid_s' :(1.0e-5, 1.0e+0),
                 'particle_density_on_grid_p' :(1.0e+0, 1.0e+4),
                 'velocity_magnitude_s'       :(1.0e+1, 1.5e+2),
                 'resolution_size_s'          :(1.0e+0, 1.0e+2),
                 'resolution_size_p'          :(1.0e+0, 1.0e+2),
                 ('all',      'particle_mass'):(1.0e+1, 1.0e+5),
                 ('Halo',     'particle_mass'):(1.0e+4, 1.0e+6),
                 ('Disk',     'particle_mass'):(3.0e+1, 1.0e+3),
                 ('new_star', 'particle_mass'):(3.0e-1, 3.0e+1),
                 ('exp_SNII', 'particle_mass'):(1.0e+1, 1.0e+3),
               }
zoomed_width = {
                 'a': 450.0,
                 'b':  30.0,
                 'c':  16.0,
                 'm':   6.0,
                 'n':   3.0,
               }
dpi          = 150

# options
hasDust   = True
plotPar   = False
printGrid = False

yt.enable_parallelism()

# load the dataset
ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)

WLMDwarfGalaxy_derived_fields.set_particle_types(code)

# main loop
for ds in ts.piter():

   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds, hasDust)

#  decide output fields
   fields_list  = []
   fields_list.append( 'density'                  )
#   fields_list.append( 'T'                        )
#   fields_list.append( 'kinetic_energy_density'   )
#   fields_list.append( 'velocity_magnitude'       )
#   fields_list.append( 'resolution_size'          )
#   fields_list.append( 'particle_density_on_grid' ) if code == 'GAMER' else None
#   if  ('gas', 'dust'                     ) in ds.derived_field_list: fields_list.append('dust2gas'                  )
#   if  ('gas', 'magnetic_energy_density'  ) in ds.derived_field_list: fields_list.append('magnetic_energy_density'   )
#   if  ('gas', 'plasma_beta'              ) in ds.derived_field_list: fields_list.append('plasma_beta'               )
#   if  ('gas', 'cosmic_ray_energy_density') in ds.derived_field_list: fields_list.append('cosmic_ray_energy_density' )

   pfields_list = []
#   pfields_list.append( ('all',      'particle_mass') )
#   pfields_list.append( ('Halo',     'particle_mass') )
#   pfields_list.append( ('Disk',     'particle_mass') )
#   pfields_list.append( ('new_star', 'particle_mass') ) if ('new_star', 'particle_mass') in ds.derived_field_list else None
#   pfields_list.append( ('exp_SNII', 'particle_mass') ) if code == 'GAMER'                                        else None

   # decide the center
   center = ds.domain_center

   for zoom_mode in ['m']:

      # zoom in
      width_kpc = zoomed_width[zoom_mode]

      for direction in ['x', 'z']:
         for field in fields_list:

            # slices
            s = yt.SlicePlot( ds, direction, field, center=center, width=(width_kpc, 'kpc'), buff_size=(1024, 1024) )
            s.set_axes_unit( 'kpc' )
            s.set_unit( field, field_unit[field] )
            s.set_cmap( field, colormap[field] )
            s.set_zlim( field, zlim[field+'_s'][0], zlim[field+'_s'][1] )
            if field == 'dust2gas':
               s.set_colorbar_label(field, r"$\rho_{\rm d}/\rho_{\rm g}$")
            if field == 'velocity_magnitude':
               s.set_log ( field, False )
               if direction == 'z':
                  s.annotate_quiver('velocity_x', 'velocity_y', factor=16)
            if field == 'magnetic_energy_density' and ('gas', 'magnetic_energy_density') in ds.derived_field_list:
               #s.annotate_magnetic_field(normalize=True)
               s.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), color='black') # linewidth=("gas", "magnetic_field_strength")
               #s.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"))

            s.annotate_timestamp( time_unit='Myr', corner='upper_right' )
            s.annotate_text( (0.02, 0.88), '%s'%(field), coord_system='axis', text_args={'color':'w', 'path_effects':[patheffects.withStroke(linewidth=2, foreground='k')]} )
            s.save( './imgs_a/fig_%s_%s_Slice_%s_%s.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
            if code == 'GAMER' and printGrid:
               s.annotate_grids( periodic=False )
               s.save( './imgs_a/fig_%s_%s_Slice_%s_%s_withgrids.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )


            if field == 'velocity_magnitude':
               continue

            weight_field =  'density_square'          if field == 'T'           else \
                            'density'                 if field == 'dust2gas'    else \
                            'magnetic_energy_density' if field == 'plasma_beta' else None
            project_unit = '' if field == 'T' or field == 'resolution_size' or field == 'dust2gas' or field == 'plasma_beta' else '*pc'
            proj_method  = 'min' if field == 'resolution_size' else 'integrate'


            # projections
            if field == 'T' or field == 'dust2gas' or field == 'magnetic_energy_density' or field == 'plasma_beta':
               width_x = 1 if direction == 'x' else width_kpc
               width_z = 1 if direction == 'z' else width_kpc
               box  = ds.box( center - 0.5*ds.arr([width_x, width_kpc, width_z], 'kpc'),
                              center + 0.5*ds.arr([width_x, width_kpc, width_z], 'kpc') )
               p = yt.ProjectionPlot( ds, direction, field, data_source=box, center=center, width=(width_kpc, 'kpc'), method=proj_method, weight_field=weight_field, buff_size=(1024, 1024) )
            else:
               p = yt.ProjectionPlot( ds, direction, field, center=center, width=(width_kpc, 'kpc'), method=proj_method, weight_field=weight_field, buff_size=(1024, 1024) )
            p.set_background_color( field )
            p.set_axes_unit( 'kpc' )
            p.set_unit( field, field_unit[field]+project_unit )
            p.set_cmap( field, colormap[field] )
            p.set_zlim( field, zlim[field+'_p'][0], zlim[field+'_p'][1] )
            if field == 'dust2gas':
               p.set_colorbar_label(field, r"$\Sigma_{\rm d}/\Sigma_{\rm g}$")
            p.annotate_timestamp( time_unit='Myr', corner='upper_right' )
            p.annotate_text( (0.02, 0.88), '%s'%(field), coord_system='axis', text_args={'color':'w', 'path_effects':[patheffects.withStroke(linewidth=2, foreground='k')]} )
            if field == 'magnetic_energy_density' and ('gas', 'magnetic_energy_density') in ds.derived_field_list:
               #p.annotate_magnetic_field(normalize=True)
               p.annotate_streamlines(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), color='black')
               #p.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), lim=(0.5, 0.65))
            try:
               p.save( './imgs_a/fig_%s_%s_Projection_%s_%s.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
               if code == 'GAMER' and printGrid:
                  p.annotate_grids( periodic=False )
                  p.save( './imgs_a/fig_%s_%s_Projection_%s_%s_withgrids.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                  p.clear_annotations( index=-1 )
               if field == 'density' and plotPar:
                  if ('new_star', 'particle_mass') in ds.derived_field_list:
                     p.annotate_particles( (width_kpc, 'kpc'), ptype='new_star', p_size=1, col='w', alpha=1.0, marker='o' )
                     p.save( './imgs_a/fig_%s_%s_Projection_%s_%s_withStars.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                  if code == 'GAMER':
                     p.annotate_particles( (width_kpc, 'kpc'), ptype='exp_SNII', p_size=1, col='r', alpha=1.0, marker='o' )
                     p.save( './imgs_a/fig_%s_%s_Projection_%s_%s_withSNe.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                     p.clear_annotations( index=-1 )
                  if ('young_star', 'particle_mass') in ds.derived_field_list:
                     p.clear_annotations( index=-1 )
                     p.annotate_particles( (width_kpc, 'kpc'), ptype='young_star', p_size=10, col='w', alpha=1.0, marker='o' )
                     p.save( './imgs_a/fig_%s_%s_Projection_%s_%s_withYStars.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                  if code == 'GAMER':
                     p.annotate_particles( (width_kpc, 'kpc'), ptype='young_SNII', p_size=10, col='r', alpha=1.0, marker='o' )
                     p.save( './imgs_a/fig_%s_%s_Projection_%s_%s_withYSNe.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
            except Exception as e:
               print( e )
               pass


         # particle plots
         for pfield in pfields_list:

            p = yt.ParticleProjectionPlot( ds, direction, pfield, center=center, width=(width_kpc, 'kpc'), depth=(width_kpc, 'kpc') )
            p.set_unit( pfield, 'Msun' )
            p.set_zlim( pfield, zlim[pfield][0], zlim[pfield][1] )
            p.set_cmap( pfield, colormap['particle'] )
            p.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
            p.save( './imgs_a/fig_%s_%s_Particles_%s_%s.png'%(ds, zoom_mode, direction, pfield[0]), mpl_kwargs={'dpi':dpi} )

            p = yt.ParticlePlot( ds, (pfield[0], 'particle_velocity_'+direction ), (pfield[0], 'particle_velocity_y'), (pfield[0], 'particle_mass') )
            p.set_unit( pfield, 'Msun' )
            p.set_unit( (pfield[0], 'particle_velocity_'+direction), 'km/s' )
            p.set_unit( (pfield[0], 'particle_velocity_y'),          'km/s' )
            p.set_xlim( -100.0, 100.0 )
            p.set_ylim( -100.0, 100.0 )
            p.set_zlim( pfield, zlim[pfield][0], zlim[pfield][1] )
            p.set_cmap( pfield, colormap['particle'] )
            p.annotate_text( xpos=-96, ypos=76, text='%s\n%s'%(pfield), color='w', path_effects=[patheffects.withStroke(linewidth=2, foreground='k')] )
            p.save( './imgs_a/fig_%s_%s_Particles_v%s_vy_%s_%s.png'%(ds, zoom_mode, direction, pfield[0], pfield[1]), mpl_kwargs={'dpi':dpi} )
