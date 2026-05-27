import yt
import numpy as np
from matplotlib import patheffects


TRangePhase = {
                "all"      :[   0.0, 1.0e10],
                "hot"      :[1.0e+5, 1.0e10],
                "warm-hot" :[4.0e+4, 1.0e+5],
                "warm-cool":[1.0e+4, 4.0e+4],
                "cool"     :[   0.0, 1.0e+4],
                "cold"     :[   0.0, 5.0e+2],
              }


DRangePhase = {
                "all"      :[4.0e-32, 1.0],
                "hot"      :[4.0e-32, 1.0],
                "warm-hot" :[4.0e-32, 1.0],
                "warm-cool":[4.0e-32, 1.0],
                "cool"     :[4.0e-32, 1.0],
                "cold"     :[4.0e-32, 1.0],
              }

def get_inoutflow_region( ds, bound_z_kpc, phase, isIn=False ):

   # define the region
   cen            = ds.domain_center
   upper_corner_L = cen - 0.2*ds.domain_width + ds.arr( [0.0, 0.0, 0.2*ds.domain_width[2].in_units('kpc').d + bound_z_kpc], 'kpc' )
   upper_corner_R = cen + 0.2*ds.domain_width + ds.arr( [0.0, 0.0, 0.3*ds.domain_width[2].in_units('kpc').d],               'kpc' )
   lower_corner_L = cen - 0.2*ds.domain_width - ds.arr( [0.0, 0.0, 0.3*ds.domain_width[2].in_units('kpc').d],               'kpc' )
   lower_corner_R = cen + 0.2*ds.domain_width - ds.arr( [0.0, 0.0, 0.2*ds.domain_width[2].in_units('kpc').d + bound_z_kpc], 'kpc' )

   if ds.dataset_type == 'gamer':

      FlowDir = "obj['gas', 'outflow_z_velocity'].d > 0.0" if not isIn else "obj['gas', 'outflow_z_velocity'].d < 0.0"
      upper_region   = ds.box( upper_corner_L, upper_corner_R ).cut_region( [FlowDir] )
      lower_region   = ds.box( lower_corner_L, lower_corner_R ).cut_region( [FlowDir] )
      GasPhase = "(obj['gas', 'T'].in_units('K').d > %.1e) & (obj['gas', 'T'].in_units('K').d < %.1e) & (obj['gas', 'density'].in_units('g/cm**3').d > %.1e)"%(TRangePhase[phase][0],TRangePhase[phase][1],DRangePhase[phase][0])
      inoutflow_region = upper_region.cut_region( [GasPhase] ) + lower_region.cut_region( [GasPhase] )
      inoutflow_p_type = 'gas'

   elif ds.dataset_type == 'gadget_hdf5':

      FlowDir = +1 if not isIn else -1
      def inoutflow_region_particle( pfilter, data ):
         filter = ( (
                       (
                         (data[ 'gas', 'x' ] > upper_corner_L[0]) & (data[ 'gas', 'x' ] < upper_corner_R[0]) &
                         (data[ 'gas', 'y' ] > upper_corner_L[1]) & (data[ 'gas', 'y' ] < upper_corner_R[1]) &
                         (data[ 'gas', 'z' ] > upper_corner_L[2]) & (data[ 'gas', 'z' ] < upper_corner_R[2]) &
                         (FlowDir*data[ 'gas', 'outflow_z_velocity' ] > 0.0)
                       ) |
                       (
                         (data[ 'gas', 'x' ] > lower_corner_L[0]) & (data[ 'gas', 'x' ] < lower_corner_R[0]) &
                         (data[ 'gas', 'y' ] > lower_corner_L[1]) & (data[ 'gas', 'y' ] < lower_corner_R[1]) &
                         (data[ 'gas', 'z' ] > lower_corner_L[2]) & (data[ 'gas', 'z' ] < lower_corner_R[2]) &
                         (FlowDir*data[ 'gas', 'outflow_z_velocity' ] > 0.0)
                       )
                    ) &
                    (
                       (data['gas', 'T'].in_units('K').d > TRangePhase[phase][0]) &
                       (data['gas', 'T'].in_units('K').d < TRangePhase[phase][1]) &
                       (data['gas', 'density'].in_units('g/cm**3').d > DRangePhase[phase][0])
                    )
                  )
         return filter
      p_type_name = 'outflow_region_particle' if not isIn else 'inflow_region_particle'
      yt.add_particle_filter( p_type_name, function=inoutflow_region_particle, filtered_type='gas' )
      ds.add_particle_filter( p_type_name )
      inoutflow_region = ds.all_data()
      inoutflow_p_type = p_type_name

   return inoutflow_region, inoutflow_p_type


def plot_inoutflow_region(ds, inoutflow_region, field_type, suffix, isIn=False):

   FlowDir = 'out' if not isIn else 'in'

   def plot_inoutflow_region_slice(direction, field, width):
      s = yt.SlicePlot( ds, direction, field, center=ds.domain_center, width=width, buff_size=(1024, 1024), data_source=inoutflow_region )
      s_frb = s.frb[field].d
      has_valid = np.any(~np.isnan(s_frb))
      if not has_valid:
         return

      s.set_axes_unit( 'kpc' )
      if field[1] == 'T':
         s.set_cmap( field, 'magma'  )
         s.set_zlim( field, 2.0e3, 2.0e5  )
      elif field[1] == 'density':
         s.set_cmap( field, 'viridis' )
         s.set_zlim( field, 1.0e-31, 1.0e-27 )
      s.set_log ( field, True )
      s.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
      if direction == 'z':
         s.annotate_quiver(('gas','velocity_x'), ('gas','velocity_y'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
      if direction == 'x':
         s.annotate_quiver(('gas','velocity_y'), ('gas','velocity_z'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
      try:
         s.save( './imgs_o/fig_%s_galactic_%sflow_global_%s_gas_Slice_%s_%s.png'%(ds, FlowDir, suffix, direction, field[1]), mpl_kwargs={'dpi':150} )
      except Exception as e:
         print( e )
         pass

   def plot_inoutflow_region_particle_projection(direction, field, width):
      if ds.dataset_type != 'gadget_hdf5':
         return
      p = yt.ParticleProjectionPlot( ds, direction, field, weight_field=(field_type,'ones'), center=ds.domain_center, width=width, depth=width )
      p.set_axes_unit( 'kpc' )
      if field[1] == 'T':
         p.set_cmap( field, 'magma' )
         p.set_zlim( field, 2.0e3, 2.0e5  )
      elif field[1] == 'density':
         p.set_cmap( field, 'viridis' )
         p.set_zlim( field, 1.0e-34, 1.0e-27 )
      p.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
      p.save( './imgs_o/fig_%s_galactic_%sflow_global_%s_gas_Projection_%s_%s.png'%(ds, FlowDir, suffix, direction, field[1]), mpl_kwargs={'dpi':150} )

   width = (6.0*(1.25**16), 'kpc')
   for direction in ['x']:
      for field in [(field_type,'density'), (field_type,'T')]:
         plot_inoutflow_region_slice(direction, field, width)
         plot_inoutflow_region_particle_projection(direction, field, width)


def plot_cold_dense_gas(ds, inoutflow_region, field_type, suffix, isIn=False):

   FlowDir = 'out' if not isIn else 'in'

   if ds.dataset_type == 'gamer':
      min_Tsqr_over_rho, *location = inoutflow_region.cut_region( ["obj['gas', 'density'].in_units('g/cm**3').d > 3.0e-31"] ).quantities.min_location( (field_type, 'Tsqr_over_rho') )
   elif ds.dataset_type == 'gadget_hdf5':
      mask = (inoutflow_region[(field_type, 'density')].in_units('g/cm**3').d > 3.0e-31)
      if len(inoutflow_region[(field_type, 'Tsqr_over_rho')][mask]) <= 0:
         return
      idx = np.arange(inoutflow_region[(field_type, 'Tsqr_over_rho')].shape[0])[mask][inoutflow_region[(field_type, 'Tsqr_over_rho')][mask].argmin()]
      Dens_par =   inoutflow_region[(field_type, 'density')][idx]
      T_par    =   inoutflow_region[(field_type, 'T')][idx]
      location = [ inoutflow_region[(field_type, 'x')][idx],
                   inoutflow_region[(field_type, 'y')][idx],
                   inoutflow_region[(field_type, 'z')][idx] ]
   Dens_point = ds.point(location)[(field_type,'density')]
   T_point    = ds.point(location)[(field_type,'T')]

   if len(T_point.d) <= 0 or len(Dens_point.d) <= 0:
      return

   annotated_str = r'$T_\mathrm{min}$ = %.1f K, $\rho$ = %.1e g cm$^{-3}$'%(T_point.in_units('K').d[0], Dens_point.in_units('g/cm**3').d[0])+\
                     '\n at [ %.2f, %.2f, %.2f ] kpc'%((location[0]-ds.domain_center[0]).in_units('kpc').d,
                                                       (location[1]-ds.domain_center[1]).in_units('kpc').d,
                                                       (location[2]-ds.domain_center[2]).in_units('kpc').d)

   def plot_cold_dens_gas_slice(direction, field, width, zoom_out):
      s = yt.SlicePlot( ds, direction, field, center=location, width=width, buff_size=(1024, 1024), data_source=inoutflow_region )
      s.set_axes_unit( 'kpc' )
      cmap = 'viridis' if field[1] == 'density' else 'magma'
      s.set_cmap( field, cmap )
      if field[1] == 'T':
         s.set_cmap( field, 'magma'  )
         s.set_zlim( field, 2.0e3, 2.0e5  )
      elif field[1] == 'density':
         s.set_cmap( field, 'viridis' )
         s.set_zlim( field, 1.0e-31, 1.0e-27 )
      s.set_log ( field, True )
      s.annotate_timestamp( time_unit='Myr', corner='upper_right' )
      if direction == 'z':
         s.annotate_quiver(('gas','velocity_x'), ('gas','velocity_y'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
      if direction == 'x':
         s.annotate_quiver(('gas','velocity_y'), ('gas','velocity_z'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
      s.annotate_text( (0.02, 0.08), annotated_str,
                       coord_system='axis', text_args={'color':'w', 'path_effects':[patheffects.withStroke(linewidth=2, foreground='k')]} )
      try:
         s.save( './imgs_o/fig_%s_galactic_%sflow_local_%s_gas_Slice_%s_%s_%02d.png'%(ds, FlowDir, suffix, direction, field[1], zoom_out), mpl_kwargs={'dpi':150} )
      except Exception as e:
         print( e )
         pass

   def plot_cold_dens_gas_particle_projection(direction, field, width, zoom_out):
      if ds.dataset_type != 'gadget_hdf5':
         return
      p = yt.ParticleProjectionPlot( ds, direction, field, weight_field=(field_type,'ones'), center=location, width=width, depth=width )
      p.set_axes_unit( 'kpc' )
      if field[1] == 'T':
         p.set_cmap( field, 'magma' )
         p.set_zlim( field, 2.0e3, 2.0e5  )
      elif field[1] == 'density':
         p.set_cmap( field, 'viridis' )
         p.set_zlim( field, 1.0e-34, 1.0e-27 )
      p.annotate_timestamp( time_unit='Myr', corner='upper_right' )
      p.annotate_text( (0.02, 0.08), annotated_str,
                       coord_system='axis', text_args={'color':'w', 'path_effects':[patheffects.withStroke(linewidth=2, foreground='k')]} )
      p.save( './imgs_o/fig_%s_galactic_%sflow_local_%s_gas_Projection_%s_%s_%02d.png'%(ds, FlowDir, suffix, direction, field[1], zoom_out), mpl_kwargs={'dpi':150} )

   for direction in ['x']:
      zoom_out_max = 1 if direction == 'x' else 0
      for field in [(field_type,'density'), (field_type,'T')]:
         for zoom_out in range(0, zoom_out_max+1, 1):
            width = (3.0*(2.0**zoom_out), 'kpc')
            plot_cold_dens_gas_slice(direction, field, width, zoom_out)
            #plot_cold_dens_gas_particle_projection(direction, field, width, zoom_out)
