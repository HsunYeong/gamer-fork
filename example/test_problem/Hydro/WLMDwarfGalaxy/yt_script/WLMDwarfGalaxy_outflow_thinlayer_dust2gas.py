import yt

def plot_outflow_thinlayer_projection_dust2gas(ds):

   width_x =  1.0
   width_y =  5.0
   width_z = 12.5

   center_thinlayer = ds.domain_center + ds.arr([0.0, 0.0, 0.5*width_z], 'kpc')
   width_thinlayer  = ((width_y, 'kpc'),(width_z, 'kpc'))

   if ds.dataset_type == 'gamer':

      box  = ds.box( center_thinlayer-0.5*ds.arr([1.0*width_x, 1.1*width_y, 1.1*width_z], 'kpc'),
                     center_thinlayer+0.5*ds.arr([1.0*width_x, 1.1*width_y, 1.1*width_z], 'kpc') )

      field = ('gas','dust2gas')

      p = yt.ProjectionPlot( ds, 'x', field,
                             data_source=box,
                             center=center_thinlayer,
                             origin='center-domain',
                             width=width_thinlayer,
                             method='integrate',
                             weight_field=('gas','density'),
                             buff_size=(1024, 1024) )

   p.set_axes_unit( 'kpc' )
   p.set_zlim( field, 5.0e-4, 2.0e-3 )
   p.set_cmap( field, 'viridis' )
   p.set_colorbar_label(field, r"$\Sigma_{\rm d}/\Sigma_{\rm g}$")
   p.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   p.save( './imgs_tl/fig_%s_galactic_outflow_thinlayer_Projection_x_dust2gas.png'%(ds), mpl_kwargs={'dpi':150} )

   if ds.dataset_type == 'gamer':
      p.annotate_grids( periodic=False )

   elif ds.dataset_type == 'gadget_hdf5':
      mask = ds.all_data()[('outflow_thinlayer_particle','particle_position_z')].in_units('kpc').d > 2.5
      for px, pr in zip(ds.all_data()[('outflow_thinlayer_particle','particle_position')][mask].in_units('kpc'),ds.all_data()[('outflow_thinlayer_particle','smoothing_length')][mask].in_units('kpc')):
         p.annotate_sphere(center=px, radius=pr/25.0, circle_args={'color':'grey','alpha':0.7})

   #p.save( './imgs_tl/fig_%s_galactic_outflow_thinlayer_Projection_x_dust2gas_withgrids.png'%(ds), mpl_kwargs={'dpi':150} )
