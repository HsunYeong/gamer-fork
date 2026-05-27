import yt

def plot_outflow_patch_projection(ds):

   width_x =  1.0
   width_y =  5.0
   width_z = 12.5

   center_patch = ds.domain_center + ds.arr([0.0, 0.0*width_y, 0.5*width_z], 'kpc')
   width_patch  = ((width_y, 'kpc'),(width_z, 'kpc'))

   box  = ds.box( center_patch-0.55*ds.arr([width_x, width_y, width_z], 'kpc'),
                  center_patch+0.55*ds.arr([width_x, width_y, width_z], 'kpc') )

   p = yt.ProjectionPlot( ds, 'x', 'density', data_source=box, center=center_patch, origin='center-domain',
                          width=width_patch, method='integrate' , weight_field=('index','ones'), buff_size=(1024, 1024) )
   p.set_axes_unit( 'kpc' )
   p.set_unit( 'density', 'g/cm**3' )
   p.set_zlim( 'density', 2.0e-30, 2.0e-27 )
   p.set_cmap( 'density', 'viridis' )
   p.annotate_quiver(('gas','velocity_y'), ('gas','velocity_z'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
   p.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   p.annotate_grids( periodic=False )
   p.save( 'fig_%s_galactic_outflow_patch_Projection_x_density.png'%(ds), mpl_kwargs={'dpi':150} )

for i in range(70, 75+1, 1):
   ds = yt.load('../Data_%06d'%i)

   plot_outflow_patch_projection(ds)
