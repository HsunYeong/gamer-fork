import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
from mpl_toolkits.axes_grid1 import AxesGrid


###################################################################################################
plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
colormap    = 'inferno'
dpi         = 150
###################################################################################################

def log_bin_width(centers):
   log_centers = np.log10( centers )

   log_edges = np.concatenate( ([log_centers[0]  - 0.5*(log_centers[1]   - log_centers[0] )],
                                                   0.5*(log_centers[:-1] + log_centers[1:]),
                                [log_centers[-1] + 0.5*(log_centers[-1]  - log_centers[-2])]) )

   return np.diff( log_edges )


def create_PDF( ds, obj, prefix, field_type, field_x, weighting_type, lim_x_min, lim_x_max, nbin, units_x ):

   assert weighting_type == 'energy' or weighting_type == 'mass' or weighting_type == 'dust_mass', 'Incorrect weighting type'

   prof = yt.create_profile( obj, (field_type, field_x), fields=(field_type, weighting_type),
                             weight_field=None,
                             n_bins=nbin,
                             logs={(field_type, field_x): True, (field_type, weighting_type): True},
                             units={(field_type, field_x): units_x},
                             extrema={(field_type, field_x): (lim_x_min, lim_x_max)},
                             fractional=False,
                           )
   prof.x.tofile('./bins/%s%s_%s_weighted_%s_profile_x.bin'%(ds, prefix, weighting_type, field_x.capitalize()))
   prof[(field_type, weighting_type)].tofile('./bins/%s%s_%s_weighted_%s_profile_y.bin'%(ds, prefix, weighting_type, field_x.capitalize()))

   arr_x = prof.x.in_units(units_x).d
   arr_y = prof[(field_type, weighting_type)].d

   log_bin_w = log_bin_width(arr_x)
   pdf_y = arr_y / (sum(arr_y) * log_bin_w) if sum(arr_y) > 0 else np.zeros_like(arr_y)

   np.savetxt( './pdfs/%s%s_%s_weighted_%s_PDF'%(ds, prefix, weighting_type, field_x.capitalize()),
               np.column_stack( (arr_x, pdf_y) ),
               fmt='%24.8e',
               header='%22s %24s'%( field_x.capitalize()+'_in_'+units_x, 'PDF') )


def create_PhaseDiagram( ds, obj, prefix, field_type, norm, weighting_type, nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max ):

   assert weighting_type == 'energy' or weighting_type == 'mass' or weighting_type == 'dust_mass', 'Incorrect weighting type'

   v_field = 'outflow_radial_velocity' if norm == 'r' else 'outflow_z_velocity'

   units_x = 'km/s'
   units_y = 'K'
   units_w = 'Msun' if weighting_type == 'mass' or weighting_type == 'dust_mass' else 'Msun*km**2/s**2'

   temp_vel_profile = yt.create_profile( obj,
                                         [(field_type, v_field), (field_type, 'T')],
                                         (field_type, weighting_type),
                                         logs={(field_type, v_field): True, (field_type, 'T'): True, (field_type, weighting_type): True},
                                         units={(field_type, v_field): units_x, (field_type, 'T'): units_y, (field_type, weighting_type): units_w},
                                         extrema={(field_type, v_field): (x_lim_min, x_lim_max), (field_type, 'T'): (y_lim_min, y_lim_max)},
                                         weight_field=None, n_bins=[nbin,nbin], fractional=False )

   temp_vel_profile[(field_type, weighting_type)].tofile('./bins/%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram_z.bin'%(ds, prefix, norm, weighting_type))
   temp_vel_profile.x.tofile('./bins/%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram_x.bin'%(ds, prefix, norm, weighting_type))
   temp_vel_profile.y.tofile('./bins/%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram_y.bin'%(ds, prefix, norm, weighting_type))

   create_PDF( ds, obj, prefix, field_type, v_field, weighting_type, x_lim_min, x_lim_max, nbin, units_x )
   create_PDF( ds, obj, prefix, field_type, 'T',     weighting_type, y_lim_min, y_lim_max, nbin, units_y )


def plot_PDF( fig, pos_in_fig, indices, code, prefix, field_x, weighting_type, lim_x_min, lim_x_max, nbin, units_x ):

   assert weighting_type == 'energy' or weighting_type == 'mass' or weighting_type == 'dust_mass', 'Incorrect weighting type'

   for idx in indices[::-1]:

      if code == 'GAMER':
         filename = 'Data_%06d'%(idx)
      elif code == 'GIZMO':
         filename = 'snap_%03d'%(idx)

      if idx == indices[-1]:
         arr_x  = np.fromfile('./bins/%s%s_%s_weighted_%s_profile_x.bin'%(filename, prefix, weighting_type, field_x.capitalize()))
         arr_y  = np.fromfile('./bins/%s%s_%s_weighted_%s_profile_y.bin'%(filename, prefix, weighting_type, field_x.capitalize()))
      else:
         try:
            arr_y += np.fromfile('./bins/%s%s_%s_weighted_%s_profile_y.bin'%(filename, prefix, weighting_type, field_x.capitalize()))
         except Exception as e:
            print( e )
            pass


   log_bin_w = log_bin_width(arr_x)
   pdf_y = arr_y / (sum(arr_y) * log_bin_w) if sum(arr_y) > 0 else np.zeros_like(arr_y)

   if pos_in_fig == 'x':
      rect = [0.32,0.74,0.588,0.21]
   elif pos_in_fig == 'y':
      rect = [0.03,0.1,0.21,0.6]

   ax = fig.add_axes(rect=rect)

   W = 'E' if weighting_type == 'volume' else 'M'

   if pos_in_fig == 'x':
      ax.loglog( arr_x, pdf_y )

      ax.set_title( 'PDF weighted by %s'%weighting_type )
      #ax.set_xlabel( field_x.capitalize()+' (%s)'%units_x.replace("**", "^") )
      ax.set_ylabel( r'$\frac{1}{%s}\,\frac{\mathrm{d}\,%s}{\mathrm{d}\,\log_{10}{(\mathrm{\frac{%s}{%s}})}}$'%(W, W, field_x.capitalize().replace("_", "\ "), units_x.replace("**", "^")), fontsize='x-large')

      ax.set_xlim( lim_x_min, lim_x_max )
      ax.set_ylim( 1e-8, 2e1 )

   if pos_in_fig == 'y':
      ax.loglog( pdf_y, arr_x )

      #ax.set_ylabel( field_x.capitalize()+' (%s)'%units_x.replace("**", "^"), rotation=90, labelpad=8.0 )
      ax.yaxis.tick_left()
      ax.yaxis.set_label_position("left")
      ax.set_xlabel( r'$\frac{1}{%s}\,\frac{\mathrm{d}\,%s}{\mathrm{d}\,\log_{10}{(\mathrm{\frac{%s}{%s}})}}$'%(W, W, field_x.capitalize(), units_x.replace("**", "^")), rotation=0, fontsize='x-large')

      plt.setp(ax.get_xticklabels(), rotation=0, va="top",    ha="center")
      plt.setp(ax.get_yticklabels(), rotation=0, va="center", ha="right"  )

      ax.set_ylim( lim_x_min, lim_x_max )
      ax.set_xlim( 1e-8, 2e1 )
      ax.invert_xaxis()

   # log scale
   ax.set_xscale('log')
   ax.set_yscale('log')

   minor_ticks = LogLocator( base=10.0, subs=np.arange(1.0, 10.0)*0.1, numticks=10 )
   ax.xaxis.get_ticklocs( minor=True )
   ax.yaxis.get_ticklocs( minor=True )
   ax.minorticks_on()
   ax.xaxis.set_minor_locator( minor_ticks )
   ax.xaxis.set_minor_formatter( NullFormatter() )
   ax.xaxis.set_ticks_position( 'both' )
   ax.yaxis.set_ticks_position( 'both' )
   ax.tick_params( which='both', direction='in' )

   ax.grid()


def plot_PhaseDiagram( indices, code, prefix, norm, weighting_type, nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max, title_time, fileout ):

   assert weighting_type == 'energy' or weighting_type == 'mass' or weighting_type == 'dust_mass', 'Incorrect weighting type'

   v_field = 'outflow_radial_velocity' if norm == 'r' else 'outflow_z_velocity'

   units_x = 'km/s'
   units_y = 'K'
   units_w = 'Msun' if weighting_type == 'mass' or weighting_type == 'dust_mass' else 'Msun*km**2/s**2'

#  plot
   fig = plt.figure(figsize=(10.75,9.2))

   grid = AxesGrid( fig, rect=[0.32, 0.1, 0.6, 0.6], nrows_ncols=(1, 1),
                    axes_pad=0.05, label_mode="L", share_all=True,
                    cbar_location="right", cbar_mode="single", cbar_size="2%", cbar_pad="0%",
                    aspect=False,
                  )

   for idx in indices[::-1]:

      if code == 'GAMER':
         filename = 'Data_%06d'%(idx)
      elif code == 'GIZMO':
         filename = 'snap_%03d'%(idx)

      if idx == indices[-1]:
         xf  = np.fromfile('./bins/%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram_x.bin'%(filename, prefix, norm, weighting_type))
         yf  = np.fromfile('./bins/%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram_y.bin'%(filename, prefix, norm, weighting_type))
         zf  = np.fromfile('./bins/%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram_z.bin'%(filename, prefix, norm, weighting_type)).reshape(len(xf),len(yf))
      else:
         try:
            zf += np.fromfile('./bins/%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram_z.bin'%(filename, prefix, norm, weighting_type)).reshape(len(xf),len(yf))
         except Exception as e:
            print( e )
            pass

   w_total_obj = zf.sum()
   zf = zf / w_total_obj if w_total_obj > 0 else np.zeros_like(zf)

   FONT_SIZE = 18

   im = grid[0].axes.pcolormesh(xf, yf, zf.T, shading='nearest', vmin=5.0e-10, vmax=5.0e-3, cmap=colormap, norm='log' )
   grid[0].axes.set_xlim( x_lim_min, x_lim_max )
   grid[0].axes.set_ylim( y_lim_min, y_lim_max )
   grid[0].axes.set_xscale('log')
   grid[0].axes.set_yscale('log')
   grid[0].axes.set_xlabel(r'Outflow %s Velocity  $(\mathrm{km\ s^{-1}})$'%norm, fontsize=FONT_SIZE)
   grid[0].axes.set_ylabel('T  (K)', fontsize=FONT_SIZE)
   grid[0].axes.axes.tick_params( which='both', left=True, right=True, bottom=True, top=True )
   grid[0].axes.tick_params( which='both', direction='in', labelsize=FONT_SIZE )
   cbar = fig.colorbar(im, ax=grid[0].axes, cax=grid.cbar_axes[0], pad=0.0)
   cbar.set_label(weighting_type.capitalize()+' Fraction', fontsize=FONT_SIZE)
   cbar.ax.tick_params( which='both', direction='in', labelsize=FONT_SIZE )

   grid[0].axes.annotate( title_time,
                          (x_lim_min*10**(0.70*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                           y_lim_min*10**(0.90*(np.log10(y_lim_max)-np.log10(y_lim_min)))),
                          color='black', fontsize=FONT_SIZE )
   grid[0].axes.annotate( 'Total {:s} = {:.2e} {:s}'.format( weighting_type, w_total_obj, units_w.replace("**", "^") ),
                          (x_lim_min*10**(0.70*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                           y_lim_min*10**(0.85*(np.log10(y_lim_max)-np.log10(y_lim_min)))),
                          color='black', fontsize='small' )

   plot_PDF( fig, 'x', indices, code, prefix, v_field, weighting_type, x_lim_min, x_lim_max, nbin, units_x )
   plot_PDF( fig, 'y', indices, code, prefix, 'T',     weighting_type, y_lim_min, y_lim_max, nbin, units_y )
   fig.savefig( './imgs_tv/fig_%s%s_Temp-Vel%s_%s_weighted_PhaseDiagram.png'%(fileout, prefix, norm, weighting_type) )
   plt.close()
