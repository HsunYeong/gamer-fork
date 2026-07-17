import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import WLMDwarfGalaxy_load_datasets
import WLMDwarfGalaxy_derived_fields
from matplotlib.pyplot import cm
from matplotlib.ticker import LogLocator, NullFormatter
from matplotlib.patches import Rectangle

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot various gas profiles' )

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


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
code        = args.code

plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
dpi         = 150
FONT_SIZE=24
TICK_SIZE=20

disk_normal = [0.0, 0.0, 1.0]
width_kpc   = 9
nbin        = 300
radius      = 8.0
height      = 3.0

markersize  = 4.0


yt.enable_parallelism()

# load the dataset
ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)

WLMDwarfGalaxy_derived_fields.set_particle_types(code)

yt_radius   = ('index', 'cylindrical_radius')
yt_theta    = ('index', 'cylindrical_theta')
yt_tan_vel  = 'velocity_cylindrical_theta'
yt_type     = 'cell'
yt_mass     = ('gas', 'cell_mass')
yt_vol      = ('gas', 'cell_volume')

my_storage = {}

# loop over all datasets
for sto, ds in ts.piter(storage=my_storage):

   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds)
   def _abs_z( field, data ):
      return (data[('gas', 'z')]**2)**0.5
   ds.add_field( ("gas", "abs_z"), function=_abs_z, sampling_type="cell", units="kpc" )
   def _Br_square( field, data ):
      return data[('gas', 'magnetic_field_cylindrical_radius')]**2
   ds.add_field( ("gas", "Br_square"), function=_Br_square, sampling_type="cell", units="G**2" )
   def _Bp_square( field, data ):
      return data[('gas', 'magnetic_field_cylindrical_radius')]**2
   ds.add_field( ("gas", "Bp_square"), function=_Bp_square, sampling_type="cell", units="G**2" )
   def _Bz_square( field, data ):
      return data[('gas', 'magnetic_field_cylindrical_z')]**2
   ds.add_field( ("gas", "Bz_square"), function=_Bz_square, sampling_type="cell", units="G**2" )

   cen     = ds.domain_center
#  only include the data within a sphere with a radius of 0.5*width_kpc
   sp_gas = ds.disk( center='c', normal=disk_normal, radius=(radius, 'kpc'), height=(height, 'kpc') )
   '''
   sp_gas  = ds.sphere( cen, (0.5*width_kpc, 'kpc') ).cut_region( ["obj['gas', 'density'].in_units('g/cm**3') > 1.0e-30"] )
   sp_gas.set_field_parameter( 'normal', disk_normal )
   sp_disk = ds.sphere( cen, (0.5*width_kpc, 'kpc') )
   sp_disk.set_field_parameter( 'normal', disk_normal )
   sp_halo = ds.sphere( cen, (90.0, 'kpc') )
   '''

   idx  = ds.parameters["DumpID"]
   time = np.round(ds.quan( ds.parameters["Time"][0], 'code_time' ).in_units('Myr').d)

#  (1) gas surface density
   prof       = yt.ProfilePlot( sp_gas, yt_radius, yt_mass, weight_field=None,
                                n_bins=nbin, x_log=False, accumulation=False )
   gas_dens   = prof.profiles[0][yt_mass].in_units('Msun').d
   gas_radius = prof.profiles[0].x.in_units('kpc').d

#  convert mass to surface density in Msun/pc^2
   dr = gas_radius[1] - gas_radius[0] # assuming linear bin
   for b in range( len(gas_radius) ):
      area         = np.pi*( (gas_radius[b]+0.5*dr)**2 - (gas_radius[b]-0.5*dr)**2 )
      gas_dens[b] /= area*1.0e6


#  (2) gas energy
   prof     = yt.ProfilePlot( sp_gas, yt_radius,
                              [ ('gas', 'magnetic_energy_density'),
                                ('gas', 'kinetic_energy_density'),
                                ('gas', 'dual_internal_energy_density') ],
                              weight_field=yt_vol,
                              n_bins=nbin, x_log=False, accumulation=False )
   Emag = prof.profiles[0][('gas', 'magnetic_energy_density')]
   Ekin = prof.profiles[0][('gas', 'kinetic_energy_density')]
   Eint = prof.profiles[0][('gas', 'dual_internal_energy_density')]
   beta  = (Eint/Emag).d
   betak = (Ekin/Emag).d
   meanB = (Emag.in_units('G**2').d*8*np.pi)**0.5

#  (3) gas rotational velocity
#  consider only dense enough gas in order to exclude the gaseous halo
#  --> follow the AGORA analysis script: https://bitbucket.org/mornkr/agora-analysis-script/
   prof     = yt.ProfilePlot( sp_gas, yt_radius,  ('gas', yt_tan_vel),
                              weight_field=yt_mass, n_bins=nbin, x_log=False )
   gas_vrot = prof.profiles[0][yt_tan_vel].in_units('km/s').d


#  (4) gas velocity dispersion
#  --> follow the AGORA analysis script: https://bitbucket.org/mornkr/agora-analysis-script/
   def _local_rotational_velocity_x( field, data ):
      vx = np.zeros( data[('gas', 'velocity_x')].shape )
      for r, vrot in zip(gas_radius, gas_vrot):
         idx = np.where( (data[yt_radius].in_units('kpc') >= (r - 0.5*dr)) &
                         (data[yt_radius].in_units('kpc') <  (r + 0.5*dr)) )
         vx[idx] = -np.sin( data[yt_theta][idx] ) * vrot
      return data.ds.arr( vx, 'km/s' ).in_base( data.ds.unit_system.name )
   ds.add_field( ('gas', 'local_rotational_velocity_x'), function=_local_rotational_velocity_x,
                 sampling_type=yt_type, take_log=False, units='km/s' )

   def _local_rotational_velocity_y( field, data ):
      vy = np.zeros( data[('gas', 'velocity_y')].shape )
      for r, vrot in zip(gas_radius, gas_vrot):
         idx = np.where( (data[yt_radius].in_units('kpc') >= (r - 0.5*dr)) &
                         (data[yt_radius].in_units('kpc') <  (r + 0.5*dr)) )
         vy[idx] =  np.cos( data[yt_theta][idx] ) * vrot
      return data.ds.arr( vy, 'km/s' ).in_base( data.ds.unit_system.name )
   ds.add_field( ('gas', 'local_rotational_velocity_y'), function=_local_rotational_velocity_y,
                 sampling_type=yt_type, take_log=False, units='km/s' )

   def _velocity_minus_local_rotational_velocity_squared( field, data ):
      return ( data[('gas', 'velocity_x')] - data[('gas', 'local_rotational_velocity_x')] )**2 + \
             ( data[('gas', 'velocity_y')] - data[('gas', 'local_rotational_velocity_y')] )**2 + \
             ( data[('gas', 'velocity_z')]                                                )**2
   ds.add_field( ('gas', 'velocity_minus_local_rotational_velocity_squared'), function=_velocity_minus_local_rotational_velocity_squared,
                 sampling_type=yt_type, take_log=False, units='km**2/s**2' )

   prof     = yt.ProfilePlot( sp_gas, yt_radius,  ('gas', 'velocity_minus_local_rotational_velocity_squared'),
                              weight_field=yt_mass, n_bins=nbin, x_log=False )
   gas_vdis = np.sqrt( prof.profiles[0]['velocity_minus_local_rotational_velocity_squared'] ).in_units('km/s').d

#  (2) magnetic field
   prof     = yt.ProfilePlot( sp_gas, yt_radius,
                              [ ('gas', 'Br_square'),
                                ('gas', 'Bp_square'),
                                ('gas', 'Bz_square') ],
                              weight_field=yt_vol,
                              n_bins=nbin, x_log=False, accumulation=False )
   Br2 = prof.profiles[0][('gas', 'Br_square')].in_units('G**2').d
   Bp2 = prof.profiles[0][('gas', 'Bp_square')].in_units('G**2').d
   Bz2 = prof.profiles[0][('gas', 'Bz_square')].in_units('G**2').d

   sto.result = {
        "idx"     : idx,
        "time"    : time,
        "mean_r"  : gas_radius,
        "mean_vp" : gas_vrot,
        "sigma_v" : gas_vdis,
        "surdens" : gas_dens,
        "beta"    : beta,
        "betak"   : betak,
        "mean_B"  : meanB,
        "mean_Br2": Br2,
        "mean_Bp2": Bp2,
        "mean_Bz2": Bz2
   }

if yt.is_root():
   results = list(my_storage.values())
   def plot_fig(field, title, ylabel, ymin, ymax, ylog):
      plt.figure(figsize=(8.5,7.0), dpi = dpi)
      for r in results:
         idx  = r["idx"]
         time = r["time"]
         color = cm.Blues(0.3+0.6*(idx-idx_start)/(idx_end-idx_start)) if idx_start != idx_end else cm.Blues(0.6)
         plt.plot(r["mean_r"], r["%s"%field], color = color, label ='t=%3d Myr'%(np.round(time)))
         np.save('Data_Disk_%06d_%s'%(idx, title), np.array([r["mean_r"], r["%s"%field]]))

      plt.xlim((0, radius))
      if ymin != None and ymax != None:
         plt.ylim((ymin, ymax))
      if ylog:
         plt.yscale('log')
      plt.grid(ls='--')
      plt.legend(bbox_to_anchor=(1.3, 0), loc='lower right', borderaxespad=0, shadow=True, fontsize=0.8*TICK_SIZE)
      plt.xlabel(r"$R$ (kpc)", fontsize=FONT_SIZE)
      plt.ylabel( ylabel, fontsize=FONT_SIZE )
      plt.xticks(fontsize=TICK_SIZE)
      plt.yticks(fontsize=TICK_SIZE)
      plt.savefig("fig__%s.png"%title, dpi = dpi, bbox_inches="tight", pad_inches=0.05)
      plt.close()

   plot_fig( 'mean_vp',  'rotation_curve',  r"$v_{\rm cir}$ (km/s)",                        0,   80, False )
   plot_fig( 'surdens',  'surface_density', r"$\Sigma$ (${\rm M}_{\odot}/{\rm kpc}^2$)",  1e2,  1e8, True  )
   plot_fig( 'sigma_v',  'sigma_v',         r"$\sigma$ (${\rm km}/{\rm s}$)",               0,   30, False )
   plot_fig( 'beta',     'plasma_beta',     r"$\beta$",                                  None, None, True  )
   plot_fig( 'betak',    'plasma_betak',    r"$\beta_k$",                                None, None, True  )
   plot_fig( 'mean_B',   'mean_B',          r"$\langle B \rangle$ (G)",                  None, None, True  )

   fig = plt.figure(figsize=(8.5,7.0), dpi = dpi)
   axs = fig.add_subplot(111)
   label_rows = [ ['time (Myr)'  , r'$B_r$ (G)', r'$B_\phi$ (G)', r'$B_z$ (G)'] ]
   extra = Rectangle((0, 0), 1, 1.0, fc="w", fill=False, edgecolor='none', linewidth=0)
   legend_handle = [ [extra, extra, extra, extra] ]

   for r in results:
      idx  = r["idx"]
      time = r["time"]
      Br2  = r["mean_Br2"]
      Bp2  = r["mean_Bp2"]
      Bz2  = r["mean_Bz2"]

      color1 = cm.Blues (0.3+0.6*(idx-idx_start)/(idx_end-idx_start)) if idx_start != idx_end else cm.Blues (0.6)
      color2 = cm.Reds  (0.3+0.6*(idx-idx_start)/(idx_end-idx_start)) if idx_start != idx_end else cm.Reds  (0.6)
      color3 = cm.Greens(0.3+0.6*(idx-idx_start)/(idx_end-idx_start)) if idx_start != idx_end else cm.Greens(0.6)

      l1, = axs.plot(r["mean_r"], Br2, color = color1, label =r'$B_r^2$'   )
      l2, = axs.plot(r["mean_r"], Bp2, color = color2, label =r'$B_\phi^2$')
      l3, = axs.plot(r["mean_r"], Bz2, color = color3, label =r'$B_z^2$'   )
      label_rows.append( ['%3d'%(np.round(time)), '', '', ''])
      legend_handle.append( [extra, l1, l2, l3] )
   legend_labels = np.array( label_rows ).flatten('F')
   legend_handle  = np.array(legend_handle).flatten('F')
   axs.legend(legend_handle, legend_labels, loc='lower left',fontsize=9, ncol = 4, handletextpad = -2.5, handlelength=2.5, handleheight=1.5, columnspacing=0.6, labelspacing=0.4)
   axs.set_xlim((0, radius))
   axs.set_yscale('log')
   axs.set_xlabel(r"$R$ (kpc)", fontsize=FONT_SIZE)
   axs.set_ylabel(r"$B^2$ (G$^2$)", fontsize=FONT_SIZE )
   plt.xticks(fontsize=TICK_SIZE)
   plt.yticks(fontsize=TICK_SIZE)
   plt.savefig("fig__B2.png", dpi = dpi, bbox_inches="tight", pad_inches=0.05)
   plt.close()


