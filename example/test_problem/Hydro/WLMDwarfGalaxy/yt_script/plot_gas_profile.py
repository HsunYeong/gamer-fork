import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import WLMDwarfGalaxy_load_datasets
import WLMDwarfGalaxy_derived_fields

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

disk_normal = [0.0, 0.0, 1.0]
width_kpc   = 9
nbin        = 50
markersize  = 4.0
dpi         = 150


yt.enable_parallelism()

# load the dataset
ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)

WLMDwarfGalaxy_derived_fields.set_particle_types(code)

yt_radius   = ('index', 'cylindrical_radius') if code == 'GAMER' else ('gas', 'cylindrical_radius')
yt_theta    = ('index', 'cylindrical_theta')  if code == 'GAMER' else ('gas', 'cylindrical_theta')
yt_tan_vel  = 'velocity_cylindrical_theta'
yt_type     = 'cell' if code == 'GAMER' else 'particle'
yt_mass     = ('gas', 'mass')

# loop over all datasets
for ds in ts.piter():

   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds)

   cen     = ds.domain_center

#  only include the data within a sphere with a radius of 0.5*width_kpc
   sp_gas  = ds.sphere( cen, (0.5*width_kpc, 'kpc') ).cut_region( ["obj['gas', 'density'].in_units('g/cm**3') > 1.0e-30"] )
   sp_gas.set_field_parameter( 'normal', disk_normal )
   sp_disk = ds.sphere( cen, (0.5*width_kpc, 'kpc') )
   sp_disk.set_field_parameter( 'normal', disk_normal )
   sp_halo = ds.sphere( cen, (90.0, 'kpc') )


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


#  (2) gas temperature
   prof     = yt.ProfilePlot( sp_gas, yt_radius, ('gas', 'T'), weight_field='density',
                              n_bins=nbin, x_log=False, accumulation=False )
   gas_temp = prof.profiles[0]['T'].in_units('K').d


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


#  (5) Disk surface density
   Disk_prof     = yt.ProfilePlot( sp_disk, ('Disk', 'particle_position_cylindrical_radius'), ('Disk','particle_mass'), weight_field=None,
                                   n_bins=nbin, x_log=False, accumulation=False )
   Disk_dens     = Disk_prof.profiles[0][('Disk','particle_mass')].in_units('Msun').d
   Disk_radius   = Disk_prof.profiles[0].x.in_units('kpc').d

#  convert mass to surface density in Msun/pc^2
   Disk_dr = Disk_radius[1] - Disk_radius[0] # assuming linear bin
   for b in range( len(Disk_radius) ):
      area         = np.pi*( (Disk_radius[b]+0.5*Disk_dr)**2 - (Disk_radius[b]-0.5*Disk_dr)**2 )
      Disk_dens[b] /= area*1.0e6


#  (6) Disk rotational velocity
   Disk_prof  = yt.ProfilePlot( sp_disk, ('Disk', 'particle_position_cylindrical_radius'), ('Disk', 'particle_velocity_cylindrical_theta'),
                                weight_field=('Disk','particle_mass'), n_bins=nbin, x_log=False )
   Disk_vrot = Disk_prof.profiles[0][('Disk', 'particle_velocity_cylindrical_theta')].in_units('km/s').d


#  (7) Halo density
   Halo_prof     = yt.ProfilePlot( sp_halo, ('Halo', 'particle_position_spherical_radius'), ('Halo','particle_mass'), weight_field=None,
                                   n_bins=nbin, x_log=False, accumulation=False )
   Halo_dens     = Halo_prof.profiles[0][('Halo','particle_mass')].in_units('Msun').d
   Halo_radius   = Halo_prof.profiles[0].x.in_units('kpc').d

#  convert mass to volume density in Msun/pc^3
   Halo_dr = Halo_radius[1] - Halo_radius[0] # assuming linear bin
   for b in range( len(Halo_radius) ):
      volume       = 4.0/3.0*np.pi*( (Halo_radius[b]+0.5*Halo_dr)**3 - (Halo_radius[b]-0.5*Halo_dr)**3 )
      Halo_dens[b] /= volume*1.0e9


#  (8) Halo velocity
   Halo_prof  = yt.ProfilePlot( sp_halo, ('Halo', 'particle_position_spherical_radius'), ('Halo', 'particle_velocity_magnitude'),
                                weight_field=('Halo','particle_mass'), n_bins=nbin, x_log=False )
   Halo_vtot = Halo_prof.profiles[0][('Halo', 'particle_velocity_magnitude')].in_units('km/s').d


#  plot
#  ==================================================================================
   f, ax = plt.subplots( 2, 4, figsize=(12.8, 4.8) )
   f.subplots_adjust( wspace=0.4 )
   markersize  = 4.0

   ax[1][0].set_xlabel( '$\mathrm{Cylindrical\ radius\ [kpc]}$', fontsize='large' )
   ax[1][1].set_xlabel( '$\mathrm{Cylindrical\ radius\ [kpc]}$', fontsize='large' )
   ax[1][2].set_xlabel( '$\mathrm{Cylindrical\ radius\ [kpc]}$', fontsize='large' )
   ax[1][3].set_xlabel( '$\mathrm{Spherical\ radius\ [kpc]}$',   fontsize='large' )

#  (1) gas surface density
   ax[0][0].plot( gas_radius,  gas_dens,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[0][0].set_yscale( 'log', nonpositive='clip' )
   ax[0][0].set_xlim(    0.0,   7.0 )
   ax[0][0].set_ylim( 1.0e-2, 3.0e1 )
   ax[0][0].set_ylabel( '$\mathrm{\Sigma_{gas}\ [M_{\odot}/pc^2]}$', fontsize='large' )
   #ax[0][0].legend()

#  (2) gas temperature
   ax[0][1].plot( gas_radius,  gas_temp,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[0][1].set_yscale( 'log', nonpositive='clip' )
   ax[0][1].set_xlim(    0.0,  7.0 )
   ax[0][1].set_ylim( 1.0e3, 1.0e5 )
   ax[0][1].set_ylabel( '$\mathrm{T\ [K]}$', fontsize='large' )

#  (3) gas rotational velocity
   ax[1][0].plot( gas_radius,  gas_vrot,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[1][0].set_xlim( 0.0,  7.0 )
   ax[1][0].set_ylim( 0.0, 60.0 )
   ax[1][0].set_ylabel( '$\mathrm{v_{rot,gas}\ [km/s]}$', fontsize='large' )

#  (4) gas velocity dispersion
   ax[1][1].plot( gas_radius,  gas_vdis,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[1][1].set_xlim( 0.0,  7.0 )
   ax[1][1].set_ylim( 0.0,  5.0 )
   ax[1][1].set_ylabel( '$\mathrm{\sigma_{gas}\ [km/s]}$', fontsize='large' )

#  (5) Disk surface density
   ax[0][2].plot( Disk_radius,  Disk_dens,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[0][2].set_yscale( 'log', nonpositive='clip' )
   ax[0][2].set_xlim(    0.0,   7.0 )
   ax[0][2].set_ylim( 1.0e-3, 3.0e0 )
   ax[0][2].set_ylabel( '$\mathrm{\Sigma_{Disk}\ [M_{\odot}/pc^2]}$', fontsize='large' )

#  (6) Disk rotational velocity
   ax[1][2].plot( Disk_radius,  Disk_vrot,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[1][2].set_xlim( 0.0,  7.0 )
   ax[1][2].set_ylim( 0.0, 50.0 )
   ax[1][2].set_ylabel( '$\mathrm{v_{rot,Disk}\ [km/s]}$', fontsize='large' )

#  (7) Halo density
   ax[0][3].plot( Halo_radius,  Halo_dens,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[0][3].set_yscale( 'log', nonpositive='clip' )
   ax[0][3].set_xlim(    0.0,   90.0 )
   ax[0][3].set_ylim( 1.0e-7, 1.0e-1 )
   ax[0][3].set_ylabel( r'$\mathrm{\rho_{Halo}\ [M_{\odot}/pc^3]}$', fontsize='large' )

#  (8) Halo velocity
   ax[1][3].plot( Halo_radius,  Halo_vtot,  'r-s',  lw=1, mec='none', ms=markersize, label=''  )
   ax[1][3].set_xlim( 0.0, 90.0 )
   ax[1][3].set_ylim( 0.0, 55.0 )
   ax[1][3].set_ylabel( '$\mathrm{v_{Halo}\ [km/s]}$', fontsize='large' )

#  add title
   time = ds.current_time.in_units('Myr')
   f.suptitle( "t = %6.2f %s"%(time.d, time.units), fontsize='large' )

#  show/save figure
   f.savefig( 'fig_%s_gas_profile.png'%ds, bbox_inches='tight', pad_inches=0.05, dpi=dpi )

