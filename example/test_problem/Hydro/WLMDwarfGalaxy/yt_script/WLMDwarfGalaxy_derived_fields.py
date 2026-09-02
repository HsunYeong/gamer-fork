import yt
import numpy as np


def set_particle_types(code):
   if code == 'GAMER':
      # the halo particles
      def Halo( pfilter, data ):
         filter = data[ 'all', 'ParType' ] == 2
         return filter
      yt.add_particle_filter( 'Halo', function=Halo, filtered_type='all', requires=['ParType'] )

      # the disk particles
      def Disk( pfilter, data ):
         filter = data[ 'all', 'ParType' ] == 3
         return filter
      yt.add_particle_filter( 'Disk', function=Disk, filtered_type='all', requires=['ParType'] )

      # the newly formed stars
      def new_star( pfilter, data ):
         filter = data[ 'all', 'ParCreTime' ] > 0
         return filter
      yt.add_particle_filter( 'new_star', function=new_star, filtered_type='all', requires=['ParCreTime'] )

      # the exploded SNe
      def exp_SNII( pfilter, data ):
         filter = data[ 'all', 'ParSNIITime' ] <= 0
         return filter
      yt.add_particle_filter( 'exp_SNII', function=exp_SNII, filtered_type='all', requires=['ParSNIITime'] )

      def young_star( pfilter, data ):
         filter = (data[ 'new_star', 'ParCreTime' ] > data.ds.current_time - data.ds.quan(2.5, 'Myr'))
         return filter
      yt.add_particle_filter( 'young_star', function=young_star, filtered_type='new_star', requires=['ParCreTime'] )

      def young_SNII( pfilter, data ):
         filter = (-1.0*data[ 'exp_SNII', 'ParSNIITime' ]*data.ds.units.code_time > data.ds.current_time - data.ds.quan(2.5, 'Myr'))
         return filter
      yt.add_particle_filter( 'young_SNII', function=young_SNII, filtered_type='exp_SNII', requires=['ParSNIITime'] )

   elif code == 'GIZMO':
      def Halo( pfilter, data ):
         filter = data[ 'PartType1', 'ParticleIDs' ] > 0
         return filter
      yt.add_particle_filter( 'Halo', function=Halo, filtered_type='PartType1' )

      def Disk( pfilter, data ):
         filter = data[ 'PartType2', 'ParticleIDs' ] > 0
         return filter
      yt.add_particle_filter( 'Disk', function=Disk, filtered_type='PartType2' )

      def new_star( pfilter, data ):
         filter = data[ 'PartType4', 'ParticleIDs' ] > 0
         return filter
      yt.add_particle_filter( 'new_star', function=new_star, filtered_type='PartType4' )

      def young_star( pfilter, data ):
         filter = (data[ 'new_star', 'StellarFormationTime' ] > data.ds.current_time - data.ds.quan(2.5, 'Myr'))
         return filter
      yt.add_particle_filter( 'young_star', function=young_star, filtered_type='new_star', requires=['StellarFormationTime'] )


def set_derived_fields(ds):

   sampling_type = 'cell' if ds.dataset_type == 'gamer' else 'particle'
   gamma = 5.0/3.0

   # basic fields for hydrodynamics
   if ds.dataset_type == 'gadget_hdf5':

      def _specific_thermal_energy( field, data ):
         return data[('PartType0', 'specific_thermal_energy')]
      ds.add_field( ('gas', 'specific_thermal_energy'), function=_specific_thermal_energy, sampling_type=sampling_type, units='km**2/s**2' )

      def _specific_kinetic_energy( field, data ):
         return 0.5*data[('gas', 'velocity_magnitude')]**2
      ds.add_field( ('gas', 'specific_kinetic_energy'), function=_specific_kinetic_energy, sampling_type=sampling_type, units='km**2/s**2' )

      def _specific_total_energy( field, data ):
         return 0.5*data[('gas', 'velocity_magnitude')]**2
      ds.add_field( ('gas', 'specific_total_energy'), function=_specific_total_energy, sampling_type=sampling_type, units='km**2/s**2' )

      def _pressure( field, data ):
         return (gamma-1.0) * data[('gas', 'specific_thermal_energy')] * data[('gas', 'density')]
      ds.add_field( ('gas', 'pressure'), function=_pressure, sampling_type=sampling_type, units='erg/cm**3' )

      def _volume( field, data ):
         return data[('gas', 'mass')] / data[('gas', 'density')]
      ds.add_field( ('gas', 'volume'), function=_volume, sampling_type=sampling_type, units='pc**3' )

   def _energy( field, data ):
      return data[('gas', 'mass')] * data[('gas', 'specific_total_energy')]
   ds.add_field( ('gas', 'energy'), function=_energy, sampling_type=sampling_type, units='Msun*km**2/s**2' )

   if ('gamer', 'Dust') in ds.derived_field_list:
      def _dust( field, data ):
         return data[('gamer', 'Dust')] * data.ds.quan(1, 'code_density')
      ds.add_field( ('gas', 'dust'), function=_dust, sampling_type=sampling_type, units='g/cm**3' )

      def _dust_mass( field, data ):
         return data[('gas', 'dust')] * data[('gas', 'volume')]
      ds.add_field( ('gas', 'dust_mass'), function=_dust_mass, sampling_type=sampling_type, units='g' )

      def _dust2gas( field, data ):
         return data[('gas', 'dust')] / data[('gas', 'density')]
      ds.add_field( ('gas', 'dust2gas'), function=_dust2gas, sampling_type=sampling_type, units='dimensionless' )

   # define density^2 for calculating the weighted temperature
   def _density_square( field, data ):
      return data[('gas','density')]**2
   ds.add_field( ('gas', 'density_square'), function=_density_square, sampling_type=sampling_type, units='g**2/cm**6' )

   def _resolution_size( field, data ):
      return data[('gas','volume')]**(1./3.)
   ds.add_field( ('gas', 'resolution_size'), function=_resolution_size, sampling_type=sampling_type, units='pc' )


   def _dual_internal_energy_density( field, data ):
      Enth = data[('gas', 'total_energy_density')] - data[('gas', 'thermal_energy_density')]
      mask = data[('gas', 'thermal_energy_density')]/Enth > 2.0e-1
      dual_pres = data['Dual']*(data['Dens'].in_units('code_density').d)**(gamma - 1.0)
      return np.where( mask, data[('gas', 'thermal_energy_density')].in_units('code_pressure').d, dual_pres / (gamma - 1.0) ) * data.ds.quan(1, 'code_pressure')
   ds.add_field( ('gas', 'dual_internal_energy_density'), function=_dual_internal_energy_density, sampling_type=sampling_type, units='erg/cm**3' )

   def _dual_specific_thermal_energy( field, data ):
      return data[('gas', 'dual_internal_energy_density')]/data[('gas', 'density')]
   ds.add_field( ('gas', 'dual_specific_thermal_energy'), function=_dual_specific_thermal_energy, sampling_type=sampling_type, units='km**2/s**2' )

   # re-define temperature with a fixed mean molecular weight
   field_u = ('gas', 'dual_specific_thermal_energy') if ds.dataset_type == 'gamer' else ('PartType0', 'InternalEnergy')
   ds.mu   = 0.588235294117647  # Fully-ionized, 1/( 2.0*0.76/1 + 3.0*(1-0.76)/4 )

   def _T( field, data ):
      return data[field_u] * (gamma - 1.0) * data.ds.mu * data.ds.units.proton_mass / data.ds.units.boltzmann_constant
   ds.add_field( ('gas', 'T'), function=_T, sampling_type=sampling_type, units='K' )


   def _specific_entropy( field, data ):
      return data[('gas', 'pressure')] * ( data[('gas', 'density')]/(data.ds.mu * data.ds.units.proton_mass) )**(-5.0/3.0)
   ds.add_field( ('gas', 'specific_entropy'), function=_specific_entropy, sampling_type=sampling_type, units='erg*cm**2' )


   # outflow-related fields
   def _outflow_radius( field, data ):
      return np.sqrt( (data[('gas','x')]-data.ds.domain_center[0])**2 + (data[('gas','y')]-data.ds.domain_center[1])**2 + (data[('gas','z')]-data.ds.domain_center[2])**2 )
   ds.add_field( ('gas', 'outflow_radius'), function=_outflow_radius, sampling_type=sampling_type, units='kpc' )

   def _outflow_radial_velocity( field, data ):
      return ( data[('gas','velocity_x')]*(data[('gas','x')]-data.ds.domain_center[0]) + data[('gas','velocity_y')]*(data[('gas','y')]-data.ds.domain_center[1]) + data[('gas','velocity_z')]*(data[('gas','z')]-data.ds.domain_center[2]) )/data[('gas','outflow_radius')]
   ds.add_field( ('gas', 'outflow_radial_velocity'), function=_outflow_radial_velocity, sampling_type=sampling_type, units='km/s' )

   def _outflow_z( field, data ):
      return np.abs( data[('gas','z')] - data.ds.domain_center[2] )
   ds.add_field( ('gas', 'outflow_z'), function=_outflow_z, sampling_type=sampling_type, units='kpc' )

   def _outflow_z_velocity( field, data ):
      return data[('gas','velocity_z')] * np.sign( data[('gas','z')] - data.ds.domain_center[2] )
   ds.add_field( ('gas', 'outflow_z_velocity'), function=_outflow_z_velocity, sampling_type=sampling_type, units='km/s' )

   def _cell_mass_radial_outflow_flux( field, data ):
      return data[('gas', 'mass')] * data[('gas', 'outflow_radial_velocity')]
   ds.add_field( ('gas', 'cell_mass_radial_outflow_flux'), function=_cell_mass_radial_outflow_flux, sampling_type=sampling_type, units='Msun*km/s' )

   def _cell_energy_radial_outflow_flux( field, data ):
      return data[('gas', 'mass')]*( 0.5*data[('gas', 'velocity_magnitude')]**2 + gamma*data[('gas', 'specific_thermal_energy')] )*data[('gas', 'outflow_radial_velocity')]
   ds.add_field( ('gas', 'cell_energy_radial_outflow_flux'), function=_cell_energy_radial_outflow_flux, sampling_type=sampling_type, units='erg*km/s' )

   def _cell_mass_z_outflow_flux( field, data ):
      return data[('gas', 'mass')] * data[('gas', 'outflow_z_velocity')]
   ds.add_field( ('gas', 'cell_mass_z_outflow_flux'), function=_cell_mass_z_outflow_flux, sampling_type=sampling_type, units='Msun*km/s' )

   def _cell_energy_z_outflow_flux( field, data ):
      return data[('gas', 'mass')]*( 0.5*data[('gas', 'velocity_magnitude')]**2 + gamma*data[('gas', 'specific_thermal_energy')] )*data[('gas', 'outflow_z_velocity')]
   ds.add_field( ('gas', 'cell_energy_z_outflow_flux'), function=_cell_energy_z_outflow_flux, sampling_type=sampling_type, units='erg*km/s' )

   if ('gas', 'magnetic_energy_density') in ds.derived_field_list:
      def _cell_magnetic_energy_z_outflow_flux( field, data ):
         return data[('gas', 'volume')]*data[('gas', 'magnetic_energy_density')]*data[('gas', 'outflow_z_velocity')]
      ds.add_field( ('gas', 'cell_magnetic_energy_z_outflow_flux'), function=_cell_magnetic_energy_z_outflow_flux, sampling_type=sampling_type, units='erg*km/s' )

      def _plasma_beta_k(field, data):
         return data[('gas', 'kinetic_energy_density')] / data[('gas', 'magnetic_energy_density')]
      ds.add_field( ('gas', 'plasma_beta_k'), function=_plasma_beta_k, sampling_type=sampling_type, units='dimensionless' )

      def _plasma_beta_dual(field, data):
         return data[('gas', 'dual_internal_energy_density')] / data[('gas', 'magnetic_energy_density')]
      ds.add_field( ('gas', 'plasma_beta_dual'), function=_plasma_beta_dual, sampling_type=sampling_type, units='dimensionless' )

   if ('gamer', 'Dust') in ds.derived_field_list:
      def _cell_dust_mass_radial_outflow_flux( field, data ):
         return data[('gas', 'mass')] * data[('gas', 'outflow_radial_velocity')] * data[('gas', 'dust2gas')]
      ds.add_field( ('gas', 'cell_dust_mass_radial_outflow_flux'), function=_cell_dust_mass_radial_outflow_flux, sampling_type=sampling_type, units='Msun*km/s' )

      def _cell_dust_mass_z_outflow_flux( field, data ):
         return data[('gas', 'mass')] * data[('gas', 'outflow_z_velocity')] * data[('gas', 'dust2gas')]
      ds.add_field( ('gas', 'cell_dust_mass_z_outflow_flux'), function=_cell_dust_mass_z_outflow_flux, sampling_type=sampling_type, units='Msun*km/s' )

   if ('gamer', 'GrackleTCool') in ds.derived_field_list
      def _grackle_Tcool( field, data ):
         return data[('gamer', 'GrackleTCool')] * data.ds.quan(1, 'code_time')
      ds.add_field( ('gas', 'Tcool'), function=_grackle_Tcool, sampling_type=sampling_type, units='Myr' )


   # auxiliary fields
   def _Tsqr_over_rho( field, data ):
      return data[('gas', 'T')]**2 / data[('gas', 'density')]
   ds.add_field( ('gas', 'Tsqr_over_rho'), function=_Tsqr_over_rho, sampling_type=sampling_type, units='K**2*cm**3/g' )


   # add the particle filter
   ds.add_particle_filter( 'Halo'     )
   ds.add_particle_filter( 'Disk'     )
   ds.add_particle_filter( 'new_star' )
   if ('new_star', 'particle_mass') in ds.derived_field_list:
      ds.add_particle_filter( 'young_star' )
   if ds.dataset_type == 'gamer':
      ds.add_particle_filter( 'exp_SNII' )
      ds.add_particle_filter( 'young_SNII' )
