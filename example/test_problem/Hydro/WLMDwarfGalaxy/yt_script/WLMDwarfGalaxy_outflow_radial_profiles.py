import yt
import os
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTArray

markersize  = 4.0

# Analytic Mode
SFR         = 3.0e-3        # Msun/yr

E_dot       = 9.3e38        # erg/s
alpha       = E_dot/9.3e38  # dimensionless

M_dot       = 5.0e-3 *5.0   # Msun/yr
beta        = M_dot/SFR     # dimensionless

R           = 0.15   *5.0   # kpc

rho_0       = 7.1e-28 *(alpha**-0.5) *(beta**1.5)  *(R/0.15)**-2 # g/cm^3
u_inf       = 9.9e2   *(alpha**0.5)  *(beta**-0.5)               # km/s
P_0         = 6.6e-13 *(alpha**0.5)  *(beta**0.5)  *(R/0.15)**-2 # erg/cm^3
T_0         = 6.8e6   *(alpha)       *(beta**-1.)                # K
cs_0        = 3.9e2   *(alpha**0.5)  *(beta**-0.5)               # km/s
Mach_0      = 2.5                                                # dimensionless

ls          = {
                "all"      :"r-s",
                "hot"      :"m-",
                "warm-hot" :"y-",
                "warm-cool":"g-",
                "cool"     :"c-",
                "cold"     :"b-",
              }


def create_profiles(ds, outflow_region, outflow_p_type, phase):

   nbin       = 64

   min_radius = ds.quan(  2.5, 'kpc' )
   max_radius = ds.quan( 45.0, 'kpc' )

   prof_v     = yt.create_profile( outflow_region, (outflow_p_type, 'outflow_radius'), fields=[(outflow_p_type, 'density'),(outflow_p_type, 'pressure')],
                                   weight_field=(outflow_p_type, 'volume'),
                                   n_bins=nbin,
                                   units={(outflow_p_type, 'outflow_radius'): 'kpc',(outflow_p_type, 'density'): 'g/cm**3'},
                                   extrema = {(outflow_p_type, 'outflow_radius'): (min_radius, max_radius)},
                                   logs={(outflow_p_type, 'outflow_radius'): False} )

   prof_m     = yt.create_profile( outflow_region, (outflow_p_type, 'outflow_radius'), fields=[(outflow_p_type, 'T'),(outflow_p_type, 'specific_entropy'),(outflow_p_type, 'outflow_radial_velocity')],
                                   weight_field=(outflow_p_type, 'mass'),
                                   n_bins=nbin,
                                   units={(outflow_p_type, 'outflow_radius'): 'kpc',(outflow_p_type, 'T'): 'K', (outflow_p_type, 'outflow_radial_velocity'):'km/s'},
                                   extrema = {(outflow_p_type, 'outflow_radius'): (min_radius, max_radius)},
                                   logs={(outflow_p_type, 'outflow_radius'): False} )

   prof_n     = yt.create_profile( outflow_region, (outflow_p_type, 'outflow_radius'), fields=[(outflow_p_type, 'cell_mass_radial_outflow_flux'), (outflow_p_type, 'cell_energy_radial_outflow_flux')],
                                   weight_field=None,
                                   n_bins=nbin,
                                   units={(outflow_p_type, 'outflow_radius'): 'kpc', (outflow_p_type, 'cell_mass_radial_outflow_flux'): 'Msun*km/s', (outflow_p_type, 'cell_energy_radial_outflow_flux'): 'erg*km/s'},
                                   extrema = {(outflow_p_type, 'outflow_radius'): (min_radius, max_radius)},
                                   logs={(outflow_p_type, 'outflow_radius'): False} )

   gas_dens   = prof_v[(outflow_p_type, 'density')].in_units('g/cm**3').d
   gas_pres   = prof_v[(outflow_p_type, 'pressure')].in_units('erg/cm**3').d
   gas_temp   = prof_m[(outflow_p_type, 'T')].in_units('K').d
   gas_entr   = prof_m[(outflow_p_type, 'specific_entropy')].in_units('erg*cm**2').d
   gas_velr   = prof_m[(outflow_p_type, 'outflow_radial_velocity')].in_units('km/s').d
   radius_v   = prof_v.x.in_units('kpc').d
   radius_m   = prof_m.x.in_units('kpc').d

   dr         = prof_n.x[1] - prof_n.x[0]  # assuming linear bin
   gas_mout   = ( prof_n[(outflow_p_type, 'cell_mass_radial_outflow_flux')] / dr ).in_units('Msun/yr').d
   gas_eout   = ( prof_n[(outflow_p_type, 'cell_energy_radial_outflow_flux')] / dr ).in_units('erg/s').d
   radius_n   = prof_n.x.in_units('kpc').d

   gas_cs     = (YTArray( (5./3.)*np.divide(gas_pres, gas_dens, out=np.zeros_like(gas_dens), where=gas_dens>0), "erg/g" )**0.5).in_units('km/s').d
   gas_mach   = np.divide(gas_velr, gas_cs, out=np.zeros_like(gas_cs), where=gas_cs>0)


   # save to file
   np.savetxt( './profs/%s_GalacticOutflow_Profiles_%s'%(ds, phase),
               np.column_stack( (radius_v, gas_dens, gas_pres, radius_m, gas_temp, gas_entr, gas_velr, gas_cs, gas_mach, radius_n, gas_mout, gas_eout) ),
               fmt='%19.8e',
               header='%17s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s'%('r_v', 'dens', 'pres', 'r_m', 'temp', 'entr', 'velr', 'cs', 'mach', 'r_n', 'mout', 'eout') )

   return radius_v, gas_dens, gas_pres, radius_m, gas_temp, gas_entr, gas_velr, gas_cs, gas_mach, radius_n, gas_mout, gas_eout


def plot_outflow_radial_profiles_initial():
   f_orp, ax_orp = plt.subplots( 2, 4, figsize=(14.0, 4.8) )
   f_orp.subplots_adjust( wspace=0.4 )
   return f_orp, ax_orp

def plot_outflow_radial_profiles(ax_orp, code, indices, phase):

   gas_dens    = 0
   gas_pres    = 0
   gas_temp    = 0
   gas_entr    = 0
   gas_velr    = 0
   gas_cs      = 0
   gas_mach    = 0
   gas_mout    = 0
   gas_eout    = 0
   N_avg_skip0 = 0
   N_avg_total = 0

   for idx in indices[::-1]:
      if code == 'GAMER':
         filename = './profs/Data_%06d_GalacticOutflow_Profiles_%s'%(idx, phase)
      elif code == 'GIZMO':
         filename = './profs/snap_%03d_GalacticOutflow_Profiles_%s'%(idx, phase)

      if not os.path.exists(filename):   continue

      radius_v_i, gas_dens_i, gas_pres_i, radius_m_i, gas_temp_i, gas_entr_i, gas_velr_i, gas_cs_i, gas_mach_i, radius_n_i, gas_mout_i, gas_eout_i = np.loadtxt( filename, skiprows=1, unpack=True )

      if idx == indices[-1]:
         radius_v = radius_v_i
         radius_m = radius_m_i
         radius_n = radius_n_i

      mask = (gas_dens_i > 0)
      gas_dens    += np.where(mask, gas_dens_i, 0)
      gas_pres    += np.where(mask, gas_pres_i, 0)
      gas_temp    += np.where(mask, gas_temp_i, 0)
      gas_entr    += np.where(mask, gas_entr_i, 0)
      gas_velr    += np.where(mask, gas_velr_i, 0)
      gas_cs      += np.where(mask, gas_cs_i  , 0)
      gas_mach    += np.where(mask, gas_mach_i, 0)
      gas_mout    += np.where(mask, gas_mout_i, 0)
      gas_eout    += np.where(mask, gas_eout_i, 0)
      N_avg_skip0 += np.where(mask, np.ones_like(gas_dens_i), 0)
      N_avg_total += np.ones_like(gas_dens_i)

   mask = N_avg_skip0 > 0
   gas_dens       /= N_avg_total
   gas_pres       /= N_avg_total
   gas_temp[mask] /= N_avg_skip0[mask]
   gas_entr[mask] /= N_avg_skip0[mask]
   gas_velr[mask] /= N_avg_skip0[mask]
   gas_cs  [mask] /= N_avg_skip0[mask]
   gas_mach[mask] /= N_avg_skip0[mask]
   gas_mout       /= N_avg_total
   gas_eout       /= N_avg_total

   ax_orp[0][0].plot( radius_v,  gas_dens,                 ls[phase],  lw=1, mec='none', ms=markersize, label=phase )
   ax_orp[0][1].plot( radius_v,  gas_pres,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_orp[0][2].plot( radius_m,  gas_temp,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_orp[0][3].plot( radius_m,  gas_velr,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_orp[1][0].plot( radius_n,  gas_mout,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_orp[1][1].plot( radius_n,  gas_eout,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_orp[1][2].plot( radius_m,  gas_entr,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_orp[1][3].plot( radius_m,  gas_mach,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )

   if phase != 'all':
      return

   # density
   ax_orp[0][0].plot( radius_v,  rho_0*(radius_v/R)**(-2.),    'k--',  lw=1, mec='none', ms=markersize, label=r'$r^{-2}$'            )

   # pressure
   ax_orp[0][1].plot( radius_v,  P_0*(radius_v/R)**(-10./3.),  'k--',  lw=1, mec='none', ms=markersize, label=r'$r^{-10/3}$'         )

   # temperature
   ax_orp[0][2].plot( radius_m,  T_0*(radius_m/R)**(-4./3.),   'k--',  lw=1, mec='none', ms=markersize, label=r'$r^{-4/3}$'          )

   # velocity
   ax_orp[0][3].plot( radius_m,  u_inf*(radius_m/R)**(0.0),    'k--',  lw=1, mec='none', ms=markersize, label=r'$r^{0}$'             )

   # mass outflow rate
   ax_orp[1][0].plot( radius_n,  M_dot*(radius_n/R)**(0.0),    'k--',  lw=1, mec='none', ms=markersize, label=r'$\dot{M}_{\rm inj}$' )
   ax_orp[1][0].plot( radius_n,  SFR*(radius_n/R)**(0.0),      'k-.',  lw=1, mec='none', ms=markersize, label=r'SFR'                 )
   ax_orp[1][0].axvline( x=R, c='k', ls=':', lw=1 )

   # energy outflow rate
   ax_orp[1][1].plot( radius_n,  E_dot*(radius_n/R)**(0.0),    'k--',  lw=1, mec='none', ms=markersize, label=r'$\dot{E}_{\rm inj}$' )
   ax_orp[1][1].axvline( x=R, c='k', ls=':', lw=1 )

   # entropy
   mumH = 0.588235294117647*1.007825*1.660539040e-24
   K_0  = P_0*(rho_0/mumH)**(-5./3.)
   ax_orp[1][2].plot( radius_v,  K_0*(radius_v/R)**(0.0),      'k--',  lw=1, mec='none', ms=markersize, label=r'$r^{0}$'             )

   #   Mach number
   ax_orp[1][3].plot( radius_m,  Mach_0*(radius_m/R)**(2./3.), 'k--',  lw=1, mec='none', ms=markersize, label=r'$r^{2/3}$'           )


def plot_outflow_radial_profiles_final(f_orp, ax_orp, code, indices, title, fileout):
   ax_orp[0][0].set_ylim( 5.0e-32, 5.0e-28 )
   ax_orp[0][0].set_ylabel( r'$\rho_\mathrm{gas}\ \mathrm{[g\ cm^{-3}]}$',          fontsize='large' )
   ax_orp[0][1].set_ylim( 3.0e-20, 3.0e-15 )
   ax_orp[0][1].set_ylabel( r'$P_\mathrm{gas}\ \mathrm{[erg\ cm^{-3}]}$',           fontsize='large' )
   ax_orp[0][2].set_ylim( 2.0e+3, 4.0e+5 )
   ax_orp[0][2].set_ylabel( r'$T_\mathrm{gas}\ \mathrm{[K]}$',                      fontsize='large' )
   ax_orp[0][3].set_ylim( 8.0e+0, 4.0e+2 )
   ax_orp[0][3].set_ylabel( r'$v_\mathrm{r}\ \mathrm{[km/s]}$',                     fontsize='large' )
   ax_orp[1][0].set_ylim( 1.0e-6, 3.0e-2 )
   ax_orp[1][0].set_ylabel( r'$\dot{M}_\mathrm{out}\ \mathrm{[M_\odot\ yr^{-1}]}$', fontsize='large' )
   ax_orp[1][1].set_ylim( 2.0e34, 2.0e39 )
   ax_orp[1][1].set_ylabel( r'$\dot{E}_\mathrm{out}\ \mathrm{[erg\ s^{-1}]}$',      fontsize='large' )
   ax_orp[1][2].set_ylim( 1.0e-10, 2.0e-6 )
   ax_orp[1][2].set_ylabel( r'$K_\mathrm{gas}\ \mathrm{[erg\ cm^2]}$',              fontsize='large' )
   ax_orp[1][3].set_ylim( 4.0e-1, 4.0e+1 )
   ax_orp[1][3].set_ylabel( r'$\mathcal{M}$',                                       fontsize='large' )
   for i in range(2):
      for j in range(4):
         ax_orp[i][j].legend(loc=1+3*i, labelspacing=0.3)
         ax_orp[i][j].set_xlim( 0.0, 45.0 )
         ax_orp[i][j].set_xscale( 'linear' )
         ax_orp[i][j].set_yscale( 'log', nonpositive='clip' )
         ax_orp[1][j].set_xlabel( '$r\ \mathrm{[kpc]}$', fontsize='large' )

   f_orp.suptitle( title, fontsize='large' )
   f_orp.savefig( './imgs_o/fig_'+fileout+'_outflow_profiles.png', bbox_inches='tight', pad_inches=0.05, dpi=150 )
   plt.close(f_orp)
