import yt
import os
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTArray

markersize  = 4.0

# Analytic Mode
SFR         = 3.0e-3        # Msun/yr

E_dot       = 9.3e38        # erg/s

Z           = 0.15   *5.0   # kpc

ls          = {
                "all"      :"r-s",
                "hot"      :"m-",
                "warm-hot" :"y-",
                "warm-cool":"g-",
                "cool"     :"c-",
                "cold"     :"b-",
              }


def create_profiles(ds, outflow_region, outflow_p_type, phase):

   nbin = 32

   min_z = ds.quan(  2.5, 'kpc' )
   max_z = ds.quan( 45.0, 'kpc' )

   prof_v     = yt.create_profile( outflow_region, (outflow_p_type, 'outflow_z'), fields=[(outflow_p_type, 'density'),(outflow_p_type, 'pressure')],
                                   weight_field=(outflow_p_type, 'volume'),
                                   n_bins=nbin,
                                   units={(outflow_p_type, 'outflow_z'): 'kpc',(outflow_p_type, 'density'): 'g/cm**3'},
                                   extrema = {(outflow_p_type, 'outflow_z'): (min_z, max_z)},
                                   logs={(outflow_p_type, 'outflow_z'): False} )

   prof_m     = yt.create_profile( outflow_region, (outflow_p_type, 'outflow_z'), fields=[(outflow_p_type, 'T'),(outflow_p_type, 'specific_entropy'),(outflow_p_type, 'outflow_z_velocity')],
                                   weight_field=(outflow_p_type, 'mass'),
                                   n_bins=nbin,
                                   units={(outflow_p_type, 'outflow_z'): 'kpc',(outflow_p_type, 'T'): 'K', (outflow_p_type, 'outflow_z_velocity'):'km/s'},
                                   extrema = {(outflow_p_type, 'outflow_z'): (min_z, max_z)},
                                   logs={(outflow_p_type, 'outflow_z'): False} )

   prof_n     = yt.create_profile( outflow_region, (outflow_p_type, 'outflow_z'), fields=[(outflow_p_type, 'cell_mass_z_outflow_flux'), (outflow_p_type, 'cell_energy_z_outflow_flux')],
                                   weight_field=None,
                                   n_bins=nbin,
                                   units={(outflow_p_type, 'outflow_z'): 'kpc', (outflow_p_type, 'cell_mass_z_outflow_flux'): 'Msun*km/s', (outflow_p_type, 'cell_energy_z_outflow_flux'): 'erg*km/s'},
                                   extrema = {(outflow_p_type, 'outflow_z'): (min_z, max_z)},
                                   logs={(outflow_p_type, 'outflow_z'): False} )

   gas_dens   = prof_v[(outflow_p_type, 'density')].in_units('g/cm**3').d
   gas_pres   = prof_v[(outflow_p_type, 'pressure')].in_units('erg/cm**3').d
   gas_temp   = prof_m[(outflow_p_type, 'T')].in_units('K').d
   gas_entr   = prof_m[(outflow_p_type, 'specific_entropy')].in_units('erg*cm**2').d
   gas_velr   = prof_m[(outflow_p_type, 'outflow_z_velocity')].in_units('km/s').d
   z_v        = prof_v.x.in_units('kpc').d
   z_m        = prof_m.x.in_units('kpc').d

   dz         = prof_n.x[1] - prof_n.x[0]  # assuming linear bin
   gas_mout   = ( prof_n[(outflow_p_type, 'cell_mass_z_outflow_flux')] / dz ).in_units('Msun/yr').d
   gas_eout   = ( prof_n[(outflow_p_type, 'cell_energy_z_outflow_flux')] / dz ).in_units('erg/s').d
   z_n        = prof_n.x.in_units('kpc').d

   gas_cs     = (YTArray( (5./3.)*np.divide(gas_pres, gas_dens, out=np.zeros(len(gas_dens)), where=gas_dens>0), "erg/g" )**0.5).in_units('km/s').d
   gas_mach   = np.divide(gas_velr, gas_cs, out=np.zeros(len(gas_cs)), where=gas_cs>0)


   # save to file
   np.savetxt( './profs/%s_GalacticOutflow_Z_Profiles_%s'%(ds, phase),
               np.column_stack( (z_v, gas_dens, gas_pres, z_m, gas_temp, gas_entr, gas_velr, gas_cs, gas_mach, z_n, gas_mout, gas_eout) ),
               fmt='%19.8e',
               header='%17s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s'%('z_v', 'dens', 'pres', 'z_m', 'temp', 'entr', 'velr', 'cs', 'mach', 'z_n', 'mout', 'eout') )

   return z_v, gas_dens, gas_pres, z_m, gas_temp, gas_entr, gas_velr, gas_cs, gas_mach, z_n, gas_mout, gas_eout


def plot_outflow_z_profiles_initial():
   f_ozp, ax_ozp = plt.subplots( 2, 4, figsize=(14.0, 4.8) )
   f_ozp.subplots_adjust( wspace=0.4 )
   return f_ozp, ax_ozp

def plot_outflow_z_profiles(ax_ozp, code, indices, phase):

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
         filename = './profs/Data_%06d_GalacticOutflow_Z_Profiles_%s'%(idx, phase)
      elif code == 'GIZMO':
         filename = './profs/snap_%03d_GalacticOutflow_Z_Profiles_%s'%(idx, phase)

      if not os.path.exists(filename):   continue

      z_v_i, gas_dens_i, gas_pres_i, z_m_i, gas_temp_i, gas_entr_i, gas_velr_i, gas_cs_i, gas_mach_i, z_n_i, gas_mout_i, gas_eout_i = np.loadtxt( filename, skiprows=1, unpack=True )

      if idx == indices[-1]:
         z_v = z_v_i
         z_m = z_m_i
         z_n = z_n_i

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

   ax_ozp[0][0].plot( z_v,  gas_dens,                 ls[phase],  lw=1, mec='none', ms=markersize, label=phase )
   ax_ozp[0][1].plot( z_v,  gas_pres,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_ozp[0][2].plot( z_m,  gas_temp,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_ozp[0][3].plot( z_m,  gas_velr,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_ozp[1][0].plot( z_n,  gas_mout,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_ozp[1][1].plot( z_n,  gas_eout,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_ozp[1][2].plot( z_m,  gas_entr,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )
   ax_ozp[1][3].plot( z_m,  gas_mach,                 ls[phase],  lw=1, mec='none', ms=markersize, label=''    )

   if phase != 'all':
      return

   # mass outflow rate
   ax_ozp[1][0].plot( z_n,  SFR*(z_n/Z)**(0.0),      'k-.',  lw=1, mec='none', ms=markersize, label=r'SFR'                 )

   # energy outflow rate
   ax_ozp[1][1].plot( z_n,  E_dot*(z_n/Z)**(0.0),    'k--',  lw=1, mec='none', ms=markersize, label=r'$\dot{E}_{\rm inj}$' )


def plot_outflow_z_profiles_final(f_ozp, ax_ozp, code, indices, title, fileout):
   ax_ozp[0][0].set_ylim( 5.0e-32, 5.0e-28 )
   ax_ozp[0][0].set_ylabel( r'$\rho_\mathrm{gas}\ \mathrm{[g\ cm^{-3}]}$',          fontsize='large' )
   ax_ozp[0][1].set_ylim( 3.0e-20, 3.0e-15 )
   ax_ozp[0][1].set_ylabel( r'$P_\mathrm{gas}\ \mathrm{[erg\ cm^{-3}]}$',           fontsize='large' )
   ax_ozp[0][2].set_ylim( 2.0e+3, 4.0e+5 )
   ax_ozp[0][2].set_ylabel( r'$T_\mathrm{gas}\ \mathrm{[K]}$',                      fontsize='large' )
   ax_ozp[0][3].set_ylim( 8.0e+0, 4.0e+2 )
   ax_ozp[0][3].set_ylabel( r'$v_\mathrm{z}\ \mathrm{[km/s]}$',                     fontsize='large' )
   ax_ozp[1][0].set_ylim( 1.0e-6, 3.0e-2 )
   ax_ozp[1][0].set_ylabel( r'$\dot{M}_\mathrm{out}\ \mathrm{[M_\odot\ yr^{-1}]}$', fontsize='large' )
   ax_ozp[1][1].set_ylim( 2.0e34, 2.0e39 )
   ax_ozp[1][1].set_ylabel( r'$\dot{E}_\mathrm{out}\ \mathrm{[erg\ s^{-1}]}$',      fontsize='large' )
   ax_ozp[1][2].set_ylim( 1.0e-10, 2.0e-6 )
   ax_ozp[1][2].set_ylabel( r'$K_\mathrm{gas}\ \mathrm{[erg\ cm^2]}$',              fontsize='large' )
   ax_ozp[1][3].set_ylim( 4.0e-1, 4.0e+1 )
   ax_ozp[1][3].set_ylabel( r'$\mathcal{M}$',                                       fontsize='large' )
   for i in range(2):
      for j in range(4):
         ax_ozp[i][0].legend(loc=1+3*i)
         ax_ozp[i][j].set_xlim( 0.0, 45.0 )
         ax_ozp[i][j].set_xscale( 'linear' )
         ax_ozp[i][j].set_yscale( 'log', nonpositive='clip' )
         ax_ozp[1][j].set_xlabel( '$z\ \mathrm{[kpc]}$', fontsize='large' )

   f_ozp.suptitle( title, fontsize='large' )
   f_ozp.savefig( './imgs_o/fig_'+fileout+'_outflow_z_profiles.png', bbox_inches='tight', pad_inches=0.05, dpi=150 )
   plt.close(f_ozp)
