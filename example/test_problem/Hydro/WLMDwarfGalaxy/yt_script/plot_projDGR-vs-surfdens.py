import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
from mpl_toolkits.axes_grid1 import AxesGrid
import argparse
import sys
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

plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
dpi         = 150

idx_min     = 30 if code == 'GAMER' else 150
idx_sta     = max( idx_start, idx_min ) if idx_end > idx_min else idx_start
didx_avg    = 1 if code == 'GAMER' else 5
didx_avg    = max( didx_avg, didx )

# load the dataset
ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)

WLMDwarfGalaxy_derived_fields.set_particle_types(code)

nbin = 1024
dmin = 3.0e-3
dmax = 3.0e+3
ratio = (dmax/dmin)**(1.0/(nbin - 1))

dust2gas_sum_all = np.zeros(nbin)
count_all        = np.zeros(nbin)
surfdens_arr = dmin * ratio**np.arange(nbin)
surfdens_arr_boundary = (surfdens_arr*ratio**0.5)[:-1]

FONT_SIZE=24
TICK_SIZE=20

# main loop
for idx in range(idx_sta, idx_end+1, didx_avg):
   ds = yt.load( '../Data_%06d'%idx )

   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds)

   center = ds.domain_center

   p = yt.ProjectionPlot( ds, 'z', 'density',  center=center, width=(6.0, 'kpc'), method='integrate', weight_field= None,     buff_size=(1024, 1024) )
   d = yt.ProjectionPlot( ds, 'z', 'dust2gas', center=center, width=(6.0, 'kpc'), method='integrate', weight_field='density', buff_size=(1024, 1024) )


   surfdens = p.frb['density' ].in_units('Msun/pc**2').d
   dust2gas = d.frb['dust2gas'].d

   dust2gas_arr = np.zeros(nbin)

   surfdens_flat = surfdens.ravel()
   dust2gas_flat = dust2gas.ravel()

   sigma_idx = np.searchsorted(surfdens_arr_boundary, surfdens_flat)

   dust2gas_sum = np.bincount(sigma_idx, weights=dust2gas_flat, minlength=nbin)
   counts       = np.bincount(sigma_idx, minlength=nbin)

   dust2gas_sum_all += dust2gas_sum
   count_all        += counts
   mask = counts > 0
   dust2gas_arr[mask] = dust2gas_sum[mask] / counts[mask]

   plt.figure(figsize=(8.5,7.0), dpi = dpi)
   plt.plot(surfdens_arr, dust2gas_arr)
   plt.xlabel(r"$\Sigma_\mathrm{g} (\mathrm{M}_{\odot}\mathrm{pc}^{-2})$", fontsize=FONT_SIZE)
   plt.ylabel(r"$Z_\mathrm{d}^\mathrm{proj}$", fontsize=FONT_SIZE)
   plt.xticks(fontsize=TICK_SIZE)
   plt.yticks(fontsize=TICK_SIZE)
   plt.xscale('log')
   plt.yscale('log')
   plt.ylim(1e-4, 1e-2)
   plt.xlim(dmin, dmax)
   plt.savefig("imgs_DGR/fig_Data_%06d_Surface_Density_vs_Dust2Gas.png"%idx, bbox_inches="tight", pad_inches=0.05, dpi = dpi)
   plt.close()

# get average value
mask = count_all > 0
dust2gas_arr = np.zeros(nbin)
dust2gas_arr[mask] = dust2gas_sum_all[mask] / count_all[mask]

plt.figure(figsize=(8.5,7.0), dpi = dpi)
#plt.title("Gas Surface Density vs Projected Dust2Gas Ratio")
plt.plot(surfdens_arr, dust2gas_arr)
plt.xlabel(r"$\Sigma_\mathrm{g} (\mathrm{M}_{\odot}\mathrm{pc}^{-2})$", fontsize=FONT_SIZE)
plt.ylabel(r"$Z_\mathrm{d}^\mathrm{proj}$", fontsize=FONT_SIZE)
plt.xticks(fontsize=TICK_SIZE)
plt.yticks(fontsize=TICK_SIZE)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-4, 1e-2)
plt.xlim(dmin, dmax)
plt.savefig("imgs_DGR/fig_Time-Averaged_Surface_Density_vs_Dust2Gas.png", bbox_inches="tight", pad_inches=0.05, dpi = dpi)
plt.close()

