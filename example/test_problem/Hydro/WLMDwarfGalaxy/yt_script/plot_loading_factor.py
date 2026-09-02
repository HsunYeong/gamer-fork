import yt
import argparse
import sys
import os
import gc
import WLMDwarfGalaxy_load_datasets
import WLMDwarfGalaxy_derived_fields
import WLMDwarfGalaxy_galactic_inoutflow
import WLMDwarfGalaxy_outflow_radial_profiles
import WLMDwarfGalaxy_outflow_z_profiles
import WLMDwarfGalaxy_outflow_rates
import WLMDwarfGalaxy_TempDens_Phase_and_PDF
import WLMDwarfGalaxy_outflow_TV_relation
import WLMDwarfGalaxy_outflow_thinlayer
import WLMDwarfGalaxy_outflow_thinlayer_dust2gas
from yt.utilities.parallel_tools.parallel_analysis_interface import communication_system

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the galactic outflow rate' )

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

nbin        = 300
x_lim_min   = 1.0e-33
x_lim_max   = 2.0e-26
y_lim_min   = 5.0e1
y_lim_max   = 2.0e7
v_lim_min   = 8.0e-1
v_lim_max   = 3.0e+3

outflow_z_kpc  =  1.0
outflow_dz_kpc =  0.1*outflow_z_kpc

list_phases = ["hot", "warm-hot", "warm-cool", "cool", "all"]

yt.enable_parallelism()

# load data
ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)

WLMDwarfGalaxy_derived_fields.set_particle_types(code)

hasDust = True

for ds in ts.piter():

   idx = int(str(ds)[5:11]) if code == 'GAMER' else int(str(ds)[5:11])
   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds)

   for phase in list_phases:

      # outflow rates at a speific z
      WLMDwarfGalaxy_outflow_rates.calculate_galactic_outflow_rate( ds, outflow_z_kpc, outflow_dz_kpc, phase, hasDust )


comm = communication_system.communicators[-1]
comm.barrier()


if yt.is_root():
   idx_min     = 30 if code == 'GAMER' else 150
   idx_sta     = max( idx_start, idx_min ) if idx_end > idx_min else idx_start
   didx_avg    = 1 if code == 'GAMER' else 5
   didx_avg    = max( didx_avg, didx )
   indices_avg = range(idx_sta, idx_end+1, didx_avg)

   for phase in list_phases:

      for idx in indices_avg:
         prefix = 'Data_%06d'%idx if code == 'GAMER' else 'snap_%03d'%idx
         if idx == idx_start:
            os.system("head -1 %s >> %s"%('./tables/%s_Galactic_Outflow_Rate_%s_z_%02d_kpc'%(prefix, phase, int(outflow_z_kpc)), './tables/Galactic_Outflow_Rate_%s_z_%02d_kpc'%(phase, int(outflow_z_kpc))))
         os.system("tail -1 %s >> %s"%('./tables/%s_Galactic_Outflow_Rate_%s_z_%02d_kpc'%(prefix, phase, int(outflow_z_kpc)), './tables/Galactic_Outflow_Rate_%s_z_%02d_kpc'%(phase, int(outflow_z_kpc))))

   # outflow rates time-evolution
   WLMDwarfGalaxy_outflow_rates.plot_galactic_outflow_rate_evolution(outflow_z_kpc, list_phases, hasDust)


