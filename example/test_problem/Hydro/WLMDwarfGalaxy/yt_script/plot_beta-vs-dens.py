import argparse
import sys
import yt
import WLMDwarfGalaxy_BetaDens_Phase_and_PDF
import WLMDwarfGalaxy_load_datasets
import WLMDwarfGalaxy_derived_fields
from yt.utilities.parallel_tools.parallel_analysis_interface import communication_system

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas density-temperature phase diagram' )

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


width_kpc   = 9
nbin        = 300

x_lim_min   = 9.0e-32
x_lim_max   = 2.0e-18
y_lim_min   = 1.0e0
y_lim_max   = 1.0e12

hasDust = True

yt.enable_parallelism()

ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)

for ds in ts.piter():

   idx = int(str(ds)[5:11]) if code == 'GAMER' else int(str(ds)[5:11])
   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds, hasDust)

#  only include the data within a sphere with a radius of width_kpc
   sp = ds.sphere( ds.domain_center, (0.5*width_kpc, 'kpc') )

   WLMDwarfGalaxy_BetaDens_Phase_and_PDF.create_PhaseDiagram( ds, sp, '', 'gas', 'mass', nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )

   WLMDwarfGalaxy_BetaDens_Phase_and_PDF.plot_PhaseDiagram( range(idx, idx+1, didx), code, '', 'mass', nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max, '$t$ = {:.1f} {:s}'.format( ds.current_time.in_units('Myr').d, 'Myr' ), ds )

comm = communication_system.communicators[-1]
comm.barrier()

if yt.is_root():
   idx_min     = 30 if code == 'GAMER' else 150
   idx_sta     = max( idx_start, idx_min ) if idx_end > idx_min else idx_start
   didx_avg    = 1 if code == 'GAMER' else 5
   didx_avg    = max( didx_avg, didx )
   indices_avg = range(idx_sta, idx_end+1, didx_avg)

   WLMDwarfGalaxy_BetaDens_Phase_and_PDF.plot_PhaseDiagram( indices_avg, code, '', 'mass', nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max, 'Time-Averaged', 'Time-Averaged' )


