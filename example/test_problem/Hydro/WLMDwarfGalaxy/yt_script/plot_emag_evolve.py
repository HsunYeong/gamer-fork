import argparse
from matplotlib import patheffects
import sys
import yt
import matplotlib.pyplot as plt
import WLMDwarfGalaxy_load_datasets

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

dpi          = 150

yt.enable_parallelism()


# load the dataset
ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)
my_storage = {}

# main loop
for sto, ds in ts.piter(storage=my_storage):

   width_x =  6.0
   width_y =  6.0
   width_z =  3.0

   # decide the center
   center = ds.domain_center
   box    = ds.box( center-0.5*ds.arr([1.0*width_x, 1.0*width_y, 1.0*width_z], 'kpc'),
                    center+0.5*ds.arr([1.0*width_x, 1.0*width_y, 1.0*width_z], 'kpc') )

   emag = box.quantities.weighted_average_quantity( 'magnetic_energy_density', 'cell_volume' ).in_units('G**2').d
   time = ds.current_time.in_units('Myr').d
   sto.result = {
        "time" : time,
        "emag" : emag,
   }

# plot
if yt.is_root():
   time_all = []
   emag_all = []
   for idx, val in sorted(my_storage.items()):
      time_all.append( val["time"] )
      emag_all.append( val["emag"] )

   plt.plot( time_all, emag_all )
   plt.yscale('log')
   #plt.xlim( 0.0, 825 )
   #plt.ylim( 3.0e-5, 2.0e-2 )
   plt.xlabel( '$\mathrm{t\ [Myr]}$',                 fontsize='large' )
   plt.ylabel( '$\mathrm{magnetic\ energy\ [G^2]}$', fontsize='large' )

   # save figure
   plt.savefig( 'fig__emag.png', bbox_inches='tight', pad_inches=0.05, dpi=dpi )
