import argparse
import sys
import yt
import numpy as np
import WLMDwarfGalaxy_load_datasets
import WLMDwarfGalaxy_derived_fields

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the star mass distribution' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = '../'
code        = 'GAMER'

colormap    = 'algae'
dpi         = 150


yt.enable_parallelism()

ts = WLMDwarfGalaxy_load_datasets.load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx)

WLMDwarfGalaxy_derived_fields.set_particle_types(code)

x_lim_min = 0.0
x_lim_max = 8.25e2
y_lim_min = 0.0
y_lim_max = 3.0e1

for ds in ts.piter():

   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds)

#  check particle
   if len( ds.all_data()['new_star', 'particle_ones'] ) == 0:
      print( 'WARNING: There is no particle in %s !!'%ds )
      continue

#  plot
   p = yt.PhasePlot( ds, ('new_star', 'ParCreTime'), ('new_star', 'particle_mass'), ('new_star', 'particle_ones'),
                     weight_field=None, x_bins=256, y_bins=256 )

   p.set_unit( ('new_star', 'ParCreTime'),    'Myr'  )
   p.set_unit( ('new_star', 'particle_mass'), 'Msun' )
   p.set_log(  ('new_star', 'ParCreTime'),     False )
   p.set_log(  ('new_star', 'particle_mass'),  False )
   p.set_xlim( x_lim_min, x_lim_max )
   p.set_ylim( y_lim_min, y_lim_max )
   p.set_cmap( ('new_star', 'particle_ones'), colormap )
   p.save( 'fig_%s_star_mass_distribution.png'%(ds), mpl_kwargs={'dpi':dpi} )
