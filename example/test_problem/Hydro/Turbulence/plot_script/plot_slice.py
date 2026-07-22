import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt
from mpl_toolkits.axes_grid1 import AxesGrid
import argparse
import sys

#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Get density slices' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )

'''
field       = ['velocity_magnitude', 'magnetic_field_magnitude', 'helicity'                 ]
field_unit  = ['code_velocity',      'code_magnetic',            'code_velocity/code_time'  ]
colormap    = ['magma',              'cividis',                  'RdBu'                     ]
'''
field       = ['velocity_magnitude', 'vorticity_magnitude', 'helicity'                 ]
field_unit  = ['code_velocity',      '1/code_time',         'code_velocity/code_time'  ]
colormap    = ['magma',              'cividis',             'RdBu'                     ]
field2      = ['velocity', 'vorticity']

dpi       = 150
fontsize  = 24
titlepad  = 16

yt.enable_parallelism()
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   def _helicity( field, data ):
      return data["velocity_x"]*data["vorticity_x"] + data["velocity_y"]*data["vorticity_y"] + data["velocity_z"]*data["vorticity_z"]
   ds.add_field( ("gas", "helicity"), function=_helicity, sampling_type="cell", units="code_velocity/code_time" )

   dd = ds.all_data()

   fig = plt.figure()
   fig.dpi = dpi
   grid = AxesGrid( fig, (0.1, 0.05, 3.2, 2.7), nrows_ncols=(1, 3), axes_pad=(3,0.5), label_mode="all", share_all=True, cbar_location="right", cbar_mode="each", cbar_size="2%", cbar_pad="2%")

   slc = [None]*3
   for i in range(3):

      slc[i] = yt.SlicePlot( ds, 0, fields = field[i], center = 'c')
      slc[i].set_axes_unit( 'code_length' )
      slc[i].set_unit( field[i], field_unit[i])
      slc[i].set_cmap( field[i], colormap[i] )
      slc[i].set_font( {'size':fontsize} )
      if i != 2:
#         slc[i].annotate_streamlines(("gas", field2[i]+"_x"), ("gas", field2[i]+"_y"), color='black') # linewidth=("gas", "magnetic_field_strength")
         slc[i].annotate_quiver(("gas", field2[i]+"_x"), ("gas", field2[i]+"_y"), color='black') # linewidth=("gas", "magnetic_field_strength")

      plot = slc[i].plots[field[i]]
      plot.figure = fig
      plot.axes = grid[i].axes
      plot.cax = grid.cbar_axes[i]
      slc[i]._setup_plots()
      grid[i].set_title(field[i], fontsize=fontsize, pad=titlepad)

   fig.savefig("fig_%s_Slice.png"%(ds), bbox_inches='tight',pad_inches=0.02)


