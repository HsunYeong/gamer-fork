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

field     = 'Dens'
colormap  = 'algae'
dpi       = 150
fontsize  = 24
titlepad  = 16

yt.enable_parallelism()
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
   dd = ds.all_data()
   mean_field = np.mean(dd[field].d)
   fmin, fmax = np.min(dd[field].d), np.max(dd[field].d)

   fig = plt.figure()
   fig.dpi = dpi
   grid = AxesGrid( fig, (0.1, 0.05, 3.2, 2.7), nrows_ncols=(1, 3), axes_pad=(1.2,0.5), label_mode="all", share_all=True, cbar_location="right", cbar_mode="single", cbar_size="2%", cbar_pad="2%")

   slc = [None]*3
   for i in range(3):
      slc[i] = yt.SlicePlot( ds, i, fields = field, center = 'c')
      slc[i].set_zlim( field, fmin, fmax, dynamic_range=None)
      slc[i].set_axes_unit( 'code_length' )
      slc[i].set_unit( field, 'code_density')
      slc[i].set_cmap( field, colormap )
      slc[i].annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
      slc[i].set_font( {'size':fontsize} )
      plot = slc[i].plots[field]
      plot.figure = fig
      plot.axes = grid[i].axes
      plot.cax = grid.cbar_axes[i]
      slc[i]._setup_plots()
      grid[i].set_title("Slice %s"%(chr(ord('x') + i)), fontsize=fontsize, pad=titlepad)

   fig.savefig("fig_%s_%s_Slice.png"%(ds, field), bbox_inches='tight',pad_inches=0.02)


