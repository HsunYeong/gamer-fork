import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Slice of mass density' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
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
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

field         = 'Phase'
colormap_dens = 'inferno'
center_mode   = 'c'
dpi           = 150

yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sz_dens = yt.SlicePlot( ds, 'z', field, center=center_mode )

   sz_dens.set_zlim( field, -1.0e0, +1.0e0 )
   sz_dens.set_cmap( field, colormap_dens )
   sz_dens.annotate_timestamp( time_unit='Gyr', corner='upper_right' )
   sz_dens.annotate_grids()

   sz_dens.save( mpl_kwargs={"dpi":dpi} )
