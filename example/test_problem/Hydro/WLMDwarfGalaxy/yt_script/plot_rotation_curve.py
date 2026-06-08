import matplotlib
matplotlib.use('Agg')
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import yt
import gc

# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
figure_dpi  = 140


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='An example script to plot disk/halo data from existing .npy files' )

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

def SearchIndex(x, A, N):
   i = 0
   j = N - 1
   while(i <= j):
      mid = int(i + (j - i)/2)
      if(A[mid] == x):
         i = mid
         break
      elif(A[mid] > x):
         j = mid - 1
      else: i = mid + 1
   return i


nbin = 500

yt.enable_parallelism()

# load the dataset
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
my_storage = {}

#-------------------------------------------------------------------------------------------------------------------------
# angle-averaged rotation curve
for sto, ds in ts.piter(storage=my_storage):

   idx = ds.parameters["DumpID"]

   # define domain of interest
   width   =  20.0
   width_z =  10.0
   center = ds.domain_center
   box    = ds.box( center-0.5*ds.arr([1.0*width, 1.0*width, 1.0*width_z], 'kpc'),
                    center+0.5*ds.arr([1.0*width, 1.0*width, 1.0*width_z], 'kpc') )

   CoM = box.quantities.center_of_mass()
   VCM = [ box.quantities.weighted_average_quantity(  ('gas', 'velocity_x'), ('gas', 'cell_mass') ),
           box.quantities.weighted_average_quantity(  ('gas', 'velocity_y'), ('gas', 'cell_mass') ) ]

   mass = box[('gas',  'cell_mass')]
   posx = box[('gas',          'x')] - CoM[0]
   posy = box[('gas',          'y')] - CoM[1]
   velx = box[('gas', 'velocity_x')] - VCM[0]
   vely = box[('gas', 'velocity_y')] - VCM[1]

   posr = (posx**2 + posy**2)**0.5
   velp = (posx*vely - posy*velx)/posr
   sortR  = np.sort(posr)
   indexR = np.argsort(posr)

   mean_r  = np.zeros(nbin)
   mean_vp = np.zeros(nbin)
   ndata   = np.size(posr)

   r   = ds.quan(0.0, 'kpc')
   dr  = ds.quan(width/2.0, 'kpc')/nbin
   num = 0
   for j in range(nbin):
      num_pre = num
      num     = SearchIndex( r+dr, sortR.in_units('kpc'), ndata )

      mean_r [j] = np.average( posr[ indexR[ num_pre:num ] ].in_units('kpc' ).d, weights=mass[ indexR[ num_pre:num ] ] )
      mean_vp[j] = np.average( velp[ indexR[ num_pre:num ] ].in_units('km/s').d, weights=mass[ indexR[ num_pre:num ] ] )

      r = r + dr

   sto.result = {
        "idx": idx,
        "mean_r" : mean_r,
        "mean_vp": mean_vp,
   }
   box.clear_data()
   del posx, posy, posr, velx, vely, velp
   gc.collect()

if yt.is_root():
   results = list(my_storage.values())
   plt.figure(dpi = figure_dpi)

   for r in results:
      idx = r["idx"]
      plt.plot(r["mean_r"], r["mean_vp"], color = cm.Blues(0.3+0.6*(idx-idx_start)/(idx_end+1-idx_start)), label ='DumpID=%d'%idx)

   plt.xlim((0, 10))
   plt.ylim((0, 80))
   plt.grid(ls='--')
   plt.legend(loc='lower right', shadow=True, prop={'size':8})   ## loc='best', 'upper left', 'upper right', 'lower left', 'lower right'
   plt.xlabel(r"$R$ (kpc)")
   plt.ylabel(r"$v_{\rm cir}$ (km/s)")
   plt.savefig("fig__rotation_curve.png", dpi = figure_dpi)
   plt.close()


