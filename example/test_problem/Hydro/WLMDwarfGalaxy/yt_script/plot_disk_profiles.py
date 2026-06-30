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
parser = argparse.ArgumentParser( description='An example script to plot disk data' )

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

   idx  = ds.parameters["DumpID"]
   time = ds.quan( ds.parameters["Time"][0], 'code_time' ).in_units('Myr').d
   # define domain of interest
   radius = 8.0
   height = 3.0

   disk = ds.disk( center='c', normal=[0,0,1], radius=(radius, 'kpc'), height=(height, 'kpc') ).cut_region( ["obj['gas', 'density'].in_units('g/cm**3') > 1.0e-30"] )

   CoM = disk.quantities.center_of_mass()
   VCM = [ disk.quantities.weighted_average_quantity(  ('gas', 'velocity_x'), ('gas', 'cell_mass') ),
           disk.quantities.weighted_average_quantity(  ('gas', 'velocity_y'), ('gas', 'cell_mass') ),
           disk.quantities.weighted_average_quantity(  ('gas', 'velocity_z'), ('gas', 'cell_mass') ) ]

   mass = disk[('gas', 'cell_mass'  )]
   vol  = disk[('gas', 'cell_volume')]
   posx = disk[('gas',          'x')] - CoM[0]
   posy = disk[('gas',          'y')] - CoM[1]
   velx = disk[('gas', 'velocity_x')] - VCM[0]
   vely = disk[('gas', 'velocity_y')] - VCM[1]
   velz = disk[('gas', 'velocity_z')] - VCM[2]
   eint = disk[('gas', 'thermal_energy_density' )]
   emag = disk[('gas', 'magnetic_energy_density')]

   posr = (posx**2 + posy**2)**0.5
   velr = (posx*velx + posy*vely)/posr
   velp = (posx*vely - posy*velx)/posr
   sortR  = np.sort(posr)
   indexR = np.argsort(posr)

   mean_r   = np.zeros(nbin)
   mean_vr  = np.zeros(nbin)
   mean_vp  = np.zeros(nbin)
   mean_vz  = np.zeros(nbin)
   mean_vr2 = np.zeros(nbin)
   mean_vp2 = np.zeros(nbin)
   mean_vz2 = np.zeros(nbin)
   mean_B   = np.zeros(nbin)
   sigma_vr = np.zeros(nbin)
   sigma_vp = np.zeros(nbin)
   sigma_vz = np.zeros(nbin)
   plasma_b = np.zeros(nbin)
   surdens  = np.zeros(nbin)

   ndata   = np.size(posr)

   r   = ds.quan(0.0, 'kpc')
   dr  = ds.quan(radius, 'kpc')/nbin
   num = 0
   for j in range(nbin):
      num_pre = num
      num     = SearchIndex( r+dr, sortR.in_units('kpc'), ndata )

      mean_r [j]  = np.average( posr[ indexR[ num_pre:num ] ].in_units('kpc' ).d,    weights=mass[ indexR[ num_pre:num ] ] )
      mean_vr[j]  = np.average( velr[ indexR[ num_pre:num ] ].in_units('km/s').d,    weights=mass[ indexR[ num_pre:num ] ] )
      mean_vp[j]  = np.average( velp[ indexR[ num_pre:num ] ].in_units('km/s').d,    weights=mass[ indexR[ num_pre:num ] ] )
      mean_vz[j]  = np.average( velz[ indexR[ num_pre:num ] ].in_units('km/s').d,    weights=mass[ indexR[ num_pre:num ] ] )

      mean_vr2[j] = np.average( velr[ indexR[ num_pre:num ] ].in_units('km/s').d**2, weights=mass[ indexR[ num_pre:num ] ] )
      mean_vp2[j] = np.average( velp[ indexR[ num_pre:num ] ].in_units('km/s').d**2, weights=mass[ indexR[ num_pre:num ] ] )
      mean_vz2[j] = np.average( velz[ indexR[ num_pre:num ] ].in_units('km/s').d**2, weights=mass[ indexR[ num_pre:num ] ] )

      sigma_vr[j] = (mean_vr2[j] - mean_vr[j]**2)**0.5
      sigma_vp[j] = (mean_vp2[j] - mean_vp[j]**2)**0.5
      sigma_vz[j] = (mean_vz2[j] - mean_vz[j]**2)**0.5

      enc_mass    = np.sum( mass[ indexR[ num_pre:num ] ].in_units('Msun').d )
      enc_eint    = np.sum( (eint[ indexR[ num_pre:num ] ]*vol[ indexR[ num_pre:num ] ]) )
      enc_emag    = np.sum( (emag[ indexR[ num_pre:num ] ]*vol[ indexR[ num_pre:num ] ]) )

      area        = (np.pi *((r + dr)**2 - r**2)).in_units('kpc**2').d
      surdens[j]  = enc_mass/area
      plasma_b[j] = (enc_eint/enc_emag).d
      mean_B[j]   = np.average( emag[ indexR[ num_pre:num ] ].in_units('G**2').d, weights=vol[ indexR[ num_pre:num ] ] )
      mean_B[j]   = (mean_B[j]*8*np.pi)**0.5
      r = r + dr

   sto.result = {
        "idx"     : idx,
        "time"    : time,
        "mean_r"  : mean_r,
        "mean_vp" : mean_vp,
        "sigma_vr": sigma_vr,
        "sigma_vp": sigma_vp,
        "sigma_vz": sigma_vz,
        "surdens" : surdens,
        "plasma_b": plasma_b,
        "mean_B"  : mean_B
   }
   disk.clear_data()
   del mass, vol, posx, posy, posr, velx, vely, velz, velr, velp, emag, eint
   gc.collect()


if yt.is_root():
   results = list(my_storage.values())
   def plot_fig(field, title, ylabel, ymin, ymax, ylog):
      plt.figure(dpi = figure_dpi)
      for r in results:
         idx  = r["idx"]
         time = r["time"]
         plt.plot(r["mean_r"], r["%s"%field], color = cm.Blues(0.3+0.6*(idx-idx_start)/(idx_end+1-idx_start)), label ='t=%3d Myr'%(np.round(time)))
      plt.xlim((0, radius))
      if ymin != None and ymax != None:
         plt.ylim((ymin, ymax))
      if ylog:
         plt.yscale('log')
      plt.grid(ls='--')
      plt.legend(bbox_to_anchor=(1.25, 0), loc='lower right', borderaxespad=0, shadow=True, prop={'size':8})
      plt.xlabel(r"$R$ (kpc)")
      plt.ylabel( ylabel )
      plt.savefig("fig__%s.png"%title, dpi = figure_dpi, bbox_inches="tight")
      plt.close()

   plot_fig( 'mean_vp',  'rotation_curve',  r"$v_{\rm cir}$ (km/s)",                        0,   80, False )
   plot_fig( 'surdens',  'surface_density', r"$\Sigma$ (${\rm M}_{\odot}/{\rm kpc}^2$)",  1e2,  1e8, True  )
   plot_fig( 'sigma_vz', 'sigma_vz',        r"$\sigma_z$ (${\rm km}/{\rm s}$)",             0,   30, False )
   plot_fig( 'plasma_b', 'plasma_beta',     r"$\beta$",                                  None, None, True  )
   plot_fig( 'mean_B',   'mean_B',          r"$\langle B \rangle$ (G)",                  None, None, True  )


