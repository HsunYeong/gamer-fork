import argparse
from matplotlib import patheffects
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np

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
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
my_storage = {}

# main loop
for sto, ds in ts.piter(storage=my_storage):

   def _helicity( field, data ):
      return data["velocity_x"]*data["vorticity_x"] + data["velocity_y"]*data["vorticity_y"] + data["velocity_z"]*data["vorticity_z"]
   ds.add_field( ("gas", "helicity"), function=_helicity, sampling_type="cell", units="code_velocity/code_time" )

   dd   = ds.all_data()
   dens = dd["Dens"].d

#  total momentum
   wx   = dd["MomX"].d/dens**0.5
   wy   = dd["MomY"].d/dens**0.5
   wz   = dd["MomZ"].d/dens**0.5

#  remove center of mass motion
   vx = ( dd["MomX"].d - np.mean(dd["MomX"].d) )/dens**0.5
   vy = ( dd["MomY"].d - np.mean(dd["MomY"].d) )/dens**0.5
   vz = ( dd["MomZ"].d - np.mean(dd["MomZ"].d) )/dens**0.5

#  magnetic field
   bx   = dd['CCMagX'].d
   by   = dd['CCMagY'].d
   bz   = dd['CCMagZ'].d

#  mean magnetic_field
   mean_Bx = dd.quantities.weighted_average_quantity('magnetic_field_x', 'cell_volume')
   mean_By = dd.quantities.weighted_average_quantity('magnetic_field_y', 'cell_volume')
   mean_Bz = dd.quantities.weighted_average_quantity('magnetic_field_z', 'cell_volume')

   emag   = np.mean(0.5*(bx*bx + by*by + bz*bz))
   ek_tot = np.mean(0.5*(wx*wx + wy*wy + wz*wz))
   ekin   = np.mean(0.5*(vx*vx + vy*vy + vz*vz))
   heli   = np.mean(dd['helicity'].d)

   mach = (2*ekin)**0.5
   time = ds.current_time
   sto.result = {
        "time"   : time,
        "emag"   : emag,
        "ekin"   : ekin,
        "ek_tot" : ek_tot,
        "mach"   : mach,
        "heli"   : heli,
        "meanBx" : mean_Bx,
        "meanBy" : mean_By,
        "meanBz" : mean_Bz,
   }

# plot
if yt.is_root():
   time   = []
   emag   = []
   ekin   = []
   ek_tot = []
   mach   = []
   heli   = []
   meanBx = []
   meanBy = []
   meanBz = []

   for idx, val in sorted(my_storage.items()):
      time.append(   val["time"]   )
      emag.append(   val["emag"]   )
      ekin.append(   val["ekin"]   )
      ek_tot.append( val["ek_tot"] )
      mach.append(   val["mach"]   )
      heli.append(   val["heli"]   )
      meanBx.append( val["meanBx"] )
      meanBy.append( val["meanBy"] )
      meanBz.append( val["meanBz"] )

   time   = np.array(time)
   emag   = np.array(emag)
   ekin   = np.array(ekin)
   ek_tot = np.array(ek_tot)
   mach   = np.array(mach)
   heli   = np.array(heli)
   meanBx = np.array(meanBx)
   meanBy = np.array(meanBy)
   meanBz = np.array(meanBz)

   ratio = np.divide( emag, ekin,
           out=np.zeros_like(emag, dtype=float),
           where=ekin != 0 )
   # Plot
   f, ax = plt.subplots(3, 1, figsize=(6,10))
   f.subplots_adjust(wspace=0.4)
   ax[0].plot( time[1:], ratio[1:], label = r"$E_{\rm mag}/E_{\rm kin}$" )
   ax[0].set_yscale('log')
   ax[0].set_xlim(0, 100)
   ax[0].set_ylabel( r'$E_{\rm mag}/E_{\rm kin}$', fontsize='large' )

   ax[1].plot( time, mach, label = "Mach number")
   ax[1].set_xlim(0, 100)
   ax[1].set_ylim(0, 0.4)
   ax[1].set_ylabel('Mach number', fontsize='large')

   ax[2].plot( time, heli, label = "helicity")
   ax[2].set_xlim(0, 100)
   ax[2].set_xlabel(r"$t$", fontsize="large")
   ax[2].set_ylabel(r'Helicity', fontsize='large')
   plt.savefig("fig_EnergyEvolve.png", bbox_inches='tight', pad_inches=0.05, dpi=dpi)
   plt.close()

   # save figure
   np.savetxt( 'EnergyEvolve', np.column_stack( (time, emag, ekin, ek_tot, mach, heli, meanBx, meanBy, meanBz) ),
               fmt='  %16.8e',
               header='               t               Emag               Ekin             Ek_tot               Mach           helicity            mean_Bx            mean_By            mean_Bz' )


