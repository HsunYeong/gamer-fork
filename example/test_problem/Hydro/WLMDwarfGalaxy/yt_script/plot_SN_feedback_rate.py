# ref: https://yt-project.org/docs/dev/cookbook/calculating_information.html#using-particle-filters-to-calculate-star-formation-rates

import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt
import WLMDwarfGalaxy_derived_fields
import h5py
from yt.units.yt_array import YTArray


filein  = "../Data_000075"
fileout = "fig__SN_feedback_rate"
nbin    = 100
dpi     = 150


# load data
ds = yt.load( filein )


WLMDwarfGalaxy_derived_fields.set_particle_types('GAMER')

WLMDwarfGalaxy_derived_fields.set_derived_fields(ds)

# define the particle filter for the exploded SNe
def unexploded_SNe( pfilter, data ):
   filter = (data[ "all", "ParSNIITime" ] > 0) & (data[ "all", "ParSNIITime" ] < data.ds.current_time) # ParSNIITime is set as negative value of explosion time for exploded particle
   return filter

add_particle_filter( "unexploded_SNe", function=unexploded_SNe, filtered_type="all", requires=["ParSNIITime"] )
ds.add_particle_filter( "unexploded_SNe" )


# get the creation time of the new stars and explosion time of teh SNe
ad            = ds.all_data()

star_ones     = ad[ "new_star", "particle_ones" ]
creation_time = ad[ "new_star", "ParCreTime" ].in_units( "Myr" )

SNe_ones        = ad[ "exp_SNII", "particle_ones" ]
SNe_expl_time   = -1.0*( ad[ "exp_SNII", "ParSNIITime" ]*ds.units.code_time ).in_units( "Myr" ) # ParSNIITime is set as negative value of explosion time for exploded particle
SNe_unexpl_time =    ( ad[ "unexploded_SNe", "ParSNIITime" ]*ds.units.code_time ).in_units( "Myr" ) # ParSNIITime is set as negative value of explosion time for exploded particle


# bin the data
t_start        = 0.0
t_end          = ds.current_time.in_units( "Myr" )
t_bin          = np.linspace( start=t_start, stop=t_end, num=nbin+1 )
time           = 0.5*( t_bin[:-1] + t_bin[1:] )
star_upper_idx = np.digitize( creation_time, bins=t_bin, right=True )
SNe_upper_idx  = np.digitize( SNe_expl_time, bins=t_bin, right=True )


assert np.all( star_upper_idx > 0 ) and np.all( star_upper_idx < len(t_bin) ), "incorrect star_upper_idx !!"
assert np.all( SNe_upper_idx > 0 )  and np.all( SNe_upper_idx < len(t_bin) ),  "incorrect SNe_upper_idx !!"


# calculate the star formation and SNe number rate
sfr    = np.array(  [ star_ones[star_upper_idx == j+1].sum() / ( (t_bin[j+1] - t_bin[j]) ) for j in range(len(time)) ]  )
sfr[sfr == 0] = np.nan

SNr     = np.array(  [ SNe_ones[SNe_upper_idx == j+1].sum()  / ( (t_bin[j+1] - t_bin[j]) ) for j in range(len(time)) ]  )
SNr[SNr == 0] = np.nan


# calulate the conversion factor
StarsPerSN  = 1.0/(ds.parameters['FB_ResolvedSNeII_NPerMass']*ds.parameters['SF_CreateStar_MinStarMass'])
SNDelayTime = ds.quan( ds.parameters['FB_ResolvedSNeII_DelayTime'], 'code_time' ).in_units('Myr').d

# plot
plt.plot( time,               sfr,                  label='Stars' )
plt.plot( time,               SNr,                  label='SNe'   )
#plt.plot( time,               SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f'%(StarsPerSN) )
#plt.plot( time.d-SNDelayTime, SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f, shifted %.1f Myr'%(StarsPerSN, SNDelayTime) )
plt.yscale('log')
plt.xlim( 0.0, 800 )
plt.ylim( 1.0, 2.0e+3 )
plt.legend()
plt.xlabel( "$\mathrm{t\ [Myr]}$",  fontsize="large" )
plt.ylabel( "$\mathrm{[Myr^{-1}]}$", fontsize="large" )

text_string = 'Total number of formed stars   = % 7d\n'%len(creation_time) + \
              'Total number of exploded SNe   = % 7d\n'%len(SNe_expl_time) + \
              'Total number of unexploded SNe = % 7d\n'%len(SNe_unexpl_time)
plt.text( 10.0, 4.0e+2, text_string, fontfamily='monospace' )


# show/save figure
plt.savefig( fileout+".png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )
#plt.show()
plt.close()


# get total dust mass
f = h5py.File('../Data_000075', 'r')
Unit_T = f['Info']['InputPara']['Unit_T']
Unit_M = f['Info']['InputPara']['Unit_M']
f.close()

def IdxSearch(target, arr):
   left  = 0
   right = len(arr) - 1
   while left < right:
      mid = (left + right) // 2
      if arr[mid] < target:
         left = mid + 1
      elif arr[mid] > target:
         right = mid -1
      else:
         return mid
   return left

table = np.loadtxt("../Record__Conservation")
time_all    = YTArray( table[:, 0]*Unit_T, 's').in_units('Myr').d
dust_mass   = YTArray( table[:,47]*Unit_M, 'g').in_units('Msun').d

SNe_upper_idx_time_step = np.digitize( SNe_expl_time, bins=time_all, right=True )
SNe_time_step           = np.array(  [ SNe_ones[SNe_upper_idx_time_step == j+1].sum()  for j in range(len(time_all)) ]  )
SNe_sum                 = np.cumsum(SNe_time_step)
idx_300 = IdxSearch(300, time_all)
idx_750 = IdxSearch(750, time_all)


SNe_300_750 = SNe_sum  [idx_750] - SNe_sum  [idx_300]
d_dust_mass = dust_mass[idx_750] - dust_mass[idx_300]

d_dust_per_SNe = d_dust_mass / SNe_300_750
print("Dust mass decreased during 300-750 Myr = " + str(d_dust_mass) + " M_sun")
print("Total SNe during 300-750 Myr = " + str(SNe_300_750))
print("Dust mass decrease per SN = " + str(d_dust_per_SNe) + " M_sun")


plt.plot( time_all, dust_mass,      label=r'$\rho_{\rm dust}/\rho_{\rm dust}(\mathrm{t}=0)$' )
plt.plot( time_all, SNe_time_step,  label='Number of SNe'   )
plt.yscale('log')
plt.xlim( 0.0, 800 )
plt.legend(loc='upper left')
plt.xlabel( "$\mathrm{t\ [Myr]}$",  fontsize="large" )

# show/save figure
plt.savefig( "fig_dust_vs_SNe.png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )

