# ref: https://yt-project.org/docs/dev/cookbook/calculating_information.html#using-particle-filters-to-calculate-star-formation-rates

import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import WLMDwarfGalaxy_derived_fields
import h5py
from yt.units.yt_array import YTArray
from yt.data_objects.particle_filters import filter_registry

filein1 = "File1"
filein2 = "File2"
fileout = "fig__SN_feedback_rate_compare"
nbin    = 100
dpi     = 150


title1  = 'Label1'
title2  = 'Label2'

filein  = {
             title1 : filein1,
             title2 : filein2,
          }
y_loc   = {
             title1 : 4.0e+2,
             title2 : 2.0e+2,
          }
y_loc_sfr={
             title1 : 1.0e-2,
             title2 : 8.0e-3,
          }
color   = {
             title1 : cm.Reds (np.linspace(0.5, 0.9, 2)),
             title2 : cm.Blues(np.linspace(0.5, 0.9, 2)),
          }
SFR = True

# load data
for flu in [title1, title2]:
   ds = yt.load( filein[flu] )


   WLMDwarfGalaxy_derived_fields.set_particle_types('GAMER')

   WLMDwarfGalaxy_derived_fields.set_derived_fields(ds)

   # define the particle filter for the exploded SNe
   def unexploded_SNe( pfilter, data ):
      filter = (data[ "all", "ParSNIITime" ] > 0) & (data[ "all", "ParSNIITime" ] < data.ds.current_time) # ParSNIITime is set as negative value of explosion time for exploded particle
      return filter

   if 'unexploded_SNe' not in filter_registry:
      yt.add_particle_filter( "unexploded_SNe", function=unexploded_SNe, filtered_type="all", requires=["ParSNIITime"] )
   ds.add_particle_filter( "unexploded_SNe" )


   # get the creation time of the new stars and explosion time of teh SNe
   ad            = ds.all_data()

   star_ones     = ad[ "new_star", "particle_ones" ]
   mass          = ad[ 'new_star', 'ParMass'    ].in_units( 'Msun' )
   # add back the feedback mass
   exploded      = (ad[ 'new_star', 'ParSNIITime'  ] <=0 )
   mass[exploded] += ds.quan( ds.parameters['FB_ResolvedSNeII_EjectMass'], 'code_mass' ).in_units('Msun')


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
   assert np.all( SNe_upper_idx  > 0 ) and np.all( SNe_upper_idx  < len(t_bin) ), "incorrect SNe_upper_idx !!"


   # calculate the star formation and SNe number rate
   if SFR:
      Myr2yr = 1.0e6
      sfr    = np.array(  [ mass     [star_upper_idx == j+1].sum() / ( (t_bin[j+1] - t_bin[j])*Myr2yr ) for j in range(len(time)) ]  )
   else:
      sfr    = np.array(  [ star_ones[star_upper_idx == j+1].sum() / ( (t_bin[j+1] - t_bin[j])        ) for j in range(len(time)) ]  )

   sfr[sfr == 0] = np.nan

   SNr     = np.array(  [ SNe_ones[SNe_upper_idx == j+1].sum()  / ( (t_bin[j+1] - t_bin[j]) ) for j in range(len(time)) ]  )
   SNr[SNr == 0] = np.nan


   # calulate the conversion factor
   StarsPerSN  = 1.0/(ds.parameters['FB_ResolvedSNeII_NPerMass']*ds.parameters['SF_CreateStar_MinStarMass'])
   SNDelayTime = ds.quan( ds.parameters['FB_ResolvedSNeII_DelayTime'], 'code_time' ).in_units('Myr').d

   # plot
   if SFR:
      plt.plot( time,               sfr,                  label='%-6s'%flu, color=color[flu][0] )
   else:
      plt.plot( time,               sfr,                  label='Stars, %-6s'%flu, color=color[flu][0] )
      plt.plot( time,               SNr,                  label='SNe,   %-6s'%flu, color=color[flu][1] )
      #plt.plot( time,               SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f'%(StarsPerSN) )
      #plt.plot( time.d-SNDelayTime, SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f, shifted %.1f Myr'%(StarsPerSN, SNDelayTime) )
   plt.yscale('log')
   plt.xlim( 0.0, 800 )
   if SFR:
      plt.ylim( 1.0e-4, 2.0e-2 )
      plt.ylabel( '$\mathrm{SFR\ [M_\odot yr^{-1}]}$', fontsize='large' )
      text_string = '%-6s: Total number of formed stars   = % 7d\n'%(flu, len(creation_time))
      plt.text( 10.0, y_loc_sfr[flu], text_string, fontfamily='monospace' )

   else:
      plt.ylim( 1.0, 2.0e+3 )
      plt.ylabel( "$\mathrm{[Myr^{-1}]}$", fontsize="large" )
      text_string = '%-6s: Total number of formed stars   = % 7d\n'%(flu, len(creation_time)) + \
                    '%-6s: Total number of exploded SNe   = % 7d\n'%(flu, len(SNe_expl_time))
      plt.text( 10.0, y_loc[flu], text_string, fontfamily='monospace' )

   plt.legend(loc='lower right')
   plt.xlabel( "$\mathrm{t\ [Myr]}$",  fontsize="large" )


# show/save figure
plt.savefig( fileout+".png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )
#plt.show()
plt.close()


