import yt


def load_WLMDwarfGalaxy_datasets(code, prefix, idx_start, idx_end, didx):

   if code == 'GAMER':
      ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
   elif code == 'GIZMO':

      bbox = [ [-0.5*450.0, 0.5*450.0],
               [-0.5*450.0, 0.5*450.0],
               [-0.5*450.0, 0.5*450.0] ]

      unit_base = {
                    'UnitLength_in_cm'        : 3.085678e+21,
                    'UnitMass_in_g'           : 1.989e+43,
                    'UnitVelocity_in_cm_per_s': 100000,
                  }
      ts = yt.DatasetSeries( [ prefix+'/snap_%03d.hdf5'%idx for idx in range(idx_start, idx_end+1, didx) ],
                             unit_base=unit_base, bounding_box=bbox )
   else:
      raise RuntimeError('Code %s is NOT supported  !!'%code)

   return ts
