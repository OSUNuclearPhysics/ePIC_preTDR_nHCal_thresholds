from __future__ import absolute_import, unicode_literals, print_function

#from EventStore import EventStore
from podio.EventStore import EventStore

if __name__ == '__main__':

  #FILENAME = '/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/data/nhcal_sim.edm4hep.root'
  #FILENAME = '/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/130to177_degree/neutron_gun_1GeV_130to177_degree_sci_tiles_hcal_only_version2_637486_3.edm4hep.root'
  
  FILENAME = []
  FILENAME.append('/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/130to177_degree/neutron_gun_1GeV_130to177_degree_sci_tiles_hcal_only_version2_637486_15.edm4hep.root')
  #FILENAME.append('/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/130to177_degree/neutron_gun_1GeV_130to177_degree_sci_tiles_hcal_only_version2_637486_3.edm4hep.root')
  #FILENAME.append('/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/130to177_degree/neutron_gun_1GeV_130to177_degree_sci_tiles_hcal_only_version2_637486_4.edm4hep.root')
  #FILENAME.append('/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/130to177_degree/neutron_gun_1GeV_130to177_degree_sci_tiles_hcal_only_version2_637486_5.edm4hep.root')
  #FILENAME.append('/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/130to177_degree/neutron_gun_1GeV_130to177_degree_sci_tiles_hcal_only_version2_637486_7.edm4hep.root')
  
  store = EventStore(FILENAME)  # pylint: disable=invalid-name # too strict before 2.5.0
  #store = EventStore('/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/130to177_degree/neutron_gun_1GeV_130to177_degree_sci_tiles_hcal_only_version2_637486_2.edm4hep.root')  # pylint: disable=invalid-name # too strict before 2.5.0
  for i, event in enumerate(store):
      print('reading event', i)
      hits = store.get("HcalEndcapNHits")
      print('HcalEndcapNHits size = ', hits.size())
      for hit in hits:
          print(' Referenced hit has an energy of', hit.getEnergy())