from __future__ import absolute_import, unicode_literals, print_function

from EventStore import EventStore

if __name__ == '__main__':

  FILENAME = '/gpfs02/eic/lkosarzew/Calorimetry/nHcalDev/data/nhcal_sim.edm4hep.root'
  store = EventStore([FILENAME])  # pylint: disable=invalid-name # too strict before 2.5.0
  for i in range(store.current_store.getEntries()):
      store.current_store.endOfEvent()
      print('reading event', i)
      hits = store.get("HcalEndcapNHits")
      for hit in hits:
          print(' Referenced hit has an energy of', hit.getEnergy())