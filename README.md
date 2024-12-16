Step 1:
Clone the repository and install eic-shell and epic. These can be found here along with instructions.

https://github.com/eic/epic/tree/3fceda3555c0f537d6158e0e3fc569c3d0542b42
https://eicrecon.epic-eic.org/#/tutorial/02-work-environment

Step 2:
On command line, or using whichever job scheduler is relevant, run the following set of commands:

**The main directory of the project gives an example of how to do this using slurm for job submission.**

```/path/to/your/eic-shell```
```source /path/to/your/epic/install/setup.sh```
```npsim --compactFile=/path/to/your/epic/install/share/epic/epic_backward_hcal_only.xml -N=$N_EVT --enableGun --gun.particle="neutron" --gun.energy $ENERGY*GeV --gun.thetaMin 130*deg --gun.thetaMax 177*deg --gun.distribution uniform --outputFile $OUT_FILE```

Where $N_EVT is the number of events you want to simulate, $ENERGY is the energy of the neutron beam, and $OUT_FILE is the name of the output file.
Run the command for each of the following energies energies in GeV for the gun energy: 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1, 1.5, 2, 3, 4, 5.

I did this with 50 jobs of 2000 events each for each energy.

Step 3:
Once those are complete, run the ```ePICSimDataAnalysis/readTreeSimMain``` (after compiling with make) on the list of output files for a particular energy. Run it in this way: ```./readTreeSimMain <file_list> [output_file] [nevents]```

**The bash scripts in the main directory of the project gives an example of how to do this using slurm for job submission.**

The output files should go to the same directory as ```plotting.C```. They should have the following names:
backhcal_0.1gev_batch1.root
backhcal_0.2gev_batch1.root
backhcal_0.3gev_batch1.root
backhcal_0.4gev_batch1.root
backhcal_0.5gev_batch1.root
backhcal_0.65gev_batch1.root
backhcal_0.8gev_batch1.root
backhcal_1gev_batch1.root
backhcal_1.5gev_batch1.root
backhcal_2gev_batch1.root
backhcal_3gev_batch1.root
backhcal_4gev_batch1.root
backhcal_5gev_batch1.root

Step 4:
Finally, run ```plotting.C```
