# ePICRecoDataAnalysis
Analysis of ePIC simulated data. The code is written so that it can be easily converted between 3 different data reading modes:
- ROOT trees
- `podio::ROOTFrameReader`
- `podio::EventStore`

### To run interactively:

```Sh
root -l -b -q readTreeSim.C+ | tee run.log
```

### To compile executable:

```Sh
make
```

and Run:

```Sh
./readTreeSimMain | tee run.log
```

### Batch scripts to use condor on RCF

```Sh
submitSimAnalysis.job
runSimAnaBatch.sh
```
Submit with:

```Sh
condor_submit submitSimAnalysis.job
```

### Python version

You can use `ipython3` to have interactive `python` and autocomplete feauture

```Sh
python readEvents.py
```

### TO DO
1. Convert to `ROOTFrameReader` once bugs and memory leaks within `podio` are fixed
2. Similarly fix `EventStore` to read legacy simulation campaign data (pre-`ROOTFrameReader`)
3. Add more histograms and functionality
4. Add macros for drawing histograms