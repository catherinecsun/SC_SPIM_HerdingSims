# SC_SPM_sims
Evaluating violations of independence (among indivdiuals, in their location and movement patterns) in SCR-based unmarked density estimation models.

- N (pop size): 140 
- sigma (range of movement):3
- Aggregation (group size): [1,4,10]
- Cohesion: [0,0.3, 0.67, 1]
- Detection probability: [0.05, 0.2]
- Sampling occasions (K): 4
- Camera Traps (J): 75 (15x 5)
- Statespace: 

- 24 population scenarios (but only 18 relevant ones given all combos with Cohesion===0 are the equivalent)
- 100 populations per scenario to simulate/generate data for, and to then estimate.

## Scripts
* Scripts with 1 generate and explore data (e.g., calculate Fletcher's c-hat)
* Scripts with 2 run the Spatial Count (Chandler and Royle) and Spatial Partial Identity (Augustine) models
  * uses the Microsoft Azure VM
  * Bayesian framework
      * uses informed gamma distribution for sigma prior
      * Data augmentation parameter M=400/600
      * single chain of 49998 iterations, burn in of the initial of 4998 and then split into 3 equal chains
* Scripts with 3 postprocess and summarize
* Script 4 calculates Relative Variance (Variance compared to the variance under independence), which can then be compared to c-hat if desired, although not very relevant)
* Script 5 takes subsets of iterations and plots them for visual determination of the stabiltiy of results (given that only 100 iterations/datasets were used)

## Future opportunities
* More PID values
* Missing partial identity covariates
* Vary trap spacing
* consider non-SCR based models

