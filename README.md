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
* Scripts with 1 generate and explore data
* Scripts with 2 run the Spatial Count (Chandler and Royle) and Spatial Partial Identity (Augustine) models
  * uses the Microsoft Azure VM
  * Bayesian framework
      * uses informed gamma distribution for sigma prior
      * Data augmentation parameter M=400/600
      * single chain of 49998 iterations, burn in of the initial of 4998 and then split into 3 equal chains
* Scripts with 3 postprocess and summarize

## Tasks

To do:
- for coverage plots, consider ribbons connecting the maxes to maxes and mins to mins, rather than a line segment through the median/mean
- reorder the factors for SC, SPIM: 4 , and SPIM: 3
- figure out a way to quantify the amount of ID info given pop size and partial id covs.. (read ben's paper and also L. Waits)
- connect with somebody about my PID calculations.


Finished:
- run antler +  sex +  coat +  collar
- run sex +  coat  + collar
- sex  + coat, SPIM
- collar + sex, SPIM

Run:
- collar, SPIM
- sex, SPIM
