Code organization

The repository is organized into three main components: Functions, Simulation, and Application, corresponding to method implementation, simulation studies, and real data analysis presented in the manuscript.




Functions/

This directory contains the core functions used to implement the proposed ITree-MMRM method, as well as auxiliary routines used in the simulation studies.

Functions/Functions.R
Contains the main functions implementing the proposed ITree-MMRM algorithm.

Functions/IT_code.R
Implements the Interaction Tree (IT) method proposed by Su (2009), which is used as a comparison method in Setting (vi) of the simulation study.

Functions/cal_nonnull.R
Computes the proportion of null trees obtained from simulation results.

Functions/cal_size.R
Computes tree size summaries based on simulation outputs.





Simulation/

This directory contains all simulation studies reported in the manuscript.

Simulation/Setting1-2/
Simulation code corresponding to Table 1 in the manuscript.

Simulation/Setting3-5/
Simulation code corresponding to Table 2 in the manuscript.

Simulation/Setting6/
Simulation code corresponding to Table 3 in the manuscript.

Each setting folder includes scripts to generate data, fit competing methods, and summarize performance metrics reported in the paper.





Application/

This directory contains the real data analysis based on the Alzheimer’s disease study.

Application/Application_update.R
Performs data preprocessing and subgroup analysis for the Alzheimer’s application, corresponding to results reported in Table 5 of the manuscript.

Application/Bootstrap_outcome.RData
Stores the final estimated tree structures and the corresponding pruning sequences under different values of the tuning parameter λ.

Application/merge_data2_1.sas7bdat
The data underlying this article will be shared on reasonable request to the corresponding author.
