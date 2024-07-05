# AVSCausalModel
R and Stan code used in data simulations together with the the generated data in the simulation study, as reported in "Applying causal inference and Bayesian statistics to understanding vaccine safety signals --- a simulation study" examining factors which influence the reporting of medically attended AEFI for further understanding of vaccine safety signals.

Factors influencing signal detection were investigated upon data generated under 12 scenarios differing in prevalence of a moderate to severe AEFI, survey participation and influence of a moderate to severe AEFI on survey participation and seeking medical attendance. The simulation code Data_generation_and_PPA_analysis.R includes the data generating process using parameters based upon literature and plausible assumptions as described in the manuscript and supplementary materials. The simulation code Figure_4_sim.R contains the data generating process for the single run of simulations used to visualise the contrasting proportions or Report MA between age groups and also demonstrate the effect of the differeng parameters governing the influence of survey participation and seeking medical attention. The individual scenario simulations for Figure 4 are uploaded as individual RDA files. The scenarios each differ in the event probabilities incorporated in the data simulations. The RDS file simulations.Rds contains the results of the PPA investigation involving 5,000 runs of the reference scenario (N = 50,000) and the 12 investigation scenarios (N = 4000). Please refer to the main text and supplementary materials for the details of the design of the scenarios.
