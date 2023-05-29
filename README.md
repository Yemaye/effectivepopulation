# Effective population size in simple infectious disease models

Madi Yerlanov, Piyush Agarwal, Caroline Colijn, and Jessica E Stockdale

This repository accompanies the above paper, and contains all relevant data and code to replicate the analyses within. The programming languages used are Julia and Python (for numerical computation) and R (for analysis of results).

## COVID-19 outbreaks in China - 'Chinese Cities Data'

This folder contains data and analysis results on outbreaks of COVID-19 in 53 cities in China. Data was obtained from Harvard Dataverse; we thank the owners of this data for making it publicly available, 

**Harvard Dataverse: China Data Lab, 2020, "China COVID-19 Daily Cases with Basemap", https://doi.org/10.7910/DVN/MR5IJN, Harvard Dataverse, V38, UNF:6:BmhcC5NqO9pMzMyBUcfmtQ== [fileUNF]**.  

Files:  
**City_Confirmed_0115_0816_infected.csv**  - Data: daily active case counts per city  
**City_Confirmed_0115_0816_recovered.csv** - Data: daily total recovered counts per city   
Column headers represent dates in MM/DD format. If you are applying the code in this repository to other datasets, note that the numbers here represent total active cases and removals, not incident.  
**china_pop.csv** - Data: Census population size of each city, in millions. Final column contains a combined format (City, CPS) used for plotting.  
**china_pop2.csv** - Data: Further city-level information used for plotting: census population size, distance to Wuhan, maximum daily active cases I_max and final size R_infinity.  
**China_Betas.csv** - Results: estimates of parameter beta per city
**China_Ns.csv** - Results: estimates of parameter N* per city
**China_Ns_lower.csv** - Results: lower confidence interval for parameter N* per city
**China_Ns_upper.csv** - Results: upper confidence interval for parameter N* per city

## Simulated outbreaks - 'Simulated Data'

This folder contains simulated outbreak data, for the main text analysis where 1000 outbreaks are simulated under 6 scenarios and for the supplementary comparison with the SDS approach. 

Files:  
**infec_####_##.csv** - Data: contains 1000 simulated outbreaks for the simulation scenario described by the #. The first 3-4 digits represent the true population size (500-1000-5000). The second pair of digits represent the basic reproduction number (1.5,2.5,3.5).
**true_sol.csv** - Deterministic system solution for the 6 scenarios simulated above.
**Sellke_infec_####_##.csv** -  Data: contains simulated outbreaks using SDS approach for supplementary analysis. # has the same meaning as above. Note, only the first 20 outbreaks in each csv were used.

## Analysis scripts - 'FixedGammaChina.jl/FixedGammaSimulated.jl'
These scripts perform the data (simulated or China COVID-19) analysis: fitting an SIR model using least squares to estimate the effective population size (N*) and infection rate (beta), given a fixed value of removal rate gamma. csv files containing the parameter estimates are output. 

## Plotting script - PlotAnalysis.r 
This script contains code to produce the manuscript plots (Figures 1, 3-7, S1, S2), using data and the output csv files generated in the above FixedGamma.jl scripts. 

## Simulating outbreaks - SimulatingData.r 
This script produces the simulated outbreak data in folder 'Simulated Data', for use in FixedGammaSimulated.jl. There are two types of simulations: the first uses the Gillepsie algorithm for the simulation study in the main text of this paper, the second uses the survival dynamical system approach introduced by KhudaBukhsh et al. [https://doi.org/10.1098/rsfs.2019.0048] used in the supplementary analysis.

## Mapping of China data - MappingOutbreak.r
This script creates a map for the 53 Chinese cities used in the COVID-19 analysis, where cities are coloured according to their maximum daily active cases (Imax) - figure 2 in paper.

any questions can be addressed to madi dot yerlanov dot colorado dot edu
