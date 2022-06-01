# Effective population size in SIR models

Madi Yerlanov, Caroline Colijn, and Jessica E Stockdale

This repository accompanies the above paper, and contains all relevant data and code to replicate the analyses within. The programming languages used are Julia (for numerical computation) and R (for analysis of results).

## COVID-19 outbreaks in China - 'Chinese Cities Data'

This folder contains data and analysis results on outbreaks of COVID-19 in 53 cities in China. Dta was obtained from Harvard Dataverse; we thank the owners of this data for making it publicly available, 

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
**true_sol.csv** - Results: solution to the determenistic system of equation for given settings

The data (simulated or real) in .csv is analysed using FixedGamma.jl, where relative parameters are computed using fixed gamma method. These parameters include the effective population size (N*), infection rate (beta), etc. Note that simulated data do not have any headings for rows, so there are some difference in the way how you work with data. 

MappingOutbreak.r is used to locate selected Chinese cities. The cities are clustered according to the number of maximum infected (Imax) in a day. Red denotes that outbreaked reached over 1000 currently infected in at least one day. Similarly, the other colors denote, orange: 250-999, yellow: 100-249 and green: 50-99 respectively. The code produces figure 2 in paper.

SimulatingData.r allows to produce simulated outbreak data that FixedGamma.jl can work with. There are two types of simulations. The first one used the Gillepsie algorithm and it is the main one for this paper. Given rates and the population size, the code gives the infected and removed data. We have 6 settings for which we test the fixed gamma method. The second simulation uses survival dynamical system introduced by KhudaBukhsh et al. Along with the simulation, the population size can be computed. In this manner, two methods can be compared. Simulated outbreak data is contained in folder "Simulated data".

PlotAnalysis.r contains several codes to produce various plots (Figure 1, 3-7 and S1, S2). It mainly uses data and the results of FixedGamma.jl, i.e. N* and other parameter estimates. However, when analysing the real data, one needs other relevant information. In our case, there are tables containing information on census population size of selected, Chinese cities and distances to Wuhan. For further analysis, other data points must also be present in the table format such Imax and final size (Rinfinity). The clustering according to Imax (thus matching MappingOutbreak.jl) is done whenever possible. If there are different settings involved, then each setting data must be analysed separately and then combined in a grid. The code does not produce any number, such as averages, however, these computations can be performed along the way using few commands.

any questions can be addressed to madi dot yerlanov dot colorado dot edu
