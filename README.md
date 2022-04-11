# Effective population size in SIR models

Madi Yerlanov, Caroline Colijn, and Jessica E Stockdale

The given repository is created to accomany the correspoding paper. It contains the relevant data and codes to analyse this data. The programming languages being used are Julia (for numerical computation) and R (for result analysis). The codes also include the building and application of the main method used in this paper - fixed gamma. Relevant comments are present to provide descriptions, explanations and warnings. Further details such as theory, results and discussion are in the paper.

The main data is COVID-19 outbreak in China (selected cities). "Outbreak data" folder contains files City_Confirmed_0115_0816_infected.csv (active infections) and City_Confirmed_0115_0816_recovered.csv (total recovered) that was obtained from Harvard Database. The columns represent dates in format MM/DD. The rows represent cities. Note that the numbers represent current situations and not new infections/removals. Any data that is used must be in similar format.

The data (simulated or real) in .csv is analysed using FixedGamma.jl, where relative parameters are computed using fixed gamma method. These parameters include the effective population size (N*), infection rate (beta), etc. Note that simulated data do not have any headings for rows, so there are some difference in the way how you work with data. 

MappingOutbreak.r is used to locate selected Chinese cities. The cities are clustered according to the number of maximum infected (Imax) in a day. Red denotes that outbreaked reached over 1000 currently infected in at least one day. Similarly, the other colors denote, orange: 250-999, yellow: 100-249 and green: 50-99 respectively. The code produces figure 2 in paper.

SimulatingData.r allows to produce simulated outbreak data that FixedGamma.jl can work with. There are two types of simulations. The first one used the Gillepsie algorithm and it is the main one for this paper. Given rates and the population size, the code gives the infected and removed data. We have 6 settings for which we test the fixed gamma method. The second simulation uses survival dynamical system introduced by KhudaBukhsh et al. Along with the simulation, the population size can be computed. In this manner, two methods can be compared. Note that this repository do not contain the simulated data, but anyone can generate it using this code.

PlotAnalysis.r contains several codes to produce various plots (Figure 1, 3-7 and S1, S2). It mainly uses data and the results of FixedGamma.jl, i.e. N* and other parameter estimates. However, when analysing the real data, one needs other relevant information. In our case, there are tables containing information on census population size of selected, Chinese cities and distances to Wuhan. For further analysis, other data points must also be present in the table format such Imax and final size (Rinfinity). The clustering according to Imax (thus matching MappingOutbreak.jl) is done whenever possible. If there are different settings involved, then each setting data must be analysed separately and then combined in a grid. The code does not produce any number, such as averages, however, these computations can be performed along the way using few commands.

any questions can be addressed to madi dot yerlanov dot colorado dot edu
