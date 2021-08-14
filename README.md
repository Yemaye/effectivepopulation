# effectivepopulation

The real COVID-19 outbreak in selected Chinese cities are in files City_Confirmed_0115_0816_infected.csv and City_Confirmed_0115_0816_recovered.csv. 
These data is used in RealDataAnalysis.jl to produce data-frame on computed parameter values (.csv). FixedPopulationSize.jl allows to see how fitting using the census population size
fails.

MappingOutbreak.r is used to located selected Chinese cities clustered according to maximum infected in a day.

SimpleSimulation.r and ComplexSimulation.r produce simulated outbreak data. Vary parameters to obtain various settings. The resulted .csv is used in SimulDataAnalysis.jl (primary) or
PartialDataAnalysis.jl (cuts data at a given time) to produce .csv. OutlierAnalysis.jl extends analysis by detecting and identifying outliers in parameter data-frames. 

The produced .csv from SimulDataAnalysis.jl is used in BoxPlotAnalysis.r to give box-plots for parameters with mean and median (over 1000 simulations). For .csv produced by 
RealDataAnalysis.jl, one has to use NonBoxPlotAnalysis.r to plot optimum values with 95 confidence intervals (instead of box-plots), clustered. Correlation plots are also produced by this file.
Requires data-frame on the census population and factor in consideration: density, distance etc (china_pop2.csv). 
