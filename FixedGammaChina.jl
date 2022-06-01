#Packages that are used/relevant
using Plots
using DifferentialEquations
using LsqFit
using DataFrames
using CSV
using LinearAlgebra
using Optim
using Roots
using StatsBase
using LaTeXStrings

#Real data on Chinese cities. The first two entries are the name of a city and a prefecture
df=DataFrame(CSV.File("City_Confirmed_0115_0816_infected.csv")) #I-data
df1=Matrix(df)
#dg=DataFrame(CSV.File("City_Confirmed_0115_0816_total.csv")) #I+R-data
#dg1=Matrix(dg)

#Important functions
function f(du,u,p,t) #Sets up a system of differential equations to be solved
    du[1]=dx=p[1]*(p[3]-u[1]-u[2])*u[1]-p[2]*u[1]
    du[2]=dy=p[2]*u[1]
end

function sir(datapoints) #Returns the indices when the outbreak started and ended
    ind1=findfirst(x->x>0,datapoints)
    ind2=findfirst(x->x==0,datapoints[ind1+10:length(datapoints)])
    return (ind1,ind1+9+ind2)
end

function lss(u,time, phi) #Extracts I-data from the solution of the system, created by function f
    tspan=(time[1],time[end])
    prob=ODEProblem(f,u,tspan,phi)
    sol=solve(prob,Vern9(),saveat=time)
    estimated=reduce(hcat,sol.u)
    return estimated[1,:]
end

n=53 #Number of data/number of rows
Niv=1000.0 #Initial point for the pop.size
Betaiv=0.0005 #Initial point for the beta
Gammaiv=1/6 #Intial point for the gamma (fixed in this case)
p1=[Betaiv,Gammaiv,Niv]

#Similar as above, but the data shifted so the non-numerical values are not taken
L=zeros(n,8)
p2=[p1[3],p1[1]]
for i in 1:n
    data=df1[i,3:200] #here shifts happens, change appropriately for your data
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) 
            du[1]=dx=p[2]*(p[1]-u[1]-u[2])*u[1]-1/6*u[1]
            du[2]=dy=1/6*u[1]
    end
    function lss(u,time, phi) 
        tspan=(time[1],time[end])
        prob=ODEProblem(h,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    lssk(ti,pi)=lss(idata,ti,pi)
    fir=curve_fit(lssk,tdata,ndata,p2)
    L[i,1:2]=fir.param
    L[i,3]=confidence_interval(fir,0.05)[1][1]
    L[i,4]=confidence_interval(fir,0.05)[1][2]
    L[i,5]=confidence_interval(fir,0.05)[2][1]
    L[i,6]=confidence_interval(fir,0.05)[2][2]
    L[i,7]=dot(fir.resid,fir.resid)
    L[i,8]=ind[2]
end

#Saving the results, in this case, eff. pop. size
Neff=zeros(n,1)
Neff[:,1]=L[:,1]
CSV.write("China_Ns.csv",  DataFrame(Neff), writeheader=false)

Nlow=zeros(n,1)
Nlow[:,1]=L[:,3]
CSV.write("China_Ns_lower.csv",  DataFrame(Nlow), writeheader=false)

Nup=zeros(n,1)
Nup[:,1]=L[:,4]
CSV.write("China_Ns_upper.csv",  DataFrame(Nup), writeheader=false)

Betas=zeros(n,1)
Betas[:,1]=L[:,2]
CSV.write("China_Betas.csv",  DataFrame(Betas), writeheader=false)

