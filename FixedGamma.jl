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

#Simulated data, change the path,names
di=DataFrame(CSV.File("simulations\\infec_5000_25.csv",header=false)) #data on daily infected, I-data
dr=DataFrame(CSV.File("simulations\\remov_5000_25.csv",header=false)) #data on daily removed, R-data
dim=Matrix(di) #transforming into matrices
drm=Matrix(dr)
dtm=dim.+drm #total daily affected

#Real data on Chinese cities. The first two entries are the name of a city and a prefecture
df=DataFrame(CSV.File("City_Confirmed_0115_0816_infected.csv")) #I-data
df1=Matrix(df)
dg=DataFrame(CSV.File("City_Confirmed_0115_0816_total.csv")) #I+R-data
dg1=Matrix(dg)

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

#Parameters
n=1000 #Number of data/number of rows
Niv=1000.0 #Initial point for the pop.size
Betaiv=0.0005 #Initial point for the beta
Gammaiv=0.2 #Intial point for the gamma (fixed in this case)
p1=[Betaiv,Gammaiv,Niv]

#Use for simulated data (where the columns do not contain non-numerical values such as city names)
L=zeros(n,8)
p2=[p1[3],p1[1]]
@time begin
for i in 1:n
    data=dim[i,:] #data
    ind=sir(data) #getting indices of start and end of an outbreak
    ndata=data[ind[1]:ind[2]] #cutting the data
    idata=[ndata[1],0.0] #initial data point
    tdata=collect(1.0:float(length(ndata))) #time vector
    function h(du,u,p,t)  #system set up
            du[1]=dx=p[2]*(p[1]-u[1]-u[2])*u[1]-0.2*u[1]
            du[2]=dy=0.2*u[1]
    end
    function lss(u,time, phi) #Modified function of extracting I-solution
        tspan=(time[1],time[end])
        prob=ODEProblem(h,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    lssk(ti,pi)=lss(idata,ti,pi) #Data-fixed function
    fir4=curve_fit(lssk,tdata,ndata,p3) #LSQ fitting
    L[i,1:2]=fir4.param #Returns N* and beta
    L[i,3]=confidence_interval(fir4,0.05)[1][1] #CI for N
    L[i,4]=confidence_interval(fir4,0.05)[1][2]
    L[i,5]=confidence_interval(fir4,0.05)[2][1] #CI for beta
    L[i,6]=confidence_interval(fir4,0.05)[2][2]
    L[i,7]=dot(fir4.resid,fir4.resid) #SSE
    L[i,8]=ind[2] #Last day of the outbreak
end
end

#Similar as above, but the data shifted so the non-numerical values are not taken
L=zeros(53,8)
p2=[p1[3],p1[1]]
for i in 1:53
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
    fir=curve_fit(lssk,tdata,ndata,p4)
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
CSV.write("Neff.csv",  DataFrame(Neff), writeheader=false)

#This will creates the true solution to the given setting
truet=collect(1.0:150.0)
truespan=(1.0,150.0)
trueprob=ODEProblem(f,[1,0],truespan,p1) 
truesol=solve(trueprob,Vern9(),saveat=truet)
trueplot=plot(truesol,vars=(0,1),xlims=[0.0,150.0],legend=false,color="black",linewidth=5)

#Plotting the first row in your data
data1=dim[1,:]
newdata1=float.(dim[1,sir(data1)[1]:sir(data1)[2]])
t1=collect(1.0:float(length(newdata1)))
tspan=(t1[1],t1[end])

#Plotting all data and the solution togethether
plot1=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],palette=:Dark2_8)
for j in 2:n #Adding other if needed
    datas=dim[j,:]
    newdata=float.(dim[j,sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
        plot!(plot8,t,newdata,legend=false,xlims=[0.0,150.0],palette=:Dark2_8)

end
plot1
plot!(truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",
linewidth=3, xlabel="Days since introduction",ylabel="Active infections",size=(600,400),guidefontsize=12)
#Addinf the true solution
savefig(plot1,"simulations/gamma_sim_25_5000.png")#Saving