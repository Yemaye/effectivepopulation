#Showing how using real population does not work for SIR
#Procedure is the same as in the main code

using Plots
using DifferentialEquations
using LsqFit
using DataFrames
using CSV
using LinearAlgebra
using Optim

df=DataFrame(CSV.File("City_Confirmed_0115_0816_infected.csv"))
df1=convert(Matrix,df)

dg=DataFrame(CSV.File("City_Confirmed_0115_0816_total.csv"))
dg1=convert(Matrix,dg)

dp=DataFrame(CSV.File("china_pop2.csv",header=false))
dp1=convert(Matrix,dp)

#Important functions
function f(du,u,p,t) #IR diff eq
    du[1]=dx=p[1]*(p[3]-u[1]-u[2])*u[1]-p[2]*u[1]
    du[2]=dy=p[2]*u[1]
end

function sir(datapoints) #datapoints cutter
    ind1=findfirst(x->x>0,datapoints)
    ind2=findfirst(x->x==0,datapoints[ind1+10:length(datapoints)])
    return (ind1,ind1+9+ind2)
end

#Population size and gamma are fixed
K=zeros(53,4)
p5=[0.000000005]
for i in 1:53
    data=df1[i,3:200] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]] #new cutted data
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #IR We have to specify for each city separately
            du[1]=dx=(dp1[i,2]*10^6-u[1]-u[2])*u[1]*p[1]-1/6*u[1]
            du[2]=dy=1/6*u[1]
    end
    function lss(u,time, phi) #IR extract of I
        tspan=(time[1],time[end])
        prob=ODEProblem(h,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    lssk(ti,pi)=lss(idata,ti,pi)
    fir=curve_fit(lssk,tdata,ndata,p5)
    K[i,1]=fir.param[1]
    K[i,2]=confidence_interval(fir,0.05)[1][1]
    K[i,3]=confidence_interval(fir,0.05)[1][2]
    K[i,4]=dot(fir.resid,fir.resid)
end

#Creating all data and fitted curve plots
allplots=Matrix{Any}(nothing, 53,1)
for j in 1:53
    datas=df1[j,3:200]
    totaldatas=dg1[j,3:200]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    rdata=float.(totaldatas[sir(datas)[1]:sir(datas)[2]])-newdata
    u0=[newdata[1],0.0]
    t=collect(1.0:float(length(newdata)))
    tspan=(t[1],t[end])
    prob6=ODEProblem(f,u0,tspan,[K[j,1],K[j,4],dp1[j,2]*10^6]) #combo 6
    sol6=solve(prob6,Vern9(),saveat=t)
    plot6=plot(sol6,vars=(0,1),xlims=[t[1],t[end]],labels="Solution",ylabel="I, the number of infected", xlabel="t, days")
    plot62=scatter!(plot6,t,newdata,labels="Data",legend=:topleft)
    allplots[j,:]=[plot62]
end

#Choose 1
for i in 1:1
    display(allplots[i])
end

#If you want to compare with fixed gamma method, N* is not fixed
allplots1=Matrix{Any}(nothing, 53,1)
for j in 1:53
    datas=df1[j,3:200]
    totaldatas=dg1[j,3:200]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    rdata=float.(totaldatas[sir(datas)[1]:sir(datas)[2]])-newdata
    u0=[newdata[1],0.0]
    t=collect(1.0:float(length(newdata)))
    tspan=(t[1],t[end])
    prob6=ODEProblem(f,u0,tspan,[L[j,2],1/6,L[j,1]]) #combo 6
    sol6=solve(prob6,Vern9(),saveat=t)
    plot6=plot(sol6,vars=(0,1),xlims=[t[1],t[end]],labels="Solution", ylabel="I, the number of infected", xlabel="t, days")
    plot62=scatter!(plot6,t,newdata,labels="Data",legend=:topleft)
    allplots1[j,:]=[plot62]
end

#choose the same
for i in 39:39
    display(allplots1[i])
end

#fixed gamma optimum parameters
L[39,1:2]

#fixed gamma and fixed population optimum parameter (beta)
K[39,1]
savefig(allplots[1],"tianjin_real.png")

#save it
CSV.write("Beta_actual_pop.csv",  DataFrame(K), writeheader=false)



#if you want to compare with not fixed gamma
K=zeros(53,7)
p5=[0.00000005,0.15]
for i in 1:1
    data=df1[i,3:200] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]] #new cutted data
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #IR We have to specify for each city separately
            du[1]=dx=(dp1[i,2]*10^6-u[1]-u[2])*u[1]*p[1]-p[2]*u[1]
            du[2]=dy=p[2]*u[1]
    end
    function lss(u,time, phi) #IR extract of I
        tspan=(time[1],time[end])
        prob=ODEProblem(h,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    lssk(ti,pi)=lss(idata,ti,pi)
    fir=curve_fit(lssk,tdata,ndata,p5)
    K[i,1]=fir.param[1]
    K[i,2]=confidence_interval(fir,0.05)[1][1]
    K[i,3]=confidence_interval(fir,0.05)[1][2]
    K[i,4]=fir.param[2]
    K[i,5]=confidence_interval(fir,0.05)[2][1]
    K[i,6]=confidence_interval(fir,0.05)[2][2]
    K[i,7]=dot(fir.resid,fir.resid)
end


