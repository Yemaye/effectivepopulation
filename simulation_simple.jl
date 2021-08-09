#packages 
using Plots
using DifferentialEquations
using LsqFit
using DataFrames
using CSV
using LinearAlgebra
using Optim
using Roots
using StatsBase

#simulated data, change the names!
di=DataFrame(CSV.File("simulations\\infec_1000_25.csv",header=false))
dr=DataFrame(CSV.File("simulations\\remov_1000_25.csv",header=false))

di=DataFrame(CSV.File("complex_simulations\\com_infec_9_1.csv",header=false))
dr=DataFrame(CSV.File("complex_simulations\\com_remov_9_1.csv",header=false))

dim=convert(Matrix,di)
drm=convert(Matrix,dr)
dtm=dim.+drm

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

function lss(u,time, phi) #IR extract of I
    tspan=(time[1],time[end])
    prob=ODEProblem(f,u,tspan,phi)
    sol=solve(prob,Vern9(),saveat=time)
    estimated=reduce(hcat,sol.u)
    return estimated[1,:]
end

p1=[0.5/9000,0.2,9000.0] #starting point for parameters
N=zeros(1000,12) #fitting to I, other methods follow the same procedure
for i in 1:1000
    data=dim[i,:] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]] #new cutted data
    irdata=dtm[i,:] #total cases data
    r=log(irdata[ind[1]+9])/10 #growth rate
    rnot=1+6*r #r0
    idata=[ndata[1],0.0] #initial point for the system
    tdata=collect(1.0:float(length(ndata))) #time array
    lssk(ti,pi)=lss(idata,ti,pi) #new function with fixed i.p.
    fir=curve_fit(lssk,tdata,ndata,p1) #fitting
    N[i,1:3]=fir.param #saving beta gamma N
    N[i,4]=confidence_interval(fir,0.05)[1][1]
    N[i,5]=confidence_interval(fir,0.05)[1][2]
    N[i,6]=confidence_interval(fir,0.05)[2][1]
    N[i,7]=confidence_interval(fir,0.05)[2][2]
    N[i,8]=confidence_interval(fir,0.05)[3][1]
    N[i,9]=confidence_interval(fir,0.05)[3][2]
    N[i,6]=fir.param[1]*fir.param[3]/fir.param[2] #r0 based on the fitted
    N[i,11]=rnot #r0 via growth rate
    N[i,10]=dot(fir.resid,fir.resid) #residual squared
    N[i,12]=ind[2]
end


#Fitting on 3 parameters to R
M=zeros(1000,10)
for i in 1:1000
    data=dim[i,:] 
    ind=sir(data)
    rdata=drm[i,ind[1]:ind[2]]
    function rdiff(du,u,p,t)
        du[1]=p[2]*(p[3]-u[1]-(p[3]-data[ind[1]])*exp(-u[1]*p[1]/p[2]))
    end
    function lssr(u,time, phi) 
        tspan=(time[1],time[end])
        prob=ODEProblem(rdiff,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    t=collect(1.0:float(length(rdata)))
    lssrk(ti,pi)=lssr([0.0],ti,pi)
    fir2=curve_fit(lssrk,t,rdata,p1)
    M[i,1:3]=fir2.param
    M[i,4]=confidence_interval(fir2,0.05)[1][1]
    M[i,5]=confidence_interval(fir2,0.05)[1][2]
    M[i,6]=confidence_interval(fir2,0.05)[2][1]
    M[i,7]=confidence_interval(fir2,0.05)[2][2]
    M[i,8]=confidence_interval(fir2,0.05)[3][1]
    M[i,9]=confidence_interval(fir2,0.05)[3][2]
    M[i,10]=dot(fir2.resid,fir2.resid)
end

#final size
T=zeros(1000,2)

for i in 1:1000
    data=dim[i,:] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dtm[i,:] 
    ir10=irdata[ind[1]+9]
    r=log(ir10)/10 
    rnot=1+6*r
    si(N)=rnot-N/irdata[ind[2]]*log((N-1)/(N-irdata[ind[2]]))
    T[i,1]=find_zero(si,(irdata[ind[2]]+1,irdata[ind[2]]+100000),Bisection())
    T[i,2]=rnot
end

#Imax
R=zeros(1000,1)

for i in 1:1000
    data=dim[i,:] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dtm[i,:] 
    ir10=irdata[ind[1]+9]
    r=log(ir10)/10 
    rnot=1+6*r
    si(N)=findmax(ndata)[1]+N/rnot-N/rnot*log(N/rnot)-N+N/rnot*log(N-ndata[1])
    R[i]=find_zero(si,(findmax(ndata)[1]+1,irdata[ind[2]]+100000),Bisection())    
end

#beta via rnot
S=zeros(1000,7)
p2=[p1[3],p1[2]]
#p2=[1500,0.2]
for i in 1:1000
    data=dim[i,:] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dtm[i,:] 
    r=log(irdata[ind[1]+9])/10 
    rnot=1+6*r
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) 
            du[1]=dx=rnot*p[2]*(p[1]-u[1]-u[2])*u[1]/p[1]-p[2]*u[1]
            du[2]=dy=p[2]*u[1]
    end
    function lss(u,time, phi) 
        tspan=(time[1],time[end])
        prob=ODEProblem(h,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    lssk(ti,pi)=lss(idata,ti,pi)
    fir3=curve_fit(lssk,tdata,ndata,p2)
    S[i,1:2]=fir3.param
    S[i,3]=confidence_interval(fir3,0.05)[1][1]
    S[i,4]=confidence_interval(fir3,0.05)[1][2]
    S[i,5]=confidence_interval(fir3,0.05)[2][1]
    S[i,6]=confidence_interval(fir3,0.05)[2][2]
    S[i,7]=dot(fir3.resid,fir3.resid)
end

#gamma=0.2
L=zeros(1000,8)
p3=[p1[3],p1[1]]
for i in 1:1000
    data=dim[i,:] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) 
            du[1]=dx=p[2]*(p[1]-u[1]-u[2])*u[1]-0.2*u[1]
            du[2]=dy=0.2*u[1]
    end
    function lss(u,time, phi) 
        tspan=(time[1],time[end])
        prob=ODEProblem(h,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    lssk(ti,pi)=lss(idata,ti,pi)
    fir4=curve_fit(lssk,tdata,ndata,p3)
    L[i,1:2]=fir4.param
    L[i,3]=confidence_interval(fir4,0.05)[1][1]
    L[i,4]=confidence_interval(fir4,0.05)[1][2]
    L[i,5]=confidence_interval(fir4,0.05)[2][1]
    L[i,6]=confidence_interval(fir4,0.05)[2][2]
    L[i,7]=dot(fir4.resid,fir4.resid)
    L[i,8]=ind[2]
end



#beta via rnot, gamma fixed
K=zeros(1000,4)
p4=[p1[3]]
for i in 1:1000
    data=dim[i,:] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]] #new cutted data
    irdata=dtm[i,:] #total cases data
    r=log(irdata[ind[1]+9])/10 #growth rate
    rnot=1+6*r
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #IR We have to specify for each city separately
            du[1]=dx=rnot*0.2*(p[1]-u[1]-u[2])*u[1]/p[1]-0.2*u[1]
            du[2]=dy=0.2*u[1]
    end
    function lss(u,time, phi) #IR extract of I
        tspan=(time[1],time[end])
        prob=ODEProblem(h,u,tspan,phi)
        sol=solve(prob,Vern9(),saveat=time)
        estimated=reduce(hcat,sol.u)
        return estimated[1,:]
    end
    lssk(ti,pi)=lss(idata,ti,pi)
    fir5=curve_fit(lssk,tdata,ndata,p4)
    K[i,1]=fir5.param[1]
    K[i,2]=confidence_interval(fir5,0.05)[1][1]
    K[i,3]=confidence_interval(fir5,0.05)[1][2]
    K[i,4]=dot(fir5.resid,fir5.resid)
end

#Writing N's
GrandMatrix=zeros(1000,7)
GrandMatrix[1:1000,1]=N[1:1000,3]
GrandMatrix[1:1000,2]=M[1:1000,3]
GrandMatrix[:,3]=S[:,1]
GrandMatrix[:,4]=L[:,1]
GrandMatrix[:,5]=K[:,1]
GrandMatrix[:,6]=T[:,1]
GrandMatrix[:,7]=R
CSV.write("complex_simulations_results\\complex_sim_n_9_1.csv",  DataFrame(GrandMatrix), writeheader=false)

#N CI lower
LowerNMatrix=zeros(1000,7)
LowerNMatrix[:,1]=N[:,8]
LowerNMatrix[:,2]=M[:,8]
LowerNMatrix[:,3]=S[:,3]
LowerNMatrix[:,4]=L[:,3]
LowerNMatrix[:,5]=K[:,2]
CSV.write("complex_simulations_results\\complex_sim_n_lower_9_1.csv",  DataFrame(LowerNMatrix), writeheader=false)
 #N CI upper
UpperNMatrix=zeros(1000,7)
UpperNMatrix[:,1]=N[:,9]
UpperNMatrix[:,2]=M[:,9]
UpperNMatrix[:,3]=S[:,4]
UpperNMatrix[:,4]=L[:,4]
UpperNMatrix[:,5]=K[:,3]
CSV.write("complex_simulations_results\\complex_sim_n_upper_9_1.csv",  DataFrame(UpperNMatrix), writeheader=false)


#Writing Betas, CIs
BetaMatrix=zeros(1000,9)
BetaMatrix[:,1]=N[:,4]
BetaMatrix[:,2]=N[:,1]
BetaMatrix[:,3]=N[:,5]
BetaMatrix[:,4]=M[:,4]
BetaMatrix[:,5]=M[:,1]
BetaMatrix[:,6]=M[:,5]
BetaMatrix[:,7]=L[:,5]
BetaMatrix[:,8]=L[:,2]
BetaMatrix[:,9]=L[:,6]
CSV.write("complex_simulations_results\\complex_sim_beta_9_1.csv",  DataFrame(BetaMatrix), writeheader=false)

#Writing Gammas, CIs
GammaMatrix=zeros(1000,9)
GammaMatrix[:,1]=N[:,6]
GammaMatrix[:,2]=N[:,2]
GammaMatrix[:,3]=N[:,7]
GammaMatrix[:,4]=M[:,6]
GammaMatrix[:,5]=M[:,2]
GammaMatrix[:,6]=M[:,7]
GammaMatrix[:,7]=S[:,5]
GammaMatrix[:,8]=S[:,2]
GammaMatrix[:,9]=S[:,6]
CSV.write("complex_simulations_results\\complex_sim_gamma_9_1.csv",  DataFrame(GammaMatrix), writeheader=false)

#Computing SSE
SSE=zeros(1000,7)
for j in 1:1000
datas=dim[j,:]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,[newdata[1],0],tspan,N[j,1:3]) #1 fitting to I
sol1=solve(prob1,Vern9(),saveat=t)
prob2=ODEProblem(f,[newdata[1],0],tspan,[T[j,2]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #2 beta via r0
sol2=solve(prob2,Vern9(),saveat=t)
prob3=ODEProblem(f,[newdata[1],0],tspan,[L[j,2],p1[2],L[j,1]]) #3 fixed gamma
sol3=solve(prob3,Vern9(),saveat=t)
prob4=ODEProblem(f,[newdata[1],0],tspan,M[j,1:3]) #4 fitting to r
sol4=solve(prob4,Vern9(),saveat=t)
prob5=ODEProblem(f,[newdata[1],0],tspan,[T[j,2]*p1[2]/K[j,1],p1[2],K[j,1]]) #combination
sol5=solve(prob5,Vern9(),saveat=t)
prob6=ODEProblem(f,[newdata[1],0],tspan,[p1[1],p1[2],T[j,1]]) #final size
sol6=solve(prob6,Vern9(),saveat=t)
prob7=ODEProblem(f,[newdata[1],0],tspan,[p1[1],p1[2],R[j,1]]) #imax
sol7=solve(prob7,Vern9(),saveat=t)
SSE[j,1]=sum((sol1[1,:]-newdata).^2)/L[j,8]
SSE[j,3]=sum((sol2[1,:]-newdata).^2)/L[j,8]
SSE[j,4]=sum((sol3[1,:]-newdata).^2)/L[j,8]
SSE[j,2]=sum((sol4[1,:]-newdata).^2)/L[j,8]
SSE[j,5]=sum((sol5[1,:]-newdata).^2)/L[j,8]
SSE[j,6]=sum((sol6[1,:]-newdata).^2)/L[j,8]
SSE[j,7]=sum((sol7[1,:]-newdata).^2)/L[j,8]
end
CSV.write("complex_simulations_results\\complex_sim_sse_3_17.csv",  DataFrame(SSE), writeheader=false)



# Plotting data and solution for specific data
j=111
datas=dim[j,:]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,[newdata[1],0],tspan,N[j,1:3]) #1 3p to i
sol1=solve(prob1,Vern9(),saveat=t)
prob2=ODEProblem(f,[newdata[1],0],tspan,[T[j,2]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #2 2p with R0
sol2=solve(prob2,Vern9(),saveat=t)
prob3=ODEProblem(f,[newdata[1],0],tspan,[L[j,2],0.2,L[j,1]]) #3 2p with gamma=0.2
sol3=solve(prob3,Vern9(),saveat=t)
prob4=ODEProblem(f,[newdata[1],0],tspan,M[j,1:3]) #4 3p to r
sol4=solve(prob4,Vern9(),saveat=t)
prob5=ODEProblem(f,[newdata[1],0],tspan,[T[j,2]*0.2/K[j,1],0.2,K[j,1]]) 
sol5=solve(prob5,Vern9(),saveat=t)
prob6=ODEProblem(f,[newdata[1],0],tspan,[T[j,2]*0.2/T[j,1],0.2,T[j,1]]) 
sol6=solve(prob6,Vern9(),saveat=t)
prob7=ODEProblem(f,[newdata[1],0],tspan,[T[j,2]*0.2/R[j,1],0.2,R[j,1]]) 
sol7=solve(prob7,Vern9(),saveat=t)
dataplot=scatter(t,datas,labels="Data",ylabel="Number of currently infected",size=(700,500))
plot1=plot!(dataplot,sol1,vars=(0,1),xlims=[t[1],t[end]],labels="Fitting to I",legend=:outerright,linewidth=3)
plot2=plot!(plot1,sol2,vars=(0,1),xlims=[t[1],t[end]],labels="Beta via R0",linewidth=3)
plot3=plot!(plot2,sol3,vars=(0,1),xlims=[t[1],t[end]],labels="Fixed gamma",linewidth=3)
plot4=plot!(plot3,sol4,vars=(0,1),xlims=[t[1],t[end]],labels="Fitting to R",linewidth=3)
plot5=plot!(plot3,sol5,vars=(0,1),xlims=[t[1],t[end]],labels="Combination",linewidth=3)
plot6=plot!(plot5,sol6,vars=(0,1),xlims=[t[1],t[end]],labels="Final size",linewidth=3)
plot7=plot!(plot6,sol7,vars=(0,1),xlims=[t[1],t[end]],xlabel="Days since introduction",
labels="Imax",legendfontsize=15,linewidth=3,yguidefontsize=15,xguidefontsize=15,legend=:outerright)
allplots=plot(plot1,plot2,plot3,plot4,plot5,plot6,plot7, layout=(3,3),legend=false,size=(1200,900), titlefontsize=10)
savefig(plot7,"complex_simulations_results\\example_9_1.png")