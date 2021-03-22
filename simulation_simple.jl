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

#simulated data, change it!
di=DataFrame(CSV.File("simulations\\infec_1000_25.csv",header=false))
dr=DataFrame(CSV.File("simulations\\remov_1000_25.csv",header=false))

dim=convert(Matrix,di)
drm=convert(Matrix, dr)
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

#Fitting on 3 parameters to . Note that for following methods syntax is similar
function lss(u,time, phi) #IR extract of I
    tspan=(time[1],time[end])
    prob=ODEProblem(f,u,tspan,phi)
    sol=solve(prob,Vern9(),saveat=time)
    estimated=reduce(hcat,sol.u)
    return estimated[1,:]
end

p1=[0.0005,0.2,1000.0]
N=zeros(1000,8) #The main matrix
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
    N[i,4]=findmax(ndata)[1] #Imax
    N[i,5]=findmax(ndata)[2] #index of Imax
    N[i,6]=fir.param[1]*fir.param[3]/fir.param[2] #r0 based on the fitted
    N[i,7]=rnot #r0 via growth rate
    N[i,8]=dot(fir.resid,fir.resid) #residual squared
end


#Fitting on 3 parameters to R
M=zeros(1000,3)
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
end

#final size
T=zeros(1000,1)

for i in 1:1000
    data=dim[i,:] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dtm[i,:] 
    ir10=irdata[ind[1]+9]
    r=log(ir10)/10 
    rnot=1+6*r
    si(N)=rnot-N/irdata[ind[2]]*log((N-1)/(-irdata[ind[2]]))
    T[i]=find_zero(si,(irdata[ind[2]]+1,irdata[ind[2]]+10000),Bisection())  
end

#beta via rnot
S=zeros(1000,3)
p2=[p1[3],p1[2]]
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
    S[i,3]=dot(fir3.resid,fir3.resid)
end

#gamma=0.2
L=zeros(1000,3)
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
    L[i,3]=dot(fir4.resid,fir4.resid)
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

#beta via rnot, gamma fixed
K=zeros(1000,2)
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
    K[i,2]=dot(fir5.resid,fir5.resid)
end

#Writing N's
GrandMatrix=zeros(1000,7)
GrandMatrix[1:1000,1]=N[1:1000,3]
GrandMatrix[1:1000,2]=M[1:1000,3]
GrandMatrix[:,3]=S[:,1]
GrandMatrix[:,4]=L[:,1]
GrandMatrix[:,5]=K[:,1]
GrandMatrix[:,6]=T
GrandMatrix[:,7]=R
CSV.write("simulations\\sim_fitting_25_1000_cut_40.csv",  DataFrame(GrandMatrix), writeheader=false)

#Writing Betas
BetaMatrix=zeros(1000,5)
BetaMatrix[:,1]=N[:,1]
BetaMatrix[:,2]=M[:,1]
BetaMatrix[:,3]=N[:,7].*S[:,2]./S[:,1]
BetaMatrix[:,4]=L[:,2]
BetaMatrix[:,5]=N[:,7].*p1[2]./K[:,1]
CSV.write("simulations\\beta_fitting_25_1000_cut_40.csv",  DataFrame(BetaMatrix), writeheader=false)

#Writing Gammas
GammaMatrix=zeros(1000,3)
GammaMatrix[:,1]=N[:,2]
GammaMatrix[:,2]=M[:,2]
GammaMatrix[:,3]=S[:,2]
CSV.write("simulations\\gamma_fitting_25_1000_cut_40.csv",  DataFrame(GammaMatrix), writeheader=false)

#Computing SSE
SSE=zeros(1000,7)
for j in 1:1000
datas=dim[j,:]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,[newdata[1],0],tspan,N[j,1:3]) #1 3p to i
sol1=solve(prob1,Vern9(),saveat=t)
prob2=ODEProblem(f,[newdata[1],0],tspan,[N[j,7]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #2 2p with R0
sol2=solve(prob2,Vern9(),saveat=t)
prob3=ODEProblem(f,[newdata[1],0],tspan,[L[j,2],p1[2],L[j,1]]) #3 2p with gamma=0.2
sol3=solve(prob3,Vern9(),saveat=t)
prob4=ODEProblem(f,[newdata[1],0],tspan,M[j,1:3]) #4 3p to r
sol4=solve(prob4,Vern9(),saveat=t)
prob5=ODEProblem(f,[newdata[1],0],tspan,[N[j,7]*p1[2]/K[j,1],p1[2],K[j,1]]) 
sol5=solve(prob5,Vern9(),saveat=t)
prob6=ODEProblem(f,[newdata[1],0],tspan,[p1[1],p1[2],T[j,1]]) 
sol6=solve(prob6,Vern9(),saveat=t)
prob7=ODEProblem(f,[newdata[1],0],tspan,[p1[1],p1[2],R[j,1]]) 
sol7=solve(prob7,Vern9(),saveat=t)
SSE[j,1]=sum((sol1[1,:]-newdata).^2)
SSE[j,3]=sum((sol2[1,:]-newdata).^2)
SSE[j,4]=sum((sol3[1,:]-newdata).^2)
SSE[j,2]=sum((sol4[1,:]-newdata).^2)
SSE[j,5]=sum((sol5[1,:]-newdata).^2)
SSE[j,6]=sum((sol6[1,:]-newdata).^2)
SSE[j,7]=sum((sol7[1,:]-newdata).^2)
end
CSV.write("simulations\\SSE_fitting_25_1000g.csv",  DataFrame(SSE), writeheader=false)



# Plotting data and solution for specific
j=1
datas=dim[j,:]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,[newdata[1],0],tspan,N[j,1:3]) #1 3p to i
sol1=solve(prob1,Vern9(),saveat=t)
prob2=ODEProblem(f,[newdata[1],0],tspan,[N[j,7]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #2 2p with R0
sol2=solve(prob2,Vern9(),saveat=t)
prob3=ODEProblem(f,[newdata[1],0],tspan,[L[j,2],p1[2],L[j,1]]) #3 2p with gamma=0.2
sol3=solve(prob3,Vern9(),saveat=t)
prob4=ODEProblem(f,[newdata[1],0],tspan,M[j,1:3]) #4 3p to r
sol4=solve(prob4,Vern9(),saveat=t)
prob5=ODEProblem(f,[newdata[1],0],tspan,[N[j,7]*p1[2]/K[j,1],p1[2],K[j,1]]) 
sol5=solve(prob5,Vern9(),saveat=t)
prob6=ODEProblem(f,[newdata[1],0],tspan,[p1[1],p1[2],T[j,1]]) 
sol6=solve(prob6,Vern9(),saveat=t)
prob7=ODEProblem(f,[newdata[1],0],tspan,[p1[1],p1[2],R[j,1]]) 
sol7=solve(prob7,Vern9(),saveat=t)
dataplot=scatter(t,datas)
plot1=plot!(dataplot,sol1,vars=(0,1),xlims=[t[1],t[end]],title="3 parameters, I data",legend=false)
plot2=plot!(dataplot,sol2,vars=(0,1),xlims=[t[1],t[end]],title="3 parameters, R data",legend=false)
plot3=plot!(dataplot,sol3,vars=(0,1),xlims=[t[1],t[end]],title="2 parameters, beta via R0",legend=false)
plot4=plot!(dataplot,sol4,vars=(0,1),xlims=[t[1],t[end]],title="2 parameters, Gamma=0.2",legend=false)
plot5=plot!(dataplot,sol5,vars=(0,1),xlims=[t[1],t[end]],title="1 parameter",legend=false)
plot6=plot!(dataplot,sol6,vars=(0,1),xlims=[t[1],t[end]],title="Final size",legend=false)
plot7=plot!(dataplot,sol7,vars=(0,1),xlims=[t[1],t[end]],title="Imax",legend=false)
allplots=plot(plot1,plot2,plot3,plot4,plot5,plot6,plot7, layout=(3,3),legend=false,size=(1200,900), titlefontsize=10)
savefig(allplots,"simulations\\allplots.png")

