#packages
using Plots
using DifferentialEquations
using LsqFit
using DataFrames
using CSV
using LinearAlgebra
using Optim
using Roots

#csv dataset
df=DataFrame(CSV.File("City_Confirmed_0115_0816.csv"))
df1=convert(Matrix,df)

dg=DataFrame(CSV.File("City_Confirmed_total.csv"))
dg1=convert(Matrix,dg)

population=df1[:,219]
density=df1[:,220]




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
p=[0.0005,0.25,2000.0] #initial guess for parameters

N=zeros(53,10) #The main matrix, it includes ifit method and other parameter estimation
for i in 1:53
    data=df1[i,3:200] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]] #new cutted data
    irdata=dg1[i,3:200] #total cases data
    r=log(irdata[ind[1]+9])/10 #growth rate
    rnot=1+6*r #r0
    idata=[ndata[1],0.0] #initial point for the system
    tdata=collect(1.0:float(length(ndata))) #time array
    lssk(ti,pi)=lss(idata,ti,pi) #new function with fixed i.p.
    fir=curve_fit(lssk,tdata,ndata,p) #fitting
    N[i,1:3]=fir.param #saving beta gamma N
    N[i,4]=findmax(ndata)[1] #Imax
    N[i,5]=findmax(ndata)[2] #index of Imax
    N[i,6]=fir.param[1]*fir.param[3]/fir.param[2] #r0 based on the fitted
    N[i,7]=population[i] #population
    N[i,8]=density[i] #density
    N[i,9]=rnot #r0 via growth rate
    N[i,10]=dot(fir.resid,fir.resid) #residual squared
end

#Fitting on 3 parameters to R
M=zeros(53,3)
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data)
    rdata=dg1[i,2+ind[1]:2+ind[2]]-df1[i,2+ind[1]:2+ind[2]]
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
    fir=curve_fit(lssrk,t,rdata,p)
    M[i,1:3]=fir.param
end

#final size
T=zeros(53)

for i in 1:51
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dg1[i,3:200] 
    ir10=irdata[ind[1]+9]
    r=log(ir10)/10 
    rnot=1+6*r
    si(N)=rnot-N/irdata[ind[2]]*log((N-ndata[1])/(N-irdata[ind[2]]))
    T[i]=find_zero(si,(irdata[ind[2]]+1,irdata[ind[2]]+10000),Bisection())  
end

#Imax 
R=zeros(53)
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dg1[i,3:200] 
    ir10=irdata[ind[1]+9]
    r=log(ir10)/10 
    rnot=1+6*r
    si(N)=findmax(ndata)[1]+N/rnot-N/rnot*log(N/rnot)-N+N/rnot*log(N-ndata[1])
    R[i]=find_zero(si,(irdata[ind[2]]+1,irdata[ind[2]]+10000),Bisection())   
end

#beta via rnot
S=zeros(53,3)
p2=[2000.0,0.2]
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dg1[i,3:200] 
    r=log(irdata[ind[1]+9])/10 
    rnot=1+6*r
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #we have to create seprate functions
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
    fir=curve_fit(lssk,tdata,ndata,p2)
    S[i,1:2]=fir.param
    S[i,3]=dot(fir.resid,fir.resid)
end


#gamma=0.2
L=zeros(53,3)
p3=[2000.0,0.0005]
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #this function is different from beta via r0
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
    fir=curve_fit(lssk,tdata,ndata,p3)
    L[i,1:2]=fir.param
    L[i,3]=dot(fir.resid,fir.resid)
end

#beta via rnot, gamma fixed
K=zeros(53,2)
p4=[2000.0]
for i in 1:53
    data=df1[i,3:200] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]] #new cutted data
    irdata=dg1[i,3:200] #total cases data
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
    fir=curve_fit(lssk,tdata,ndata,p4)
    K[i,1]=fir.param[1]
    K[i,2]=dot(fir.resid,fir.resid)
end

#Saving all N results
GrandMatrix=zeros(53,8)
GrandMatrix[1:53,1]=N[1:53,3] #ifit
GrandMatrix[1:53,2]=M[1:53,3] #rfit
GrandMatrix[1:53,3]=S[1:53,3] #beta via r0
GrandMatrix[1:53,4]=L[1:53,1] #fixed gamma
GrandMatrix[1:53,5]=K[1:53,1] #combination
GrandMatrix[1:53,6]=R[1:53] #imax
GrandMatrix[1:53,7]=T[1:53] #final size

CSV.write("AllMethods.csv",  DataFrame(GrandMatrix), writeheader=false)

#Plotting 3 by 3 fitting each for the each mehtods
allplots = Array{Any}(nothing, 51)
for j in 1:51
datas=df1[j,3:200]
totaldatas=dg1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
rdata=float.(totaldatas[sir(datas)[1]:sir(datas)[2]])-newdata
u0=[newdata[1],0.0]
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,u0,tspan,N[j,1:3])
sol1=solve(prob1,Vern9(),saveat=t)
plot1=plot(sol1,vars=(0,1),xlims=[t[1],t[end]],title="3 parameters, I data", legend=false)
plot12=scatter!(plot1,t,newdata)
prob2=ODEProblem(f,u0,tspan,M[j,1:3])
sol2=solve(prob2,Vern9(),saveat=t)
plot2=plot(sol2,vars=(0,1),xlims=[t[1],t[end]],title="3 parameters, R data",legend=false)
plot22=scatter!(plot2,t,newdata)
prob3=ODEProblem(f,u0,tspan,J[j,1:3])
sol3=solve(prob3,Vern9(),saveat=t)
plot3=plot(sol3,vars=(0,1),xlims=[t[1],t[end]],title="3 parameters, IR data",legend=false)
plot32=scatter!(plot3,t,newdata)
prob4=ODEProblem(f,u0,tspan,[N[j,9]*S[j,2]/S[j,1],S[j,2],S[j,1]])
sol4=solve(prob4,Vern9(),saveat=t)
plot4=plot(sol4,vars=(0,1),xlims=[t[1],t[end]],title="2 parameters, I data, beta via R0",legend=false)
plot42=scatter!(plot4,t,newdata)
prob5=ODEProblem(f,u0,tspan,[L[j,2],0.2,L[j,1]])
sol5=solve(prob5,Vern9(),saveat=t)
plot5=plot(sol5,vars=(0,1),xlims=[t[1],t[end]],title="2 parameters, I data, gamma=0.2",legend=false)
plot52=scatter!(plot5,t,newdata)
prob6=ODEProblem(f,u0,tspan,[N[j,9]*0.2/K[j,1],0.2,K[j,1]])
sol6=solve(prob6,Vern9(),saveat=t)
plot6=plot(sol6,vars=(0,1),xlims=[t[1],t[end]],title="1 parameter, I data, gamma=0.2, beta via R0",legend=false)
plot62=scatter!(plot6,t,newdata)
prob7=ODEProblem(f,u0,tspan,M[j,1:3])
sol7=solve(prob7,Vern9(),saveat=t)
plot7=plot(sol7,vars=(0,2),xlims=[t[1],t[end]],title="3 parameters, R data (points)",legend=false)
plot72=scatter!(plot7,t,rdata)
prob8=ODEProblem(f,u0,tspan,[N[j,9]*0.2/T[j,1],0.2,T[j,1]])
sol8=solve(prob8,Vern9(),saveat=t)
plot8=plot(sol8,vars=(0,1),xlims=[t[1],t[end]],title="Final Size, gamma=0.2, beta via R0",legend=false)
plot82=scatter!(plot8,t,newdata)
prob9=ODEProblem(f,u0,tspan,[N[j,9]*0.2/R[j,1],0.2,R[j,1]])
sol9=solve(prob9,Vern9(),saveat=t)
plot9=plot(sol9,vars=(0,1),xlims=[t[1],t[end]],title="Imax, gamma=0.2, beta via R0",legend=false)
plot92=scatter!(plot9,t,newdata)
allplots[j]=plot(plot12,plot22,plot32,plot42,plot52,plot62,plot72,plot82,plot92, layout=(3,3),legend=false,size=(1200,1200), titlefontsize=10)
end
#Display all of them in separate windows
for i in 1:51
display(allplots[i])
end
df1[37] #Name of the city
savefig(allplots[50],"Dongguan_Full.png") #Save

#Computing SSE and AICc for 3 methods that had best fitting
SSE=zeros(51,3)
AIC=zeros(51,3)
for j in 1:51
datas=df1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,[newdata[1],0],tspan,N[j,1:3]) #1 3p
sol1=solve(prob1,Vern9(),saveat=t)
prob2=ODEProblem(f,[newdata[1],0],tspan,[N[j,9]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #2 2p with R0
sol2=solve(prob2,Vern9(),saveat=t)
prob3=ODEProblem(f,[newdata[1],0],tspan,[L[j,2],0.2,L[j,1]]) #2 2p with gamma=0.2
sol3=solve(prob3,Vern9(),saveat=t)
SSE[j,1]=sum((sol1[1,:]-newdata).^2)
SSE[j,2]=sum((sol2[1,:]-newdata).^2)
SSE[j,3]=sum((sol3[1,:]-newdata).^2)
AIC[j,1]=length(t)*log(SSE[j,1]/length(t))+8+40/(length(t)-5)
AIC[j,2]=length(t)*log(SSE[j,2]/length(t))+6+24/(length(t)-4)
AIC[j,3]=length(t)*log(SSE[j,3]/length(t))+6+24/(length(t)-4)
end

CSV.write("AICc.csv",  DataFrame(AIC), writeheader=false)

#Heatmap function for SSE and the second best method
function SSEheatmap(a,b,data,j)
    t=collect(1.0:float(length(data)))
    tspan=(t[1],t[end])
    prob1=ODEProblem(f,[data[1],0],tspan,[N[j,9]*a/b,a,b])
    sol1=solve(prob1,Vern9(),saveat=t)
    return log(sum((sol1[1,:]-data).^2))
end

SSEgamma=Array{Any}(nothing, 51) #Array for the gamma
for i in 1:51 
    SSEgamma[i]=collect(S[i,2]*0.1:S[i,2]*9.9/200:S[i,2]*10)
end

SSEN=Array{Any}(nothing,51) #Array for N
for i in 1:51 
    SSEN[i]=collect(S[i,1]*0.1:S[i,1]*9.9/200:S[i,1]*10)
end

SSEheatmapplots = Array{Any}(nothing, 51)
for j in 1:51
datas=df1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
SSEheatmaptest(a,b)=SSEheatmap(a,b,newdata,j)
plot1=heatmap(SSEgamma[j],SSEN[j],SSEheatmaptest)
plot2=contour!(plot1,SSEgamma[j],SSEN[j],SSEheatmaptest,levels=[log(SSE[j,2])*1.001],color=:red,linewidth=3)
SSEheatmapplots[j]=plot2
end

for i in 1:51
    display(SSEheatmapplots[i])
    end

#Heatmap function for SSE and the first best method
function SSEheatmap(a,b,data,j)
    t=collect(1.0:float(length(data)))
    tspan=(t[1],t[end])
    prob1=ODEProblem(f,[data[1],0],tspan,[N[j,1],a,b])
    sol1=solve(prob1,Vern9(),saveat=t)
    return log(sum((sol1[1,:]-data).^2))
end

SSEgamma=Array{Any}(nothing, 51) #Array for the gamma
for i in 1:51 
    SSEgamma[i]=collect(N[i,2]*0.1:N[i,2]*9.9/200:N[i,2]*10)
end

SSEN=Array{Any}(nothing,51) #Array for N
for i in 1:51 
    SSEN[i]=collect(N[i,3]*0.1:N[i,3]*9.9/200:N[i,3]*10)
end

SSEheatmapplots = Array{Any}(nothing, 51)
for j in 1:51
datas=df1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
SSEheatmaptest(a,b)=SSEheatmap(a,b,newdata,j)
plot1=heatmap(SSEgamma[j],SSEN[j],SSEheatmaptest)
plot2=contour!(plot1,SSEgamma[j],SSEN[j],SSEheatmaptest,levels=[log(SSE[j,1])*1.001],color=:red,linewidth=3)
SSEheatmapplots[j]=plot2
end

for i in 1:51
    display(SSEheatmapplots[i])
    end


df1[12] #Name of the city
savefig(SSEheatmapplots[12],"Anqing_SSE_big.png") #Save