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
df=DataFrame(CSV.File("City_Confirmed_0115_0816_infected.csv")) #I-data
df1=convert(Matrix,df)

dg=DataFrame(CSV.File("City_Confirmed_0115_0816_total.csv")) #I+R-data
dg1=convert(Matrix,dg)

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


ld=0 #extension of cutted data
p1=[0.0005,0.15,1000.0] #starting point for parameters
N=zeros(53,11) #Fitting to I method
for i in 1:53
    data=df1[i,3:200] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]+ld] #new cutted data
    irdata=dg1[i,3:200] #total cases data
    r=log(irdata[ind[1]+9])/10 #growth rate
    rnot=1+6*r #r0
    idata=[ndata[1],0.0] #initial point for the system
    tdata=collect(1.0:float(length(ndata))) #time array
    lssk(ti,pi)=lss(idata,ti,pi)  #new function with fixed i.p.
    fir=curve_fit(lssk,tdata,ndata,p1) #fitting
    N[i,1:3]=fir.param #saving beta gamma N #residual squared
    N[i,4]=confidence_interval(fir,0.05)[1][1] #CI's
    N[i,5]=confidence_interval(fir,0.05)[1][2]
    N[i,6]=confidence_interval(fir,0.05)[2][1]
    N[i,7]=confidence_interval(fir,0.05)[2][2]
    N[i,8]=confidence_interval(fir,0.05)[3][1]
    N[i,9]=confidence_interval(fir,0.05)[3][2]
    N[i,10]=rnot #r0 via growth rate
    N[i,11]=dot(fir.resid,fir.resid)
end
#Other methods follows similar structure


#Fitting on 3 parameters to R
M=zeros(53,10)
p2=[0.0005,0.15,1000.0]
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data)
    rdata=dg1[i,2+ind[1]:2+ind[2]+ld]-df1[i,2+ind[1]:2+ind[2]+ld]
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
    fir=curve_fit(lssrk,t,rdata,p2)
    M[i,1:3]=fir.param
    M[i,4]=confidence_interval(fir,0.05)[1][1]
    M[i,5]=confidence_interval(fir,0.05)[1][2]
    M[i,6]=confidence_interval(fir,0.05)[2][1]
    M[i,7]=confidence_interval(fir,0.05)[2][2]
    M[i,8]=confidence_interval(fir,0.05)[3][1]
    M[i,9]=confidence_interval(fir,0.05)[3][2]
    M[i,10]=dot(fir.resid,fir.resid)
end


#Final size
T=zeros(53,2)
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dg1[i,3:200] 
    ir10=irdata[ind[1]+9]
    r=log(ir10)/10 
    rnot=1+6*r
    si(N)=rnot-N/irdata[ind[2]]*log((N-ndata[1])/(N-irdata[ind[2]]))
    T[i,1]=find_zero(si,(irdata[ind[2]]+1,irdata[ind[2]]+10000),Bisection()) 
    T[i,2]=irdata[ind[2]]
end

CSV.write("final_sizes_china.csv", DataFrame(T),writeheader=false)

#Imax 
R=zeros(53,2)
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]] 
    irdata=dg1[i,3:200] 
    ir10=irdata[ind[1]+9]
    r=log(ir10)/10 
    rnot=1+6*r
    si(N)=findmax(ndata)[1]+N/rnot-N/rnot*log(N/rnot)-N+N/rnot*log(N-ndata[1])
    R[i,1]=find_zero(si,(irdata[ind[2]]+1,irdata[ind[2]]+10000),Bisection())   
    R[i,2]=findmax(ndata)[1]
end
CSV.write("imaxs_china.csv", DataFrame(R),writeheader=false)
#beta via rnot
S=zeros(53,7)
p3=[1000.0,0.15]
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]+ld] 
    irdata=dg1[i,3:200] 
    r=log(irdata[ind[1]+9])/10 
    rnot=1+6*r
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #we have to create separate functions
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
    fir=curve_fit(lssk,tdata,ndata,p3)
    S[i,1:2]=fir.param
    S[i,3]=confidence_interval(fir,0.05)[1][1]
    S[i,4]=confidence_interval(fir,0.05)[1][2]
    S[i,5]=confidence_interval(fir,0.05)[2][1]
    S[i,6]=confidence_interval(fir,0.05)[2][2]
    S[i,7]=dot(fir.resid,fir.resid)
end


#gamma=0.2
L=zeros(53,7)
p4=[1000.0,0.0005]
for i in 1:53
    data=df1[i,3:200] 
    ind=sir(data) 
    ndata=data[ind[1]:ind[2]+ld] 
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #this function is different from beta via r0
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
end

#beta via rnot, gamma fixed
K=zeros(53,4)
p5=[1000.0]
for i in 1:53
    data=df1[i,3:200] # read data
    ind=sir(data) #cut the data
    ndata=data[ind[1]:ind[2]+ld] #new cutted data
    irdata=dg1[i,3:200] #total cases data
    r=log(irdata[ind[1]+9])/10 #growth rate
    rnot=1+6*r
    idata=[ndata[1],0.0]
    tdata=collect(1.0:float(length(ndata)))
    function h(du,u,p,t) #IR We have to specify for each city separately
            du[1]=dx=rnot*1/6*(p[1]-u[1]-u[2])*u[1]/p[1]-1/6*u[1]
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

#Saving all N results
GrandMatrix=zeros(53,7)
GrandMatrix[1:53,1]=N[1:53,3] #ifit
GrandMatrix[1:53,2]=M[1:53,3] #rfit
GrandMatrix[1:53,3]=S[1:53,1] #beta via r0
GrandMatrix[1:53,4]=L[1:53,1] #fixed gamma
GrandMatrix[1:53,5]=K[1:53,1] #combination
GrandMatrix[1:53,6]=R[1:53] #imax
GrandMatrix[1:53,7]=T[1:53] #final size

CSV.write("China_Ns.csv",  DataFrame(GrandMatrix), writeheader=false)

#N lower CI
LowerNMatrix=zeros(53,7)
LowerNMatrix[1:53,1]=N[1:53,8]
LowerNMatrix[1:53,2]=M[1:53,8]
LowerNMatrix[1:53,3]=S[1:53,3]
LowerNMatrix[1:53,4]=L[1:53,3]
LowerNMatrix[1:53,5]=K[1:53,2]
CSV.write("China_Ns_lower.csv",  DataFrame(LowerNMatrix), writeheader=false)

#N upper CI
UpperNMatrix=zeros(53,7)
UpperNMatrix[1:53,1]=N[1:53,9]
UpperNMatrix[1:53,2]=M[1:53,9]
UpperNMatrix[1:53,3]=S[1:53,4]
UpperNMatrix[1:53,4]=L[1:53,4]
UpperNMatrix[1:53,5]=K[1:53,3]
CSV.write("China_Ns_upper.csv",  DataFrame(UpperNMatrix), writeheader=false)

#Betas with CI
BetaMatrix=zeros(53,9)
BetaMatrix[:,1]=N[:,4]
BetaMatrix[:,2]=N[:,1]
BetaMatrix[:,3]=N[:,5]
BetaMatrix[:,4]=M[:,4]
BetaMatrix[:,5]=M[:,1]
BetaMatrix[:,6]=M[:,5]
BetaMatrix[:,7]=L[:,5]
BetaMatrix[:,8]=L[:,2]
BetaMatrix[:,9]=L[:,6]
CSV.write("China_Betas.csv", DataFrame(BetaMatrix),writeheader=false)

#Gammas with CI
GammaMatrix=zeros(53,9)
GammaMatrix[:,1]=N[:,6]
GammaMatrix[:,2]=N[:,2]
GammaMatrix[:,3]=N[:,7]
GammaMatrix[:,4]=M[:,6]
GammaMatrix[:,5]=M[:,2]
GammaMatrix[:,6]=M[:,7]
GammaMatrix[:,7]=S[:,5]
GammaMatrix[:,8]=S[:,2]
GammaMatrix[:,9]=S[:,6]
CSV.write("China_Gammas.csv", DataFrame(GammaMatrix),writeheader=false)

#Computing SSE, transcedental are excluded, but are similar to combination
SSE=zeros(53,5)
for j in 1:53
datas=df1[j,3:200]
totaldatas=dg1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
rdata=float.(totaldatas[sir(datas)[1]:sir(datas)[2]])-newdata
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,[newdata[1],0],tspan,N[j,1:3]) #1 fitting to i
sol1=solve(prob1,Vern9(),saveat=t)
prob2=ODEProblem(f,[newdata[1],0],tspan,[N[j,10]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #2 beta via r0
sol2=solve(prob2,Vern9(),saveat=t)
prob3=ODEProblem(f,[newdata[1],0],tspan,[L[j,2],1/6,L[j,1]]) #3 fixed gamma
sol3=solve(prob3,Vern9(),saveat=t)
prob4=ODEProblem(f,[newdata[1],0],tspan,M[j,1:3]) #4 fitting to r
sol4=solve(prob4,Vern9(),saveat=t)
prob5=ODEProblem(f,[newdata[1],0],tspan,[N[j,10]/K[j,1]/6,1/6,K[j,1]])  #combination
sol5=solve(prob5,Vern9(),saveat=t)
SSE[j,1]=log(sum((sol1[1,:]-newdata).^2))
SSE[j,3]=log(sum((sol2[1,:]-newdata).^2))
SSE[j,4]=log(sum((sol3[1,:]-newdata).^2))
SSE[j,2]=log(sum((sol4[1,:]-newdata).^2))
SSE[j,5]=log(sum((sol5[1,:]-newdata).^2))
end
CSV.write("China_SSE.csv",  DataFrame(SSE), writeheader=false)

#Plotting
allplots=Matrix{Any}(nothing, 53,8)
for j in 1:53
datas=df1[j,3:200]
totaldatas=dg1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]+ld])
rdata=float.(totaldatas[sir(datas)[1]:sir(datas)[2]+ld])-newdata
u0=[newdata[1],0.0]
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,u0,tspan,N[j,1:3]) #ifit 1
sol1=solve(prob1,Vern9(),saveat=t)
plot1=plot(sol1,vars=(0,1),xlims=[t[1],t[end]],legend=false,ylabel="I, the number of infected", xlabel="t, days")
plot12=scatter!(plot1,t,newdata)
prob2=ODEProblem(f,u0,tspan,M[j,1:3]) #rfit 3
sol2=solve(prob2,Vern9(),saveat=t)
plot2=plot(sol2,vars=(0,1),xlims=[t[1],t[end]],legend=false,ylabel="I, the number of infected", xlabel="t, days")
plot22=scatter!(plot2,t,newdata)
prob4=ODEProblem(f,u0,tspan,[N[j,10]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #rnot 5
sol4=solve(prob4,Vern9(),saveat=t)
plot4=plot(sol4,vars=(0,1),xlims=[t[1],t[end]],legend=false,ylabel="I, the number of infected", xlabel="t, days")
plot42=scatter!(plot4,t,newdata)
prob5=ODEProblem(f,u0,tspan,[L[j,2],1/6,L[j,1]]) #fixed gamma 2
sol5=solve(prob5,Vern9(),saveat=t)
plot5=plot(sol5,vars=(0,1),xlims=[t[1],t[end]],legend=false,ylabel="I, the number of infected", xlabel="t, days")
plot52=scatter!(plot5,t,newdata)
prob6=ODEProblem(f,u0,tspan,[N[j,10]/K[j,1]/6,1/6,K[j,1]]) #combo 6
sol6=solve(prob6,Vern9(),saveat=t)
plot6=plot(sol6,vars=(0,1),xlims=[t[1],t[end]],legend=false,ylabel="I, the number of infected", xlabel="t, days")
plot62=scatter!(plot6,t,newdata)
prob7=ODEProblem(f,u0,tspan,M[j,1:3]) #rfit r 4
sol7=solve(prob7,Vern9(),saveat=t)
plot7=plot(sol7,vars=(0,2),xlims=[t[1],t[end]],legend=false,ylabel="R, the number of removed", xlabel="t, days")
plot72=scatter!(plot7,t,rdata)
prob8=ODEProblem(f,u0,tspan,[N[j,10]/T[j,1]/6,1/6,T[j,1]]) #7 fin size
sol8=solve(prob8,Vern9(),saveat=t)
plot8=plot(sol8,vars=(0,1),xlims=[t[1],t[end]],legend=false,ylabel="I, the number of infected", xlabel="t, days")
plot82=scatter!(plot8,t,newdata)
prob9=ODEProblem(f,u0,tspan,[N[j,10]/R[j,1]/6,1/6,R[j,1]]) #8 imax
sol9=solve(prob9,Vern9(),saveat=t)
plot9=plot(sol9,vars=(0,1),xlims=[t[1],t[end]],legend=false,ylabel="I, the number of infected", xlabel="t, days")
plot92=scatter!(plot9,t,newdata)
allplots[j,:]=[plot12,plot52,plot22,plot72,plot42,plot62,plot82,plot92]
end
#Display all of them in separate windows
k=12
savefig(allplots[k,1],"chinamaps/data_sol_12_1.png")
savefig(allplots[k,2],"chinamaps/data_sol_12_2.png")
savefig(allplots[k,3],"chinamaps/data_sol_12_3.png")
savefig(allplots[k,4],"chinamaps/data_sol_12_4.png")
savefig(allplots[k,5],"chinamaps/data_sol_12_5.png")
savefig(allplots[k,6],"chinamaps/data_sol_12_6.png")
savefig(allplots[k,7],"chinamaps/data_sol_12_7.png")
savefig(allplots[k,8],"chinamaps/data_sol_12_8.png")

panelplots = Array{Any}(nothing, 53)
for i in 1:53
panelplots[i]=plot(allplots[i,1],allplots[i,2],allplots[i,3],allplots[i,4],allplots[i,5],allplots[i,6],allplots[i,7], allplots[i,8],layout=(4,2),size=(800,1200))
end

for i in 1:53
display(panelplots[i])
end

#Computing SSE and AICc for 3 methods that had best fitting
SSE=zeros(53,3)
AIC=zeros(53,3)
for j in 1:53
datas=df1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
t=collect(1.0:float(length(newdata)))
tspan=(t[1],t[end])
prob1=ODEProblem(f,[newdata[1],0],tspan,N[j,1:3]) #1 3p
sol1=solve(prob1,Vern9(),saveat=t)
prob2=ODEProblem(f,[newdata[1],0],tspan,[N[j,10]*S[j,2]/S[j,1],S[j,2],S[j,1]]) #2 2p with R0
sol2=solve(prob2,Vern9(),saveat=t)
prob3=ODEProblem(f,[newdata[1],0],tspan,[L[j,2],1/6,L[j,1]]) #2 2p with gamma=0.2
sol3=solve(prob3,Vern9(),saveat=t)
SSE[j,1]=sum((sol1[1,:]-newdata).^2)
SSE[j,2]=sum((sol2[1,:]-newdata).^2)
SSE[j,3]=sum((sol3[1,:]-newdata).^2)
AIC[j,1]=length(t)*log(SSE[j,1]/length(t))+8+40/(length(t)-5)
AIC[j,2]=length(t)*log(SSE[j,2]/length(t))+6+24/(length(t)-4)
AIC[j,3]=length(t)*log(SSE[j,3]/length(t))+6+24/(length(t)-4)
end
CSV.write("China_AICc.csv",  DataFrame(AIC), writeheader=false)

#Heatmap function for SSE, beta via R0
function SSEheatmap(a,b,data,j)
    t=collect(1.0:float(length(data)))
    tspan=(t[1],t[end])
    prob1=ODEProblem(f,[data[1],0],tspan,[N[j,11]*a/b,a,b])
    sol1=solve(prob1,Vern9(),saveat=t)
    return log(sum((sol1[1,:]-data).^2))
end

SSEgamma=Array{Any}(nothing, 53) #Array for the gamma
for i in 1:53 
    SSEgamma[i]=collect(S[i,2]*0.1:S[i,2]*9.9/200:S[i,2]*10)
end

SSEN=Array{Any}(nothing,53) #Array for N
for i in 1:53 
    SSEN[i]=collect(S[i,1]*0.1:S[i,1]*9.9/200:S[i,1]*10)
end

SSEheatmapplots = Array{Any}(nothing, 53)
for j in 1:53
datas=df1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
SSEheatmaptest(a,b)=SSEheatmap(a,b,newdata,j)
plot1=heatmap(SSEgamma[j],SSEN[j],SSEheatmaptest)
plot2=contour!(plot1,SSEgamma[j],SSEN[j],SSEheatmaptest,levels=[log(SSE[j,2])*1.01],color=:red,linewidth=5)
plot3=scatter!(plot2,(S[j,2],S[j,1]),color=:blue,legend=false)
SSEheatmapplots[j]=plot3
end

for i in 1:53
    display(SSEheatmapplots[i])
    end

#Heatmap function for SSE, fitting to I
function SSEheatmap(a,b,data,j)
    t=collect(1.0:float(length(data)))
    tspan=(t[1],t[end])
    prob1=ODEProblem(f,[data[1],0],tspan,[N[j,1],a,b])
    sol1=solve(prob1,Vern9(),saveat=t)
    return log(sum((sol1[1,:]-data).^2))
end

SSEgamma=Array{Any}(nothing, 53) #Array for the gamma
for i in 1:53 
    SSEgamma[i]=collect(N[i,2]*0.1:N[i,2]*9.9/50:N[i,2]*10)
end

SSEN=Array{Any}(nothing,53) #Array for N
for i in 1:53
    SSEN[i]=collect(N[i,3]*0.1:N[i,3]*9.9/50:N[i,3]*10)
end

SSEheatmapplots = Array{Any}(nothing, 53)
for j in 1:53
datas=df1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
SSEheatmaptest(a,b)=SSEheatmap(a,b,newdata,j)
plot1=heatmap(SSEgamma[j],SSEN[j],SSEheatmaptest)
plot2=contour!(plot1,SSEgamma[j],SSEN[j],SSEheatmaptest,levels=[(SSE[j,1])*1.01],color=:red,linewidth=5)
SSEheatmapplots[j]=plot2
end

for i in 1:53
    display(SSEheatmapplots[i])
    end

#SSE fitting to R
function SSEheatmap(a,b,data,st,j)
    t=collect(1.0:float(length(data)))
    tspan=(t[1],t[end])
    prob1=ODEProblem(f,[st,0],tspan,[M[j,1],a,b])
    sol1=solve(prob1,Vern9(),saveat=t)
    return log(sum((sol1[2,:]-data).^2))
end

SSEgamma=Array{Any}(nothing, 53) #Array for the gamma
for i in 1:53 
    SSEgamma[i]=collect(M[i,2]*0.1:M[i,2]*4.9/200:M[i,2]*5)
end

SSEN=Array{Any}(nothing,53) #Array for N
for i in 1:53 
    SSEN[i]=collect(M[i,3]*0.1:M[i,3]*4.9/200:M[i,3]*5)
end

SSEheatmapplots = Array{Any}(nothing, 53)
for j in 1:53
datas=df1[j,3:200]
totaldatas=dg1[j,3:200]
newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
rdata=float.(totaldatas[sir(datas)[1]:sir(datas)[2]])-newdata
SSEheatmaptest(a,b)=SSEheatmap(a,b,rdata,newdata[1],j)
plot1=heatmap(SSEgamma[j],SSEN[j],SSEheatmaptest)
plot2=contour!(plot1,SSEgamma[j],SSEN[j],SSEheatmaptest,levels=[log(M[j,10])*1.01],color=:red,linewidth=5)
plot3=scatter!(plot2,(M[j,2],M[j,3]),color=:blue,legend=false)
SSEheatmapplots[j]=plot3
end

for i in 1:53
    display(SSEheatmapplots[i])
end


#sse, fixed gamma
    function SSEheatmap(a,b,data,j)
        t=collect(1.0:float(length(data)))
        tspan=(t[1],t[end])
        prob1=ODEProblem(f,[data[1],0],tspan,[a,1/6,b])
        sol1=solve(prob1,Vern9(),saveat=t)
        return log(sum((sol1[1,:]-data).^2))
    end
    
    SSEgamma=Array{Any}(nothing, 53) #Array for the gamma
    for i in 1:53 
        SSEgamma[i]=collect(L[i,2]*0.1:L[i,2]*4.9/200:L[i,2]*5)
    end
    
    SSEN=Array{Any}(nothing,53) #Array for N
    for i in 1:53 
        SSEN[i]=collect(L[i,1]*0.1:L[i,1]*4.9/200:L[i,1]*5)
    end
    
    SSEheatmapplots = Array{Any}(nothing, 53)
    for j in 1:53
    datas=df1[j,3:200]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]+ld])
    SSEheatmaptest(a,b)=SSEheatmap(a,b,newdata,j)
    plot1=heatmap(SSEgamma[j],SSEN[j],SSEheatmaptest)
    plot2=contour!(plot1,SSEgamma[j],SSEN[j],SSEheatmaptest,levels=[log(SSE[j,4])*1.01],color=:red,linewidth=5)
    plot3=scatter!(plot2,(L[j,2],L[j,1]),color=:blue,legend=false)
    SSEheatmapplots[j]=plot3
    end
    
    for i in 1:53
        display(SSEheatmapplots[i])
        end    
