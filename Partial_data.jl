#comments are the same as for main methods

#use original packages
p1=[0.0005,0.2,1000.0]
lock=40 #change accordingly

#Dataframes are the same as in the original

N=zeros(1000,8) #The main matrix
for i in 1:1000
    data=dim[i,:] # read data
    ind=sir(data)
    if ind[1]+lock>ind[2]
        index=ind[2]
    else
        index=ind[1]+lock
    end #cut the data
    ndata=data[ind[1]:index] #new cutted data
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

M=zeros(1000,3)
for i in 1:1000
    data=dim[i,:] 
    ind=sir(data)
    if ind[1]+lock>ind[2]
        index=ind[2]
    else
        index=ind[1]+lock
    end
    rdata=drm[i,ind[1]:index]
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

S=zeros(1000,3)
p2=[p1[3],p1[2]]
for i in 1:1000
    data=dim[i,:] 
    ind=sir(data)
    if ind[1]+lock>ind[2]
        index=ind[2]
    else
        index=ind[1]+lock
    end 
    ndata=data[ind[1]:index] 
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

L=zeros(1000,3)
p3=[p1[3],p1[1]]
for i in 1:1000
    data=dim[i,:] 
    ind=sir(data) 
    if ind[1]+lock>ind[2]
        index=ind[2]
    else
        index=ind[1]+lock
    end
    ndata=data[ind[1]:index] 
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

K=zeros(1000,2)
p4=[p1[3]]
for i in 1:1000
    data=dim[i,:] # read data
    ind=sir(data)
    if ind[1]+lock>ind[2]
        index=ind[2]
    else
        index=ind[1]+lock
    end #cut the data
    ndata=data[ind[1]:index] #new cutted data
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

#Ns
GrandMatrix=zeros(1000,5)
GrandMatrix[1:1000,1]=N[1:1000,3]
GrandMatrix[1:1000,2]=M[1:1000,3]
GrandMatrix[:,3]=S[:,1]
GrandMatrix[:,4]=L[:,1]
GrandMatrix[:,5]=K[:,1]
CSV.write("simulations\\sim_fitting_25_1000_lock_40.csv",  DataFrame(GrandMatrix), writeheader=false)

#Writing Betas
BetaMatrix=zeros(1000,5)
BetaMatrix[:,1]=N[:,1]
BetaMatrix[:,2]=M[:,1]
BetaMatrix[:,3]=N[:,7].*S[:,2]./S[:,1]
BetaMatrix[:,4]=L[:,2]
BetaMatrix[:,5]=N[:,7].*p1[2]./K[:,1]
CSV.write("simulations\\beta_fitting_25_1000_lock_40.csv",  DataFrame(BetaMatrix), writeheader=false)

#Gammas
GammaMatrix=zeros(1000,3)
GammaMatrix[:,1]=N[:,2]
GammaMatrix[:,2]=M[:,2]
GammaMatrix[:,3]=S[:,2]
CSV.write("simulations\\gamma_fitting_25_1000_lock_40.csv",  DataFrame(GammaMatrix), writeheader=false)
