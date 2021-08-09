#use packages from original file

#True values
truet=collect(1.0:150.0)
truespan=(1.0,150.0)
trueprob=ODEProblem(f,[1,0],truespan,p1) 
truesol=solve(trueprob,Vern9(),saveat=truet)
trueplot=plot(truesol,vars=(0,1),xlims=[0.0,150.0],legend=false,color="black",linewidth=5)

#Fixing optimum parameters (the method)
PopData=N[:,3]
ParamData=N[:,1:3]
ParamData=zeros(1000,3)
ParamData[:,1]=L[:,2]
ParamData[:,2]=[0.1 for i in 1:1000]
ParamData[:,3]=L[:,1]

#Reading the first data
data1=dim[1,:]
newdata1=float.(drm[1,sir(data1)[1]:sir(data1)[2]])
t1=collect(1.0:float(length(newdata1)))
tspan=(t1[1],t1[end])

#Plotting all data
plot8=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],palette=:Dark2_8)
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(drm[j,sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
        plot!(plot8,t,newdata,legend=false,xlims=[0.0,150.0],palette=:Dark2_8)

end
plot8 
plot!(plot8,vars=(0,1),xlims=[0.0,150.0],legend=false,color="black",
linewidth=3, xlabel="Days since introduction",ylabel="Number of removed",size=(600,400))
savefig(plot8,"complex_simulations_results/com_rdata_3_1.png")


# COMPARISON for a single data
#nonvsoutliers
outliers=findall(x->x>percentile(PopData,75)+1.5*iqr(PopData)||x<percentile(PopData,25)-1.5*iqr(PopData),PopData)

data1=dim[1,:]
newdata1=float.(data1[sir(data1)[1]:sir(data1)[2]])
t1=collect(1.0:float(length(newdata1)))
tspan=(t1[1],t1[end])
prob1=ODEProblem(f,[newdata1[1],0],tspan,ParamData[1,:]) 
sol1=solve(prob1,Vern9(),saveat=t1)

#I-solution, N-outliers
if 1 in outliers
    plot7=plot(sol1,vars=(0,1),legend=false,xlims=[0.0,150.0],color="blue",α = 0.2)
else
    plot7=plot(sol1,vars=(0,1),legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    tspan=(t[1],t[end])
    prob=ODEProblem(f,[newdata[1],0],tspan,ParamData[j,:]) 
    sol=solve(prob,Vern9(),saveat=t)
    if j in outliers
        plot!(plot7,sol,vars=(0,1),legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot7,sol,vars=(0,1),legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot7 #to put in frame 1:1
plot!(plot7,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3,title="I-solutions, N-outliers",ylabel="I, the number of infected", xlabel="t, days")

#I-data, N-outliers
if 1 in outliers
    plot8=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot8=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers
        plot!(plot8,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot8,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot8 #1:2
plot!(plot8,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="I-data, N-outliers",ylabel="I, the number of infected", xlabel="t, days")

#R-data, N-outliers
newdata1r=float.(drm[1,sir(data1)[1]:sir(data1)[2]])
if 1 in outliers
    plot9=plot(t1,newdata1r,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot9=plot(t1,newdata1r,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(drm[j,sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers
        plot!(plot9,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot9,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot9 #2:2
plot!(plot9,truesol,vars=(0,2),xlims=[0.0,75.0],legend=false,color="black",linewidth=3,title="R-data, N-outliers",ylabel="R, the number of removed", xlabel="t, days")

#R-solution, N-outliers
if 1 in outliers
    plot10=plot(sol1,vars=(0,2),legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot10=plot(sol1,vars=(0,2),legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    tspan=(t[1],t[end])
    prob=ODEProblem(f,[newdata[1],0],tspan,ParamData[j,:]) 
    sol=solve(prob,Vern9(),saveat=t)
    if j in outliers
        plot!(plot10,sol,vars=(0,2),legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot10,sol,vars=(0,2),legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot10 #to put in frame 2:1
plot!(plot10,truesol,vars=(0,2),xlims=[0.0,75.0],legend=false,color="black",linewidth=3,title="R-solutions, N-outliers",ylabel="R, the number of removed", xlabel="t, days")

#I-data, beta-outliers
outliersB=findall(x->x>percentile(ParamData[:,1],75)+1.5*iqr(ParamData[:,1])||x<percentile(ParamData[:,1],25)-1.5*iqr(ParamData[:,1]),ParamData[:,1])
if 1 in outliersB
    plot11=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot11=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    tspan=(t[1],t[end])
    prob=ODEProblem(f,[newdata[1],0],tspan,ParamData[j,:]) 
    sol=solve(prob,Vern9(),saveat=t)
    if j in outliersB
        plot!(plot11,t, newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot11,t, newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot11 #to put in frame 3:1
plot!(plot11,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3,title="I-data, beta-outliers",ylabel="I, the number of infected", xlabel="t, days")

#I-data, gamma-outliers
outliersG=findall(x->x>percentile(ParamData[:,2],75)+1.5*iqr(ParamData[:,2])||x<percentile(ParamData[:,2],25)-1.5*iqr(ParamData[:,2]),ParamData[:,2])
if 1 in outliersG
    plot12=plot(t1, newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot12=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    tspan=(t[1],t[end])
    prob=ODEProblem(f,[newdata[1],0],tspan,ParamData[j,:]) 
    sol=solve(prob,Vern9(),saveat=t)
    if j in outliersG
        plot!(plot12,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot12,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot12 #to put in frame 3:2
plot!(plot12,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3,title="I-data, gamma-outliers",ylabel="I, the number of infected", xlabel="t, days")

allindex=collect(1:1000)
outliersALL=intersect(setdiff(allindex,outliers),setdiff(allindex,outliersB),setdiff(allindex,outliersG))
#it is actually not outliers!!

#I-solution, all-outliers
if 1 in outliersALL
    plot13=plot(sol1,vars=(0,1),legend=false,xlims=[0.0,150.0],color="red",α =0.2)
else
    plot13=plot(sol1,vars=(0,1),legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    tspan=(t[1],t[end])
    prob=ODEProblem(f,[newdata[1],0],tspan,ParamData[j,:]) 
    sol=solve(prob,Vern9(),saveat=t)
    if j in outliersALL
        plot!(plot13,sol,vars=(0,1),legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    else
        plot!(plot13,sol,vars=(0,1),legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    end
end
plot13 #4:1
plot!(plot13,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3,title="I-solutions, all-outliers",ylabel="I, the number of infected", xlabel="t, days")

#I-data, all outliers
if 1 in outliersALL
    plot14=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
else
    plot14=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliersALL
        plot!(plot14,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    else
        plot!(plot14,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    end
end
plot14 #separately?
plot!(plot14,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3,title="I-data, all-outliers",ylabel="I, the number of infected", xlabel="t, days")

allplot=plot(plot7,plot8,plot10,plot9,plot11,plot12,plot13,plot14,layout=(4,2),size=(900,1200))
savefig(allplot,"simulations/plots/outliers_Ifit_25_1000.png")

##########



# comparison accross data and methods. Change appropriately the last part to include what you want
#change the indeces for parameters

#fitting to I
data1=dim[1,:] 
newdata1=float.(data1[sir(data1)[1]:sir(data1)[2]])
t1=collect(1.0:float(length(newdata1)))
outliers1=findall(x->x>percentile(N[:,1],75)+1.5*iqr(N[:,1])||x<percentile(N[:,1],25)-1.5*iqr(N[:,1]),N[:,1])
if 1 in outliers1 
    plot1=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot1=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers1
        plot!(plot1,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot1,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot1 #I fit
plot!(plot1,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="Fitting to I",titlefontsize=10,xlabel="") 

#fitting to R
outliers2=findall(x->x>percentile(M[:,1],75)+1.5*iqr(M[:,1])||x<percentile(M[:,1],25)-1.5*iqr(M[:,1]),M[:,1])

if 1 in outliers2
    plot2=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot2=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers2
        plot!(plot2,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot2,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot2 #R fit
plot!(plot2,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="Fitting to R",titlefontsize=10,xlabel="")

#fixed gamma
outliers3=findall(x->x>percentile(L[:,2],75)+1.5*iqr(L[:,1])||x<percentile(L[:,2],25)-1.5*iqr(L[:,2]),L[:,2])
if 1 in outliers3
    plot3=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot3=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers3
        plot!(plot3,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot3,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot3 #Fixed gamma
plot!(plot3,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="Fixed gamma",titlefontsize=10,xlabel="")

#beta via R0, note that you can leave it empty and skip the code if there are no results
plot4=plot(xlims=[0.0,75.0],legend=false,title="Beta via R0",titlefontsize=10,xlabel="")
outliers4=findall(x->x>percentile(S[:,2],75)+1.5*iqr(S[:,2])||x<percentile(S[:,2],25)-1.5*iqr(S[:,2]),S[:,2])
if 1 in outliers4
    plot4=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot4=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers4
        plot!(plot4,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot4,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot4 #Beta via r0
plot!(plot4,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="Beta via R0",titlefontsize=10,xlabel="")

#final size
outliers5=findall(x->x>percentile(T[:,1],75)+1.5*iqr(T[:,1])||x<percentile(T[:,1],25)-1.5*iqr(T[:,1]),T[:,1])
if 1 in outliers5
    plot5=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot5=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers5
        plot!(plot5,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot5,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot5 #final size
plot!(plot5,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="Final size",titlefontsize=15,xlabel="")

#combination
plot6=plot(xlims=[0.0,75.0],legend=false,title="Combination",titlefontsize=15,xlabel="")
outliers6=findall(x->x>percentile(K[:,1],75)+1.5*iqr(K[:,1])||x<percentile(K[:,1],25)-1.5*iqr(K[:,1]),K[:,1])
if 1 in outliers6
    plot6=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot6=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers6
        plot!(plot6,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot6,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot6 #Beta via r0
plot!(plot6,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="Combination",titlefontsize=15,xlabel="")

#imax
plot7=plot(xlims=[0.0,75.0],legend=false,title="Imax",titlefontsize=15,xlabel="")
outliers7=findall(x->x>percentile(R[:,1],75)+1.5*iqr(R[:,1])||x<percentile(R[:,1],25)-1.5*iqr(R[:,1]),R[:,1])
if 1 in outliers7
    plot7=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
else
    plot7=plot(t1,newdata1,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
end
for j in 2:1000
    datas=dim[j,:]
    newdata=float.(datas[sir(datas)[1]:sir(datas)[2]])
    t=collect(1.0:float(length(newdata)))
    if j in outliers7
        plot!(plot7,t,newdata,legend=false,xlims=[0.0,150.0],color="blue",α =0.2)
    else
        plot!(plot7,t,newdata,legend=false,xlims=[0.0,150.0],color="red",α =0.2)
    end
end
plot7 #Beta via r0
plot!(plot7,truesol,vars=(0,1),xlims=[0.0,75.0],legend=false,color="black",linewidth=3, title="Imax",titlefontsize=15,xlabel="")

#plot for 1 data only
allplot1=plot(plot1,plot2,plot3,layout=(3,1),size=(300,700),yrotation=45)
savefig(allplot1,"simulations/plots/outliers_gamma_25_1000.png") #combine accordingly
allplot1 #check

#after you worked through all option, combine to the final one
megaplot=plot(allplot1,allplot2,allplot3,allplot4,allplot5,allplot6,layout=(1,6),size=(1500,750)) #combine into one
savefig(megaplot,"simulations\\plots\\outliers_allmethods_gamma.png")