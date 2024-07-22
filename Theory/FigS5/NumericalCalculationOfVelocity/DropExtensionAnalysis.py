from matplotlib import pyplot as plt
import numpy as np
import time
import matplotlib.animation as animation
import matplotlib.cm as cm
import warnings
import matplotlib as mlib
from matplotlib import ticker
warnings.filterwarnings("error")
from scipy.optimize import curve_fit

###Plotting parameters###
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 8  # width in inches
fig_height = 8*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
#fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(16,8*golden_mean/.9))

font = {'family' : 'serif',
'color'  : 'black',
    'weight' : 'normal',
        'size'   : 32,
    }

params = {'text.usetex': True,
'axes.labelsize': 32,
    'font.size': 32,
        'legend.fontsize': 32,
            'axes.linewidth' : 1.5,
                'lines.linewidth': 1.5,
                    'lines.markersize':7,
                        'lines.markeredgewidth':1.2,
                            'figure.figsize': fig_size,
}
pp=np.array(["#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#9ACD32","#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#9ACD32"])
fmt=['+', 'x', 'o']
plt.rcParams.update(params)
valuemin=-1.125419456270195
valuemax=1.545009253060331
ScaleInX=25/300.
ScaleInT=40./137.3

###
### Parameters from the numerical calculation
###
RandomSeedVar=1
np.random.seed(RandomSeedVar)
a=0.01
kappa=np.sqrt(a)

dt1=1e-2
dt2=1e-2
###Fast growth###
t1=200.
t2=200.
tTotal=t1+t2

gamma=1.

p=1
mx=3
my=1/2
nodesx=int(128*mx)
nodesy=int(128*my)

N_it1=int(t1/dt1)
N_it2=int(t2/dt2)
N_it=N_it1+N_it2
N_frames=1000
N_frames=500
Saving_Frames=4000
Saving_Frames=2000
frame_every=int(N_it/N_frames)
frame_every_saving=int(N_it/Saving_Frames)
data_saving_every=frame_every
velocities=[]
velocities_fit=[]
error_fit=[]
mu_linear_array=np.load("MuParameters_Checking.npy") ###This loads previous values of the binding affinity to which I add +0.3 for epsilon_0 
mu_linear_array=np.hstack((mu_linear_array,np.array((-0.4,-0.3))))
x, dx = np.linspace(0,1*mx,nodesx,retstep=True)
y, dy = np.linspace(0,1*my,nodesy,retstep=True)
dts=tTotal/2000.
fig,ax = plt.subplots(figsize=(8,8))

# ax.set_xlim(0,90)
ax.tick_params(axis='x', pad=14)
ax.tick_params(axis='y', pad=14)

ax.set_xlabel(r'$\mathrm{time \ [}\mathrm{min]}$',size=32)
ax.set_ylabel(r'$\mathrm{extension \ [\mu}\mathrm{m]}$',size=32)
j_aux=0

DataForAnalysis = np.load("TrialData_Speed_StrongerNuc_ForPaper_"+str(j_aux)+".npy")

LinearAnalysisMidCoordinate=[]
for i in range(len(DataForAnalysis)):
    LinearAnalysisMidCoordinate.append(DataForAnalysis[i][32])
InitialTime=1
InitialTimeIndex=int(InitialTime/dts)
InitialTimeIndex=0

MaximimValueReachedInEveryNucleatedDrop=InitialTimeIndex


MaxValueOfTheSystem=np.amax(LinearAnalysisMidCoordinate[MaximimValueReachedInEveryNucleatedDrop])
MinValueOfTheSystem=np.amin(LinearAnalysisMidCoordinate[MaximimValueReachedInEveryNucleatedDrop])
MeanValueToTrackInterface=(MaxValueOfTheSystem+MinValueOfTheSystem)/2.
TrialIndex=MaximimValueReachedInEveryNucleatedDrop
ConcentrationArray_Trial=LinearAnalysisMidCoordinate[TrialIndex]

def FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray):
    Indices_Trial=(ConcentrationArray>MeanValueToTrackInterface).nonzero()[0]
    DiffInIndices=Indices_Trial[1::]-Indices_Trial[0:-1:]
    IndicesToTraceBackJumps=(DiffInIndices>1.1).nonzero()[0]
    NumberOfDrops=len(IndicesToTraceBackJumps)+1

    RightInterfaceIndices_Plus=np.append(Indices_Trial[IndicesToTraceBackJumps[::]],Indices_Trial[-1])
    RightInterfaceIndices_Minus=RightInterfaceIndices_Plus+1

    LeftInterfaceIndices_Plus=np.append(Indices_Trial[0],Indices_Trial[1+IndicesToTraceBackJumps[::]])
    LeftInterfaceIndices_Minus=LeftInterfaceIndices_Plus-1

    ConcentrationInRightInterfaceIndices_Plus=ConcentrationArray[RightInterfaceIndices_Plus]
    ConcentrationInRightInterfaceIndices_Minus=ConcentrationArray[RightInterfaceIndices_Minus]

    ConcentrationInLeftInterfaceIndices_Plus=ConcentrationArray[LeftInterfaceIndices_Plus]
    ConcentrationInLeftInterfaceIndices_Minus=ConcentrationArray[LeftInterfaceIndices_Minus]

    SlopesInConcentrationRight=(ConcentrationInRightInterfaceIndices_Plus-ConcentrationInRightInterfaceIndices_Minus)/(RightInterfaceIndices_Plus-RightInterfaceIndices_Minus)

    OrdinateToOriginRight=ConcentrationInRightInterfaceIndices_Plus-SlopesInConcentrationRight*RightInterfaceIndices_Plus

    SlopesInConcentrationLeft=(ConcentrationInLeftInterfaceIndices_Plus-ConcentrationInLeftInterfaceIndices_Minus)/(LeftInterfaceIndices_Plus-LeftInterfaceIndices_Minus)

    OrdinateToOriginLeft=ConcentrationInLeftInterfaceIndices_Plus-SlopesInConcentrationLeft*LeftInterfaceIndices_Plus

    GridPositionInterfaceRight=(MeanValueToTrackInterface-OrdinateToOriginRight)/SlopesInConcentrationRight
    GridPositionInterfaceLeft=(MeanValueToTrackInterface-OrdinateToOriginLeft)/SlopesInConcentrationLeft

    Extensions=GridPositionInterfaceRight-GridPositionInterfaceLeft
    return NumberOfDrops,GridPositionInterfaceLeft,GridPositionInterfaceRight,Extensions

InitialIndex=0
for i in range(InitialIndex,len(LinearAnalysisMidCoordinate),1):
    ConcentrationArray_Trial=LinearAnalysisMidCoordinate[i]
    NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial)
    if GridPositionInterfaceLeft_Trial>184:####175:
        dummy=1
    else:
        index_trial=i
        break

NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial)

IndicesAux=np.arange(index_trial,len(LinearAnalysisMidCoordinate),1)

InitialConcentrationArray=LinearAnalysisMidCoordinate[index_trial]

NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,InitialConcentrationArray)
InitialNumberOfDrops=NumberOfDrops

aux_counter=0
Time_Steps=[]
Extensions_PerDrop=Extensions
for index in IndicesAux:
    ConcentrationArray_Trial=LinearAnalysisMidCoordinate[index]
    NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial)

    aux_counter+=1
    if NumberOfDrops!=InitialNumberOfDrops or (GridPositionInterfaceLeft_Trial<1) or (GridPositionInterfaceRight_Trial>381):
        #print(aux_counter)
        break
    else:
        Time_Steps.append(index-index_trial)
        Extensions_PerDrop=np.vstack((Extensions_PerDrop,Extensions))
        pass

for l in range(InitialNumberOfDrops):
    mu_plot_number=-mu_linear_array[j_aux] -0.3
    mu_plot="%.2f" % mu_plot_number
    ScaleInX=25/300.
    ScaleInT=40./137.3
    CurrentVel=(ScaleInX/ScaleInT)*((Extensions_PerDrop[-2]-Extensions_PerDrop[1])/(dts*Time_Steps[-1]-dts*Time_Steps[1]))
    velocities.append(CurrentVel)
    #if j_aux==0:
        #dummy=1
    #else:
    plt.plot(ScaleInT*dts*np.array((Time_Steps[1:])),ScaleInX*(Extensions_PerDrop[1:-1]-Extensions_PerDrop[1]),linewidth=4,label=str(mu_plot),color=pp[j_aux])
    slope, intercept = np.polyfit(ScaleInT*dts*np.array((Time_Steps[1:])), ScaleInX*(Extensions_PerDrop[1:-1]-Extensions_PerDrop[1]), deg=1)
    def linearFun(xvar,m):
        return m*xvar
    ExtensionValues=(ScaleInX*np.array((Extensions_PerDrop[1:-1]-Extensions_PerDrop[1]))).flatten() 
    popt,pcov=curve_fit(linearFun,ScaleInT*dts*np.array((Time_Steps[1:])),ExtensionValues) ###Calculates a linear fit to extract the velocity of growth of the surface-nucleated condensate
    velocities_fit.append(popt)
    error_fit.append(np.sqrt(pcov[0]))


ax.legend(loc="upper left",frameon=False, handlelength = 1,handletextpad=0.3,fontsize=22,ncol=2,labelspacing=0.3,columnspacing=0.5,title=r'${\Delta \bar \epsilon}$')
#plt.savefig("DifferentSpeedsAsAFunctionOfMuStar.pdf",bbox_inches='tight')
plt.show()


BindingAffinities=np.array((-mu_linear_array-0.3))[::-1] ###Here I write the values of the relative binding affinities following the convention in the manuscript
SortedVelocities=velocities[::-1]
SortedVelocities=velocities_fit[::-1]
BindingAffinities2=BindingAffinities
velocities2=SortedVelocities
np.savetxt("FittingErrors.txt",error_fit[::-1])
np.savetxt("Velocity.txt",velocities2)
np.savetxt("BindingAffinities.txt",BindingAffinities2)
#plt.savefig("VelocityVsBindingAffinityvP.svg",bbox_inches='tight')
#plt.savefig("VelocityVsBindingAffinityvP.pdf",bbox_inches='tight')
