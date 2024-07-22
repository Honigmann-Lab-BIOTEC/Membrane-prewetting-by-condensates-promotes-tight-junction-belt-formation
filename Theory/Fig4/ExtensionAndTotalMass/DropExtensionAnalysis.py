from matplotlib import pyplot as plt
import numpy as np
import time
import matplotlib.animation as animation
import matplotlib.cm as cm
import warnings
import matplotlib as mlib
from matplotlib import ticker
warnings.filterwarnings("error")

###Plotting parameters###
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 8  # width in inches
fig_height = 8*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
###This defines plot parameters###
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


###
### Parameters from the numerical calculation
###
RandomSeedVar=1
np.random.seed(RandomSeedVar)
a=0.01

kappa=np.sqrt(a)
dt1=1e-2
dt2=1e-2

###NucleationAndGrowth####
t1=20.
t2=300.
###Fast growth###
t1=200.
t2=200.
tTotal=t1+t2

valuemin=-1.125419456270195
valuemax=1.545009253060331
gamma=1.

p=1
mx=3
my=1/2
nodesx=int(128*mx)
nodesy=int(128*my)

N_it1=int(t1/dt1)
N_it2=int(t2/dt2)
N_it=N_it1+N_it2
N_frames=500
Saving_Frames=2000
frame_every=int(N_it/N_frames)
frame_every_saving=int(N_it/Saving_Frames)
data_saving_every=frame_every

mu_linear_array=np.load("Mu_Linear_Array_Values.npy")
x, dx = np.linspace(0,1*mx,nodesx,retstep=True)
y, dy = np.linspace(0,1*my,nodesy,retstep=True)
dts=tTotal/2000.
fig,ax = plt.subplots(figsize=(8,8))
# ax.set_ylim(-5,25)
# ax.set_xlim(0,40)
ax.tick_params(axis='x', pad=14)
ax.tick_params(axis='y', pad=14)
# ax.set_yticks([0,5,10,15,20,25])
# ax.set_xticks([0,10,20,30,40])
ax.set_xlabel(r'$\mathrm{time \ [min]}$',size=32)
ax.set_ylabel(r'$\mathrm{extension \ [}\mu\mathrm{m]}$',size=32)

AllExtension_Array=[]
AllTime_Array=[]
velocities=[]
j_aux=0

DataForAnalysis = np.load("TrialData_Speed_StrongerNuc_ForPaper_"+str(j_aux)+".npy") ###Takes "TrialData_Speed_StrongerNuc_ForPaper_0.npy" as an input for the calculations
LinearAnalysisMidCoordinate=[]
for i in range(len(DataForAnalysis)):
    LinearAnalysisMidCoordinate.append(DataForAnalysis[i][32])

InitialTime=1
InitialTimeIndex=int(InitialTime/dts)
InitialTimeIndex=0

MaximimValueReachedInEveryNucleatedDrop=InitialTimeIndex
# t0=InitialTime+300*dts
# t1=t0+300*dts
# t2=t1+300*dts
# t3=t2+300*dts

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

InitialIndex=10

for i in range(InitialIndex,len(LinearAnalysisMidCoordinate),1):
    ConcentrationArray_Trial=LinearAnalysisMidCoordinate[i]
    NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial)
    if GridPositionInterfaceLeft_Trial>184:
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
        break
    else:
        Time_Steps.append(index-index_trial)
        Extensions_PerDrop=np.vstack((Extensions_PerDrop,Extensions))
        pass

for l in range(InitialNumberOfDrops):
    ScaleInX=32/384. 
    ScaleInT=40./137.3
    CurrentVel=(ScaleInX/ScaleInT)*((Extensions_PerDrop[-2]-Extensions_PerDrop[int(len(Extensions_PerDrop)/2)])/(dts*Time_Steps[-1]-dts*Time_Steps[int(len(Extensions_PerDrop)/2)]))
    velocities.append(CurrentVel)
    AllExtension_Array.append(ScaleInX*(Extensions_PerDrop[1:-1]-Extensions_PerDrop[1]))
    AllTime_Array.append(dts*ScaleInT*np.array((Time_Steps[1:])))
    # np.savetxt("Extension_"+str(j_aux)+".txt",(ScaleInX*(Extensions_PerDrop[1:-1]-Extensions_PerDrop[1])[0:687]))
    np.savetxt("Extension_"+str(j_aux)+".txt",(ScaleInX*(Extensions_PerDrop[1:-1]-Extensions_PerDrop[1])))                                                                                        
    plt.plot(dts*ScaleInT*np.array((Time_Steps[1:])),ScaleInX*(Extensions_PerDrop[1:-1]-Extensions_PerDrop[1]),color='g',linewidth=4,alpha=0.3)
        
# np.savetxt("ExtensionTime.txt",AllTime_Array[-1][0:687])
np.savetxt("ExtensionTime.txt",AllTime_Array[-1])




#plt.savefig("ExtensionVsTime.pdf",bbox_inches='tight')
#plt.savefig("ExtensionVsTime.svg",bbox_inches='tight')
plt.show()

