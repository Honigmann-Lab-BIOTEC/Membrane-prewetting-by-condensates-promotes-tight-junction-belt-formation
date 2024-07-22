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


###
### Parameters from the numerical calculation
###
ScaleInX=25/300.
ScaleInT=40./137.3
RandomSeedVar=1
np.random.seed(RandomSeedVar)
a=0.01
#a=0.1
kappa=np.sqrt(a)
dt1=1e-2
dt2=1e-2

t1=40.
t2=40.
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
N_frames=1000
N_frames=500
Saving_Frames=4000
Saving_Frames=2000
frame_every=int(N_it/N_frames)
frame_every_saving=int(N_it/Saving_Frames)
data_saving_every=frame_every

mu_linear_array=np.load("Mu_Linear_Array_Values.npy")
x, dx = np.linspace(0,1*mx,nodesx,retstep=True)
y, dy = np.linspace(0,1*my,nodesy,retstep=True)
dts=tTotal/2000.

AllExtension_Array=[]
AllTime_Array=[]
velocities=[]


DataForAnalysis = np.load("MultipleDrops.npy")
LinearAnalysisMidCoordinate=[]
for i in range(len(DataForAnalysis)):
    LinearAnalysisMidCoordinate.append(DataForAnalysis[i][32])

InitialTime=1
InitialTimeIndex=int(InitialTime/dts)
InitialTimeIndex=0


MaximimValueReachedInEveryNucleatedDrop=InitialTimeIndex+300
t0=InitialTime+300*dts
t1=t0+300*dts
t2=t1+300*dts
t3=t2+300*dts


# MaxValueOfTheSystem=np.amax(LinearAnalysisMidCoordinate[MaximimValueReachedInEveryNucleatedDrop])
# MinValueOfTheSystem=np.amin(LinearAnalysisMidCoordinate[MaximimValueReachedInEveryNucleatedDrop])
valuemin=-1.125419456270195
MinValueOfTheSystem=valuemin
MaxValueOfTheSystem=valuemax
MeanValueToTrackInterface=(MaxValueOfTheSystem+MinValueOfTheSystem)/2.
TrialIndex=MaximimValueReachedInEveryNucleatedDrop
# ConcentrationArray_Trial=LinearAnalysisMidCoordinate[TrialIndex]

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

InitialIndex=0 ###Wait until all initial condensates have relaxed and fused if they are close enough. -The data sample already begins after all condensates have nucleated and fused if they were close enough.

# InitialIndex=20 ###Wait until all initial condensates have relaxed and fused if they are close enough. -The data sample already begins after all condensates have nucleated and fused if they were close enough.

ConcentrationArray_Trial=LinearAnalysisMidCoordinate[InitialIndex]
NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial)
IndicesAux=np.arange(InitialIndex,len(LinearAnalysisMidCoordinate)-InitialIndex,1)
InitialNumberOfDrops=NumberOfDrops


AllExtensionsPerDrop=[]
AllTimeSteps=[]
Time_Steps=[]
index=InitialIndex
for j in range(NumberOfDrops):
    aux_counter=0
    Time_Steps=[]
    ConcentrationArray_Trial=LinearAnalysisMidCoordinate[index]
    NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial)
    IndicesAux=np.arange(index,len(LinearAnalysisMidCoordinate),1)
    InitialNumberOfDrops=NumberOfDrops
    Extensions_PerDrop=Extensions
    if index==len(LinearAnalysisMidCoordinate)-1:
        break
    for index in IndicesAux:
        ConcentrationArray_Trial=LinearAnalysisMidCoordinate[index]
        NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial)
        aux_counter+=1
        if NumberOfDrops!=InitialNumberOfDrops or (np.sum(GridPositionInterfaceLeft_Trial<1)) or ( np.sum(GridPositionInterfaceRight_Trial>381)):
            final_counter=aux_counter
            break
        else:
            Time_Steps.append((index-InitialIndex))
            Extensions_PerDrop=np.vstack((Extensions_PerDrop,Extensions))
            pass
    AllExtensionsPerDrop.append(ScaleInX*Extensions_PerDrop)
    AllTimeSteps.append(ScaleInT*dts*np.array((Time_Steps)))


fig,ax = plt.subplots(figsize=(8,8))


ax.tick_params(axis='x', pad=14)
ax.tick_params(axis='y', pad=14)


ax.set_xlabel(r'$\mathrm{time \ [min]}$',size=32)
ax.set_ylabel(r'$\mathrm{extension \ [}\mu\mathrm{m]}$',size=32)

for j in range(len(AllExtensionsPerDrop)):
    for k in range(len(AllExtensionsPerDrop[j][0])):
        np.savetxt("TimeForNumOfDrops_"+str(j)+".txt",AllTimeSteps[j])
        np.savetxt("AllExtensionsPerDrop_"+str(j)+"_"+str(k)+".txt",AllExtensionsPerDrop[j][1:,k])
        ax.plot(AllTimeSteps[j],AllExtensionsPerDrop[j][1:,k],color=pp[k])
        

#plt.savefig("MultiDrop.pdf",bbox_inches='tight')
#plt.savefig("MultiDrop.svg",bbox_inches='tight')
plt.show()
