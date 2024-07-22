from matplotlib import pyplot as plt
import numpy as np
import time
from scipy.fft import rfft2, irfft2, fftfreq, rfftfreq, rfft, irfft,fft,ifft,fft2,ifft2
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
pp=np.array(["#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#9ACD32","#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#9ACD32","#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7","#9ACD32"])
fmt=['+', 'x', 'o']
plt.rcParams.update(params)

def FindInterfaceGridPositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray,MeanValueToTrackInterface):
    Indices_Trial=(ConcentrationArray>MeanValueToTrackInterface).nonzero()[0]
    DiffInIndices=Indices_Trial[1::]-Indices_Trial[0:-1:]
    IndicesToTraceBackJumps=(DiffInIndices>1.1).nonzero()[0]
    NumberOfDrops=len(IndicesToTraceBackJumps)+1

    RightInterfaceIndices_Plus=np.append(Indices_Trial[IndicesToTraceBackJumps[::]],Indices_Trial[-1])
    RightInterfaceIndices_Minus=RightInterfaceIndices_Plus+1

    LeftInterfaceIndices_Plus=np.append(Indices_Trial[0],Indices_Trial[1+IndicesToTraceBackJumps[::]])
    LeftInterfaceIndices_Minus=LeftInterfaceIndices_Plus-1
    return RightInterfaceIndices_Plus,RightInterfaceIndices_Minus,LeftInterfaceIndices_Plus,LeftInterfaceIndices_Minus



def FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray,MeanValueToTrackInterface):
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



###
### Parameters from the numerical calculation
###
ScaleInX=25/300.
ScaleInT=40./137.3
RandomSeedVar=1
np.random.seed(RandomSeedVar)
a=0.01

kappa=np.sqrt(a)
dt1=1e-2
dt2=1e-2

t1=40.
t2=40.
tTotal=t1+t2
tTotal=t1+t2
valuemin=-1.125419456270195
valuemax=1.545009253060331
gamma=1.


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


mu_linear_array=np.load("Mu_Linear_Array_Values.npy")

x, dx = np.linspace(0,1*mx,nodesx,retstep=True)
y, dy = np.linspace(0,1*my,nodesy,retstep=True)
dts=tTotal/2000.

AllLuminiscence_Array=[]
AllTime_Array=[]

DataForAnalysis = np.load("MultipleDrops.npy")
LinearAnalysisMidCoordinate=[]
for i in range(len(DataForAnalysis)):
    LinearAnalysisMidCoordinate.append(DataForAnalysis[i][32])

InitialTime=1
InitialTimeIndex=int(InitialTime/dts)
MaximimValueReachedInEveryNucleatedDrop=InitialTimeIndex+300
t0=InitialTime+300*dts
t1=t0+300*dts
t2=t1+300*dts
t3=t2+300*dts
MaxValueOfTheSystem=valuemax
MinValueOfTheSystem=valuemin

InitialIndex=0 ###Wait until all initial condensates have relaxed and fused if they are close enough. -The data sample already begins after all condensates have nucleated and fused if they were close enough.

# InitialIndex=20 ###Wait until all initial condensates have relaxed and fused if they are close enough. -The data sample already begins after all condensates have nucleated and fused if they were close enough.

MeanValueToTrackInterface=(valuemax+valuemin)/2.
TrialIndex=InitialIndex
index_trial=InitialIndex
ConcentrationArray_Trial=LinearAnalysisMidCoordinate[TrialIndex]

PositionsOfInterface=FindInterfaceGridPositions(valuemax,valuemin,ConcentrationArray_Trial,MeanValueToTrackInterface)


ConcentrationArray_Trial=LinearAnalysisMidCoordinate[InitialIndex]


NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,ConcentrationArray_Trial,MeanValueToTrackInterface)

IndicesAux=np.arange(index_trial,len(LinearAnalysisMidCoordinate),1)

InitialConcentrationArray=LinearAnalysisMidCoordinate[index_trial]

NumberOfDrops,GridPositionInterfaceLeft_Trial,GridPositionInterfaceRight_Trial,Extensions=FindInterfacePositions(MaxValueOfTheSystem,MinValueOfTheSystem,InitialConcentrationArray,MeanValueToTrackInterface)
InitialNumberOfDrops=NumberOfDrops


ConcentrationArray2D=DataForAnalysis[index_trial]
AuxMatrix=np.zeros(ConcentrationArray2D.shape)
IndicesZerosAndOnes=np.nonzero(ConcentrationArray2D>MeanValueToTrackInterface)
AuxMatrix[IndicesZerosAndOnes]=1
InitialLuminiscence=np.sum(ConcentrationArray2D*AuxMatrix)
MeanValueToTrackInterface=0

AllDropLuminiscence=[]
AllTime=[]
for index in IndicesAux:
    ConcentrationArray_Trial=LinearAnalysisMidCoordinate[index]
    PositionsOfInterface=FindInterfaceGridPositions(valuemax,valuemin,ConcentrationArray_Trial,MeanValueToTrackInterface)
    ConcentrationArray2D=DataForAnalysis[index]
    AuxDropLuminiscence=[]
    AllTime.append(dts*ScaleInT*index)
    for DropIndex in range(len(PositionsOfInterface[0])):
        IndexLeft=PositionsOfInterface[3][DropIndex]
        IndexRight=PositionsOfInterface[0][DropIndex]
        DropValues=ConcentrationArray2D[:,IndexLeft:IndexRight]
        AuxMatrix=np.zeros(DropValues.shape)
        IndicesZerosAndOnes=np.nonzero(DropValues>MeanValueToTrackInterface)
        AuxMatrix[IndicesZerosAndOnes]=1
        LuminiscenceAux=np.sum(DropValues*AuxMatrix)
        AuxDropLuminiscence.append(LuminiscenceAux)
    AllDropLuminiscence.append(AuxDropLuminiscence)


NewDropIndices=[]
InitialNumberOfDrops
auxNumberOfDrops=InitialNumberOfDrops
for timeIndex in range(len(IndicesAux)):
    NewNumberOfDrops=len(AllDropLuminiscence[timeIndex])
    if abs(auxNumberOfDrops-NewNumberOfDrops)>0:
        NewDropIndices.append(timeIndex)
        auxNumberOfDrops=NewNumberOfDrops


fig,ax = plt.subplots(figsize=(8,8))
# ax.set_xlim(0,12)

ax.tick_params(axis='x', pad=14)
ax.tick_params(axis='y', pad=14)
ax.set_xticks([0,4,8,12])

ax.set_xlabel(r'$\mathrm{time \ [min]}$',size=32)
ax.set_ylabel(r'$\mathrm{total \ intensity \ 10^{3} \ [A.U.]}$',size=32)

FirstTwoTimeLimits=np.array((0,NewDropIndices[0]-1))
AuxTwoTimeLimits=FirstTwoTimeLimits
for j in range(len(NewDropIndices)):
    Luminiscence_j=np.array((AllDropLuminiscence[AuxTwoTimeLimits[0]:AuxTwoTimeLimits[1]]))
    Time_j=np.array((AllTime[AuxTwoTimeLimits[0]:AuxTwoTimeLimits[1]]))
    for k in range(len(AllDropLuminiscence[AuxTwoTimeLimits[0]])):
        np.savetxt("TimeForNumOfDropsLum_"+str(j)+".txt",Time_j)
        np.savetxt("AllLuminiscencesPerDrop_"+str(j)+"_"+str(k)+".txt",Luminiscence_j[:,k])
        ax.plot(Time_j,Luminiscence_j[:,k],linewidth=4,color=pp[k])
    if (j+1)!=len(NewDropIndices):
        AuxTwoTimeLimits=np.array((NewDropIndices[j],NewDropIndices[j+1]-1))
    else:
        AuxTwoTimeLimits=np.array((NewDropIndices[j],timeIndex))

Luminiscence_j=np.array((AllDropLuminiscence[AuxTwoTimeLimits[0]:AuxTwoTimeLimits[1]]))
Time_j=np.array((AllTime[AuxTwoTimeLimits[0]:AuxTwoTimeLimits[1]]))
for k in range(len(AllDropLuminiscence[AuxTwoTimeLimits[0]])):
    np.savetxt("TimeForNumOfDropsLum_"+str(j)+".txt",Time_j)
    np.savetxt("AllLuminiscencesPerDrop_"+str(j)+"_"+str(k)+".txt",Luminiscence_j[:,k])
    ax.plot(Time_j,Luminiscence_j[:,k],linewidth=4,color=pp[k])

#plt.savefig("MultipleDropIntensities_Trial.svg")
#plt.savefig("MultipleDropIntensities_Trial.pdf")
plt.show()