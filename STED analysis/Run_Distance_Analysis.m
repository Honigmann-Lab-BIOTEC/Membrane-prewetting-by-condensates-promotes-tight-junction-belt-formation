clear all; close all


%% input parameters
pix_size=15; % pixel size of image in nm


%% get filenames of image
files={};
[file,PathName,FilterIndex] = uigetfile('*.*', 'Tif stack','MultiSelect','off');
cd(PathName)

%% load image stack
disp(['loading: ',file])
out=import_multichannel_tif_stack(file);
ch2=double(out{1}); % ZO1 channel
ch1=double(out{2}); % PATJ channel

%% upsample image to 3 nm
upsampl=5; % factor for upsampling
ch2 = imresize(ch2,upsampl);
ch1 = imresize(ch1,upsampl);

%% segment tight junction ZO1
disp('segmenting ZO1 channel')
% input variables
m=9; % local window for contrast tresholding in pixel 
contrast=2; % contrast trehold
min_area=200; % minimal area of connected pixels
mean_int_tresh=55; % minimum intensity of segments
fill_holes_tresh=80; % teshold in area pixel for holes to be filled

BW3=segment_and_treshold(ch1,m,contrast,min_area,mean_int_tresh,fill_holes_tresh);

%% skeletonize segmented ZO1 channel
spur_iter=25; % iterations to remove spurios pixel of skeleton
BW4 = bwmorph(BW3,'skel',inf);
BW5 = bwmorph(BW4,'spur',spur_iter); % remove small prunes
BW5(ch1<mean_int_tresh)=0;


%% Segment PATJ clusters
disp('segmenting PATJ channel')
m=15; % local window for contrast tresholding in pixel 
contrast=3; % contrast trehold
min_area=2; % minimal area of connected pixels
mean_int_tresh=2; % minimum intensity of segments
fill_holes_tresh=80; % teshold in area pixel for holes to be filled
shrink_iter=3; % iterations to remove spurios pixel of skeleton

trs=20; % minimum mean intensity for Ecad cluster
tm=1; % 1 means median filter 0 mean_filter
grid=40; % window size in pixel for filtering image
dist_tresh=30; % distance treshold for clusters in ch2 in nm
dist_tresh=dist_tresh*upsampl/pix_size; % in pixel now

find_regional_maxima_and_filter

%% Analyze PATJ cluster distances to nearest ZO1 strand
disp('analysing PATJ distance to ZO1')

tresh_grid=0; % minimum intensity in the ZO1 channel to evaluate Ch2
rH=10; % max size of plotted cluster size

analyse_distance_to_network

%% save results
disp('saving results')
folderName = 'Analysis Results'; % Specify the name of the folder
mkdir(folderName);
CF = pwd;

ZO1_network=BW5;
PATJ_centers=cntr;
sname = [CF,'/',folderName,'/',file(1:end-4), '-Dist_Analysis.mat'];
save(sname,"results","ZO1_network","PATJ_centers")

%% plot overlay of segmentation and cluster positions
C=imfuse(BW5,ch1);
figure(1)
imagesc(C);
axis image
axis off
drawnow
hold on
plot(cntr(:,1),cntr(:,2),'r.')
