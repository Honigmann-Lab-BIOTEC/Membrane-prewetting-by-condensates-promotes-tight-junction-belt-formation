
clear all;
close all;

name='PATJ'; % protein name


% open all distance results files
files={};
[file,PathName,FilterIndex] = uigetfile('*mat*', 'mat file','MultiSelect','on');
cd(PathName)

if iscell(file)==0
    files{1}=file;
else
    files=file;
end

all_dist=[];
all_weight=[];
all_apical_dist=[];
all_TJ_width=[];
for i = 1:size(files,2)
    load (files{i})
    
    all_dist=[all_dist, results.dist];
    all_apical_dist=[all_apical_dist, results.apical_dist];
    all_weight=[all_weight, results.weight/min(results.weight)];
    all_TJ_width=[all_TJ_width, results.TJ_width];
 
end   


%% expand the distance values to account for the weights ww == with weights
all_dist_ww=[];all_apical_dist_ww=[];all_TJ_width_ww=[];
for b=1:size(all_weight,2)
    temp=all_dist(b);
    temp=repmat(temp,1,round(all_weight(b)));
    all_dist_ww=[all_dist_ww,temp];
    
    temp=all_apical_dist(b);
    temp=repmat(temp,1,round(all_weight(b)));
    all_apical_dist_ww=[all_apical_dist_ww,temp];
    
    temp=all_TJ_width(b);
    temp=repmat(temp,1,round(all_weight(b)));
    all_TJ_width_ww=[all_TJ_width_ww,temp];
    
end



h = figure(3);
histogram(all_apical_dist_ww,200,'EdgeColor',[0,0,0],'FaceColor',[0,0,0],'FaceAlpha',1)
xlim([-500,1500])
set(gca,'Fontsize',20)
xlabel('Distance to apical ZO1 strand [nm]')
ylabel('Counts')
title('Apical Distance') 
savefig(h,'Apical_Distance_Distribution.fig')

     
savename=[name '_pooled_distance_results.mat'];
save(savename,'all_apical_dist_ww','all_dist_ww','all_TJ_width_ww','name')
