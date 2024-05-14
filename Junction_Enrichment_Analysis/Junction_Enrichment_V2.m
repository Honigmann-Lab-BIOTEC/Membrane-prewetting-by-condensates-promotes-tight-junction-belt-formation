clear all; close all;

%% load data
timeres=10; % frame rate of the image stack in min

datpath=uigetdir;
cd(datpath)

seg_TJ_file='ZO1_TJ-Mask.tif';
stack_TJ=read_tiff_stack(seg_TJ_file);
stack_TJ=logical(stack_TJ);

cell_file='CellPose Masks.tif';
stack_cells=read_tiff_stack(cell_file);

ZO1_file='ZO1_MIPs.tif';
stack_ZO1=read_tiff_stack(ZO1_file);

Client_file='MAGI3_MIPs.tif';
stack_Client=read_tiff_stack(Client_file);



%% denoise segmented TJ 
min_area=5; % minimum size of TJ region (remove noise)
stack_TJ_filt=stack_TJ;
for i=1:size(stack_TJ,3)
    temp_im=stack_TJ(:,:,i);
    temp_im(temp_im>0)=1;
    figure(1)
    imagesc(temp_im)
    axis image
    
    BW2=area_treshold_binary(temp_im,min_area);
    figure(2)
    imagesc(BW2)
    axis image
    
    stack_TJ_filt(:,:,i)=BW2;
end

%% measure TJ enrichment per cell over time
se = strel('disk',4);
JEF_ZO1_TP=[];
JEF_Client_TP=[];
for i=1:size(stack_cells,3)
    % current time point of all stacks
    TP_cell=stack_cells(:,:,i);
    TP_TJ=stack_TJ_filt(:,:,i);
    TP_ZO1=stack_ZO1(:,:,i);
    TP_Client=stack_Client(:,:,i);
    
    TP_over=TP_cell.*TP_TJ; % overlap regions between Cell and TJ segmentation
    figure(3)
    imagesc(TP_over)
    axis image
    axis off
    title(sprintf('TJs at TP:%d',i))
    drawnow
    
    
    cell_IDs=unique(TP_over); % all IDs of cells with TJ
   
    
    for ii=2:size(cell_IDs,1)
        fprintf('Working on TP:%d Cell:%d \n', i,ii)
        temp=zeros(size(TP_cell));
        temp(TP_cell==cell_IDs(ii))=1;
%         figure(6)
%         imagesc(temp)
        
        %% get TJ of this cell
        J = imdilate(temp,se);
        TP_over=J.*TP_TJ; % overlap regions between Cell and TJ segmentation
%         figure(3)
%         imagesc(TP_over)
        
        if sum(TP_over(:))>0 % check that this cell has a TJ
            % measure intensity in the ZO1 channel
            TJ_ZO1=TP_over.*TP_ZO1; 
            Int_TJ_ZO1=mean(TJ_ZO1(TJ_ZO1>0));
            
            % measure intensity in the Client channel
            TJ_Client=TP_over.*TP_Client;
            Int_TJ_Client=mean(TJ_Client(TJ_Client>0));
    %         figure(9)
    %         imagesc(TJ_temp);

             %% get CP of this cell
            J = imdilate(TP_over,se);
            J = imcomplement(J);
            CP_over=J.*temp; % overlap regions between Cell and TJ segmentation
%             figure(4)
%             imagesc(CP_over)

            CP_ZO1=CP_over.*TP_ZO1;
            Int_CP_ZO1=mean(CP_ZO1(CP_ZO1>0));
            
            CP_Client=CP_over.*TP_Client;
            Int_CP_Client=mean(CP_Client(CP_Client>0));

    %         figure(10)
    %         imagesc(CP_temp);

            %% Junction vs Cytoplasm Ratio
            % ZO1 channel
            R(1)=Int_TJ_ZO1/Int_CP_ZO1;
            R(2)=i;
            JEF_ZO1_TP=[JEF_ZO1_TP;R];
            
            % Client channel
            R(1)=Int_TJ_Client/Int_CP_Client;
            R(2)=i;
            JEF_Client_TP=[JEF_Client_TP;R];
        
        end
    end
end

%% plot results 
JEF_ZO1_mean=[];
JEF_ZO1_sem=[];
JEF_Client_mean=[];
JEF_Client_sem=[];
for i=1:size(stack_cells,3)
    ind=find(JEF_ZO1_TP(:,2)==i);
    JEF_ZO1_mean(i)=nanmean(JEF_ZO1_TP(ind,1));
    JEF_ZO1_sem(i)=nanstd(JEF_ZO1_TP(ind,1))/sqrt(length(ind));
    
    JEF_Client_mean(i)=nanmean(JEF_Client_TP(ind,1));
    JEF_Client_sem(i)=nanstd(JEF_Client_TP(ind,1))/sqrt(length(ind));
end

fh(1)=figure(5);
%plot(JEF_TP(:,2),JEF_TP(:,1),'.')
hold on
x=(1:max(JEF_ZO1_TP(:,2)))*timeres;
errorbar(x, JEF_ZO1_mean,JEF_ZO1_sem,'g','LineWidth',3)
title('ZO1 Enrichment')
xlabel('Time [min]')
ylabel('Junction Enrichment Factor')
set(gca,'FontSize',18)

fh(2)=figure(6);
errorbar(x, JEF_Client_mean,JEF_Client_sem,'r','LineWidth',3)
title('Client Enrichment')
xlabel('Time [min]')
ylabel('Junction Enrichment Factor')
set(gca,'FontSize',18)

fh(3)=figure(7);
JEF_ZO1_mean_norm=JEF_ZO1_mean-min(JEF_ZO1_mean(1:3));
JEF_ZO1_mean_norm=JEF_ZO1_mean_norm/max(JEF_ZO1_mean_norm)+1;
JEF_Client_mean_norm=JEF_Client_mean-min(JEF_Client_mean(1:3));
JEF_Client_mean_norm=JEF_Client_mean_norm/max(JEF_Client_mean_norm)+1;
hold on

errorbar(x, JEF_ZO1_mean_norm,JEF_ZO1_sem,'g','LineWidth',3)
title('ZO1 Enrichment')

errorbar(x, JEF_Client_mean_norm,JEF_Client_sem,'r','LineWidth',3)
title('Normalized Enrichment')

legend('ZO1','Client','Location','SouthEast')
xlabel('Time [min]')
ylabel('Normalized Junction Enrichment')
set(gca,'FontSize',18)

%% save results
cd ..
mkdir('Analysis Results')
cd('Analysis Results')
savefig(fh(1),'JEF_ZO1.fig')
savefig(fh(2),'JEF_Client.fig')
savefig(fh(3),'JEF_Normalized.fig')
t=x;
save('Results.mat','JEF_Client_mean_norm','JEF_Client_sem','JEF_ZO1_mean_norm','JEF_ZO1_sem','JEF_Client_TP','JEF_ZO1_TP','JEF_ZO1_mean','JEF_Client_mean','t')