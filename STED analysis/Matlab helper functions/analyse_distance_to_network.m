%% Distance detection: loops through all centriods of Ecad channel

im_sz=size(ch1);
dist=[];pos=[];
ang=[];
cnt=1;
xdist=[];
ydist=[];
weight=[];
apical_dist=[];
TJ_width=[];
sum_grids=zeros(size(cntr,1),1);
for x=1:size(cntr,1)
    if cntr(x,2)+grid/2 <= im_sz(1) && cntr(x,2)-grid/2 > 0 && cntr(x,1)+grid/2 <= im_sz(2) && cntr(x,1)-grid/2 > 0
        
        SubImg1 =ch1(  round(cntr(x,2)-grid/2) : round(cntr(x,2)+grid/2) ,round(cntr(x,1)-grid/2) : round(cntr(x,1)+grid/2)); 
        
        sum_grids(x)=sum(SubImg1(:));
        if sum(SubImg1(:))>tresh_grid %check for intesity tresholds in ZO1 channel
            
            % work with the segmented ZO1 channel
            SubImg1 =BW5(  round(cntr(x,2)-grid/2) : round(cntr(x,2)+grid/2) ,round(cntr(x,1)-grid/2) : round(cntr(x,1)+grid/2)); 
            %pxL=points(x).PixelList; % Pixel List
            pxL=cntr(x,:);
            ch2_temp=zeros(im_sz(1),im_sz(2));
            ch2_temp(pxL(:,2),pxL(:,1))=ch2(pxL(:,2),pxL(:,1));
%             figure(11)
%             imagesc(ch2_temp)
            
            SubImg2 =ch2_temp(  round(cntr(x,2)-grid/2) : round(cntr(x,2)+grid/2) ,round(cntr(x,1)-grid/2) : round(cntr(x,1)+grid/2)); 
            
%             figure(22)
%             C = imfuse(SubImg1,SubImg2);
%             imagesc(C)
%            
            %% k-nearest neighbors
            ch1_points=[];
            [ch1_points(:,1),ch1_points(:,2)]=find(SubImg1>0);
            ref_point(1)=(grid+2)/2;
            ref_point(2)=(grid+2)/2;
            
            [mIdx,mD] = knnsearch(ch1_points,ref_point,'K',1);
            if isempty(mIdx)==0 % position is close to ZO1 network
                %% distance and angle 
                xdiff=ch1_points(mIdx,2)-ref_point(2);
                ydiff=ch1_points(mIdx,1)-ref_point(1);
                dist(cnt)=sqrt(xdiff^2+ydiff^2);   % distance in pixels
                ang(cnt)=atan2d(xdiff,ydiff);      % angle in degree
                pos(cnt)=x;
                xdist(cnt)=xdiff;
                ydist(cnt)=ydiff;
                weight(cnt)=ch2(cntr(x,2),cntr(x,1));
                
                
                %% distance from apical (last ZO1 signal)
                % segmented network
                x_left=cntr(x,1)-20; x_right=cntr(x,1)+20;
                if x_left<1
                    x_left=1;
                elseif x_left>size(BW5,2)
                    x_left=size(BW5,2);
                end
                
                if x_right<1
                    x_left=1;
                elseif x_left>size(BW5,2)
                    x_right=size(BW5,2);
                end
                
                slice_netw=sum(BW5(:,x_left:x_right),2);
                apic_pos=find(slice_netw,1);
                basal_pos=find(slice_netw,1,'last');
                
                try
                    apical_dist(cnt)=cntr(x,2)-apic_pos;
                    TJ_width(cnt)=basal_pos-apic_pos;
                catch
                    apical_dist(cnt)=nan;
                    TJ_width(cnt)=nan;
                end
                
                
                
%                 figure(22)
%                 hold off
%                 imagesc(BW5)
%                 hold on
%                 plot(cntr(x,1),cntr(x,2),'r+')
%                 hold on
%                 plot([x_left,x_left],[0,600],'g-')
%                 plot([x_right,x_right],[0,600],'g-')
%                 drawnow
            else % position is out of network
                
                
                %% distance from apical (last ZO1 signal)
                % segmented network
                x_left=cntr(x,1)-20; x_right=cntr(x,1)+20;
                if x_left<1
                    x_left=1;
                elseif x_left>size(BW5,2)
                    x_left=size(BW5,2);
                end
                
                if x_right<1
                    x_left=1;
                elseif x_left>size(BW5,2)
                    x_right=size(BW5,2);
                end
                
                slice_netw=sum(BW5(:,x_left:x_right),2);
                apic_pos=find(slice_netw,1);
                basal_pos=find(slice_netw,1,'last');
                
                try
                    apical_dist(cnt)=cntr(x,2)-apic_pos;
                    TJ_width(cnt)=basal_pos-apic_pos;
                catch
                    apical_dist(cnt)=nan;
                    TJ_width(cnt)=nan;
                end
                
                dist(cnt)=nan;   % distance in pixels
                ang(cnt)=nan;      % angle in degree
                pos(cnt)=x;
                xdist(cnt)=nan;
                ydist(cnt)=nan;
                weight(cnt)=ch2(cntr(x,2),cntr(x,1));
                
%                 figure(22)
%                 hold off
%                 imagesc(BW5)
%                 hold on
%                 plot(cntr(x,1),cntr(x,2),'r+')
%                 hold on
%                 plot([x_left,x_left],[0,600],'g-')
%                 plot([x_right,x_right],[0,600],'g-')
%                 drawnow
                
            end
            cnt=cnt+1;
%             w = waitforbuttonpress;
        end
    end
    
end

% convert to nm
dist=dist*pix_size/upsampl;
xdist=xdist*pix_size/upsampl;
ydist=ydist*pix_size/upsampl;
apical_dist=apical_dist*pix_size/upsampl;
TJ_width=TJ_width*pix_size/upsampl;

%% plot results;
% calculate the intensity weights for all clusters 
all_weights=zeros(size(cntr,1),1);
for uu=1:size(cntr,1)
    all_weights(uu)=ch2(cntr(uu,2),cntr(uu,1));
end

% normalize weight based intensities
weights_norm=weight/max(weight);
weights_norm2=weight/min(weight);
all_weights_norm2=all_weights/min(all_weights);
cntr_in=cntr(pos,:);
rH=rH/pix_size*upsampl;

% figure(5)
% se1 = strel('disk',3,6);
% dilatedI = imdilate(BW5,se1);
% imagesc(dilatedI)
% axis image
% map=[0 0 0];
% map(2,:)=[1,0,1];
% colormap(map)
% hold on 
% 
% % figure(5)
% % [col2,row2]=find(BW5);
% % circles(row2,col2,rH,'facecolor','magenta','edgecolor','none')
% % set(gca,'Color','k')
% % axis ij
% % axis image
% % hold on 
% 
% %% version with different cluster size representing intensity weights.
% 
% %rH=rH/pix_size*upsampl;
% cntr_plot=cntr; % plot all positions
% weights_norm2_plot=all_weights_norm2; % plot all positions
% % %cntr_plot=cntr_in; % plot only posisitions in the network
% % %weights_norm2_plot=weights_norm2; % plot only posisitions in the network
% for ii=1:size(cntr_plot,1)
%     % weight rH by intensity of pixel
%     rHw=rH*weights_norm2_plot(ii)/1.5;
% %     pos2 = [cntr_plot(ii,1)-rHw cntr_plot(ii,2)-rHw 2*rHw 2*rHw]; 
% %     rectangle('Position',pos2,'Curvature',[1 1],'FaceColor','g','EdgeColor','g')
%     circles(cntr_plot(ii,1),cntr_plot(ii,2),rHw,'facecolor','green','facealpha',0.7,'edgecolor','none')
%     hold on
% end
% axis off
% xlim([0,size(BW5,2)])
% ylim([0,size(BW5,1)])
% set(gca,'Color','k')

%% version with different local cluster number representing intensity weights
% figure(555)
% se1 = strel('disk',3,6);
% dilatedI = imdilate(BW5,se1);
% imagesc(dilatedI)
% axis image
% map=[0 0 0];
% map(2,:)=[1,0,1];
% colormap(map)
% 
% hold on 
% cntr_plot=cntr; % plot all positions
% weights_norm2_plot=all_weights_norm2; % plot all positions
% %cntr_plot=cntr_in; % plot only posisitions in the network
% %weights_norm2_plot=weights_norm2; % plot only posisitions in the network
% for ii=1:size(cntr_plot,1)
%     pos2 = [cntr_plot(ii,1)-rH cntr_plot(ii,2)-rH 2*rH 2*rH]; 
%     rectangle('Position',pos2,'Curvature',[1 1],'FaceColor','g','EdgeColor','g')
%     hold on
%     
%     % weight rH by intensity of pixel
%     for u=1:round(weights_norm2_plot(ii))-1
%         rx=(-2*rand(1)+1)*rH*1.5;
%         ry=(-2*rand(1)+1)*rH*1.5;
%         pos2 = [cntr_plot(ii,1)+rx-rH cntr_plot(ii,2)+ry-rH 2*rH 2*rH]; 
%         rectangle('Position',pos2,'Curvature',[1 1],'FaceColor','g','EdgeColor','g')
%     end
%    
%     
% end
% axis off

% fname=[filename(1:end-4),'-Segmented_Ch1+Ch2'];
% print(fname,'-dpdf','-fillpage')
% fname=[filename(1:end-4),'-Segmented_ZO1.tif'];
% imwrite(BW5*256,fname,'tif')
% fname=[filename(1:end-4),'-ZO1.tif'];
% imwrite(uint8(ch1/max(ch1(:))*256),fname,'tif')


% figure(3)
% imagesc(ch2)
% axis image
% colormap gray
% hold on 
% plot(cntr(pos,1),cntr(pos,2),'g.','MarkerSize',10)
% 
% figure(34)
% plot(dist,apical_dist,'.')
% 
% figure(35)
% hist(apical_dist,30)
% 
% figure(7) 
% hist(dist,20)
% xlabel('Distance [nm]')
% ylabel('Counts')
% set(gca,'FontSize',22)

results.dist=dist;
results.weight=weight;
results.apical_dist=apical_dist;
results.ang=ang;
results.TJ_width=TJ_width;
results.ImSize=size(ch1)*pix_size; 



