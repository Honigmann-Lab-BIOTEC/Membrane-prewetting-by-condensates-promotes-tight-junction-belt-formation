%% find regional maxima in PATJ channel
bw1 = bernsen(ch2,[m m],contrast);

% reomve small prunes
bw1 = bwmorph(bw1,'shrink',shrink_iter);
bw1 = bwmorph(bw1,'majority');

%% filter segmentation with area and intensity
bw1=logical(bw1);
s = regionprops(bw1,ch2,'Area','MeanIntensity','PixelList');
bw2=bw1;
for i=1:size(s,1)
    if s(i).MeanIntensity < trs
        pxL=s(i).PixelList; % Pixel List
        bw2(pxL(:,2),pxL(:,1))=0;
    end
end

bw3 = bwmorph(bw2,'thin',15);


%% estimate emitter positions
intensity_tresh=trs; % minimal pixel intesnity of emitter position
[col, row]=find(bw3==1);
cntr = row;
cntr(:,2)=col;

% sort pixels based on intensity
intens=zeros(size(cntr,1),1);
for iii=1:size(cntr,1)
    intens(iii)=ch2(cntr(iii,2),cntr(iii,1));
end
[Ints,Indis]=sort(intens,'descend');

ind_tresh=Ints>intensity_tresh;
cntr=cntr(Indis(ind_tresh),:);


distances = pdist2(cntr, cntr);
distances(distances==0)=1000; 
closePairs = distances < dist_tresh;

%% remove localization which are too close to each other
out_ind=1;
for ii = 1:size(closePairs,1)
    [row,col]=find(closePairs(:,ii));
    if ismember(ii,out_ind)==0  
        out_ind=[out_ind,row'];
    end
end

%delete too close neigbors
cntr(out_ind(2:end),:)=[];
figure(3)
hold off
imagesc(ch2)
colormap gray
axis image
axis off
hold on
plot(cntr(:,1),cntr(:,2),'r+')
