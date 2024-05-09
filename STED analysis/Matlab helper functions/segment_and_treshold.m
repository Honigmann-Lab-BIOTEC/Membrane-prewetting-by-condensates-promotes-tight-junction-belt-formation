
function BW3=segment_and_treshold(ch1,m,contrast,min_area,mean_int_tresh,fill_holes_tresh)
% % input variables
% ch1 = input image
% m=5; % local window for contrast tresholding in pixel 
% contrast=2; % contrast trehold
% min_area=200; % minimal area of connected pixels
% mean_int_tresh=30; % minimum intensity of segments
% fill_holes_tresh=80; % teshold in area pixel for holes to be filled

%% local intensity tresholding 
BW1 = logical(bernsen(ch1,[m m],contrast));

%% filter TJ segmentation with area and intensity
s = regionprops(BW1,ch1,'Area','MeanIntensity','PixelList');
BW2=BW1;

for i=1:size(s,1)
    if s(i).Area < min_area || s(i).MeanIntensity < mean_int_tresh
        pxL=s(i).PixelList; % Pixel List
        BW2(pxL(:,2),pxL(:,1))=0;
    end
end

SE = strel('disk',2,6);
BW2 = imdilate(BW2,SE);
BW2 = bwmorph(BW2,'close');

%% fill little whole in TJ segmentation based on area tresh
fill_tresh=fill_holes_tresh;
filled = imfill(BW2, 'holes');
holes = filled & ~BW2;

L_candidate_holes = bwlabel(holes); 
stats = regionprops(L_candidate_holes, 'Area'); 
idx = find([stats.Area]<fill_tresh); 
midsized_holes = ismember(L_candidate_holes,idx);
BW3=BW2;
BW3(midsized_holes==1)=1;


