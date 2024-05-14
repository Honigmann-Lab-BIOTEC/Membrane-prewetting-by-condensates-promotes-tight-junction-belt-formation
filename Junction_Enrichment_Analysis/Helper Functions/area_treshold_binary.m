
function BW2=area_treshold_binary(im,min_area)


s = regionprops(im,'Area','PixelList');
BW2=im;
for i=1:size(s,1)
    if s(i).Area < min_area 
        pxL=s(i).PixelList; % Pixel List
        BW2(pxL(:,2),pxL(:,1))=0;
    end
end
