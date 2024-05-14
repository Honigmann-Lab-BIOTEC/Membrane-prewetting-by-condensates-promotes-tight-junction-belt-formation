function stack=read_tiff_stack(fname)

    info=imfinfo(fname);
    W=info(1).Width;
    H=info(1).Height;
    
    stack=zeros(H,W,size(info,1));
    for i=1:size(info,1)
        stack(:,:,i)=imread(fname,i);
    end
    
end