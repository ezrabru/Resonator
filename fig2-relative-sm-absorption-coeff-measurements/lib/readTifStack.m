function stack = readTifStack(filepath)
info   = imfinfo(filepath);
width  = info(1).Width;
height = info(1).Height;
frames = size(info,1);
stack = zeros(height,width,frames);
for i=1:frames
    stack(:,:,i) = imread(filepath,i);
end
end