function [] = writeTifStack(stack,filepath)

% set up tif for intensity
t = Tiff(filepath,'w');
tagstruct.ImageLength = size(stack,1); % image height
tagstruct.ImageWidth  = size(stack,2); % image width
tagstruct.Photometric = Tiff.Photometric.MinIsBlack; % https://de.mathworks.com/help/matlab/ref/tiff.html
tagstruct.BitsPerSample = 32;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

% write away first frame (to overwrite rather than append to an existing stack)
setTag(t,tagstruct)
write(t,squeeze(uint32(stack(:,:,1))));

% read file frame by frame
for i = 2:size(stack,3)    
    writeDirectory(t);
    setTag(t,tagstruct)
    write(t,squeeze(uint32(stack(:,:,i)))); % Append
end
close(t);
end