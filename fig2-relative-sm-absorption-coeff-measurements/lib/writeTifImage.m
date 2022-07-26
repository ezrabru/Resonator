function [] = writeTifImage(image,path)
t = Tiff(path,'w');
tagstruct.ImageLength = size(image,1);
tagstruct.ImageWidth  = size(image,2); 
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
setTag(t,tagstruct); write(t,squeeze(uint32(image))); close(t);
end