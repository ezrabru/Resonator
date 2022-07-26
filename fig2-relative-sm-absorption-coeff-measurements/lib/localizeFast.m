%% localizeFast.m
% 
% This function was based on the mathworks entry 'FastPeakFind.m' by:
%   Adi Natan (natan@stanford.edu)
%   Ver 1.7 , Date: Oct 10th 2013

function locs = localizeFast(img,threshold,filter,edge,resolution)

% if 'img' is a double, convert to uint16 to speed up
if isfloat(img); img = uint16(img); end

% median filter
img = medfilt2(img,[3,3]);

% apply threshold
img = img.*uint16(img > threshold);

if any(img(:)) % proceed if at least one pixel is above threshold
    
    % smooth image
    img = conv2(single(img),filter,'same') ;
    
    switch resolution % switch between local maxima and sub-pixel methods
        
        case 1 % peak find - using the local maxima approach - 1 pixel resolution
            
            % img will be noisy on the edges, and also local maxima looks
            % for nearest neighbors so edge must be at least 1. We'll skip 'edge' pixels.
            sd = size(img);
            [x,y] = find(img(edge:sd(1)-edge,edge:sd(2)-edge));
            
            % initialize outputs
            locs = [];            
            x = x + edge - 1;
            y = y + edge - 1;
            for j=1:length(y)
                if (img(x(j), y(j)) >= img(x(j)-1, y(j)-1 )) &&...
                   (img(x(j), y(j)) >  img(x(j)-1, y(j)))    &&...
                   (img(x(j), y(j)) >= img(x(j)-1, y(j)+1))  &&...
                   (img(x(j), y(j)) >  img(x(j)  , y(j)-1))  && ...
                   (img(x(j), y(j)) >  img(x(j)  , y(j)+1))  && ...
                   (img(x(j), y(j)) >= img(x(j)+1, y(j)-1))  && ...
                   (img(x(j), y(j)) >  img(x(j)+1, y(j)))    && ...
                   (img(x(j), y(j)) >= img(x(j)+1, y(j)+1))
                    locs = [locs; [y(j) x(j)]];
                end
            end
            
        case 2 % find weighted centroids of processed image,  sub-pixel resolution.

            % get peaks areas and centroids
            stats = regionprops(logical(img),img,'Area','WeightedCentroid');
            locs = [stats.WeightedCentroid]';            
    end   
    
else % if all pixels are lower than threshold
    locs=[];
    return
end
end