function [montage,maskBorders] = generateMontageFromStack(stack,numPixBorder)


[nx,ny] = size(stack(:,:,1,1));

numColsMontage = size(stack,3);
numRowsMontage = size(stack,4);

montage = [];
for i = 1:numRowsMontage
    newRow = [];
    for j = 1:numColsMontage
        newRow = [newRow, nan(nx,numPixBorder), stack(:,:,j,i)];
    end
    montage = [montage; nan(numPixBorder,size(newRow,2)); newRow];
end

montage = montage(numPixBorder+1:end,numPixBorder+1:end);

maskBorders = isnan(montage);
montage(maskBorders) = 1;

end