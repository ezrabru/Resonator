%% Batch semi-automated analysis of 2-lobe resonator QDot data

clear all
close all
clc
addpath('lib');

directory_tif  = fullfile(pwd,'testdata_qd655'); % path to folder with tif stacks
directory_out  = directory_tif;

pathCalibration = fullfile(pwd,'testdata_qd655','QD655_roi1_B405_T40500000.tif');

%% Read data and get max intensity projections

% Create new output folder
outputdir = fullfile(directory_out,'results'); if ~exist(outputdir, 'dir'); mkdir(outputdir); end
dir_MIP   = fullfile(outputdir,'MIPs'); if ~exist(dir_MIP, 'dir'); mkdir(dir_MIP); end

% Get list of tif files in folder
filelist = dir(fullfile(directory_tif,'*.tif'));
disp(['Number of tif files found: ' num2str(length(filelist))])
for i = 1:length(filelist); disp("  "+filelist(i).name); end; disp(' ')


%% THIS STEP REQUIRES USER INPUT: selection of templates for localisation

% Generate maximum intensity projection
stack = readTifStack(pathCalibration);
MIP = max(stack,[],3);
clear variable stack

% Let user pick rois
imshow(MIP,[]); set(gcf,'Position',get(0,'Screensize'));
roi = drawrectangle;
x = floor(roi.Position(1)); y = floor(roi.Position(2));
w = ceil(roi.Position(3));  h = ceil(roi.Position(4));
roi = MIP(y:y+h,x:x+w);
imshow(roi,[]); pause(1);
psfCalibration = roi;


%% THIS STEP REQUIRES USER INPUT: Accepting/rejecting localizations
% Localizes using cross-correlation with selected experimental PSFs.
% Reopen the files one by one and let user manually accept or reject
% localisations by typing 'y' or 'n' in command window.

treshold = 0.7;

psf = psfCalibration;
c = normxcorr2(psf,MIP);

% offset found by correlation
c(c < treshold) = 0;
c = medfilt2(c,[3,3]);

% get peaks areas and centroids
stats = regionprops(logical(c),c,'Area','WeightedCentroid');
locs = [stats.WeightedCentroid]';
locs = reshape(locs',2,length(locs)/2)';
x = locs(:,1) - size(psf,2)/2;
y = locs(:,2) - size(psf,1)/2;
keep = ones(size(x));
for j=1:length(x)
    imshow(MIP,[]); hold on
    plot(x(j),y(j),'+r')
    pause(1)
    answer = input('Accept localisation (y/n)?','s');
    if strcmp(answer,'n')
        keep(j) = 0;
    end
end
coordinates.x = x(logical(keep));
coordinates.y = y(logical(keep));

x_rejected = x(~logical(keep));
y_rejected = y(~logical(keep));

% Annotate MIP with chosen rois
x = coordinates.x;
y = coordinates.y;
imshow(MIP,[]); hold on
for j = 1:length(x)
    rectangle('Position',[x(j)-size(psf,2)/2,y(j)-size(psf,1)/2,size(psf,2),size(psf,1)],'EdgeColor','y');
    text(x(j)-size(psf,2)/2,y(j),string(j),'Color','y');
end

%% THIS STEP REQUIRES USER INPUT: determine distance between lobes

% Let user crop out a psf
fig = figure;
imshow(MIP,[]); hold on
roi = drawrectangle;
x = floor(roi.Position(1)); y = floor(roi.Position(2));
w = ceil(roi.Position(3));  h = ceil(roi.Position(4));
close(fig);
roi = MIP(y:y+h,x:x+w);
fig = figure; imshow(roi,[]); pause(1);
    
% Let user localize lobes manually repeatibly
locs = ginput(inf); close(fig);

distLobesPix = [];
for i=1:2:size(locs,1)
    dx = sqrt(sum((locs(i,1) - locs(i+1,1)).^2));
    dy = sqrt(sum((locs(i,2) - locs(i+1,2)).^2));
    distLobesPix = [distLobesPix; dx dy];
end

dx = mean(distLobesPix(:,1));
dy = mean(distLobesPix(:,2));


%% Estimate background

medfiltKernel = 30;

for i = 1:length(filelist)
    
    results(i).filename = filelist(i).name;
    stack = readTifStack(fullfile(directory_tif,results(i).filename));
    bkgnd = medfilt2(median(stack,3),[medfiltKernel medfiltKernel],'symmetric');
    
    subplot(1,3,1); imshow(bkgnd,[0 65500]);
    subplot(1,3,2); imshow(stack(:,:,1),[0 65500]);
    stackMinBkgnd = stack(:,:,1) - bkgnd; stackMinBkgnd(stackMinBkgnd < 0) = 0;
    subplot(1,3,3); imshow(stackMinBkgnd,[]);
    pause(0.1)
    
    results(i).bkgnd = bkgnd;
    writeTifImage(MIP,fullfile(dir_MIP,strcat(extractBefore(filelist(i).name,'.'),'_bkgnd.tif')))

end


%% Measure intensities

windowIntensity = 7; w = ceil(windowIntensity/2);

for i = 1:length(filelist)
    
    x = coordinates.x;
    y = coordinates.y;
    numLocs = length(x);
    
    stack = readTifStack(fullfile(directory_tif,results(i).filename));
    stack = stack - results(i).bkgnd;
%     stack(stack < 0) = 0;
    
    intensitiesLobe1 = zeros(size(stack,3),numLocs);
    intensitiesLobe2 = zeros(size(stack,3),numLocs);
    for j=1:numLocs
        lobe1 = stack(round(y(j)-dy/2)-w:round(y(j)-dy/2)+w,round(x(j)-dx/2)-w:round(x(j)-dx/2)+w,:);
        lobe2 = stack(round(y(j)+dy/2)-w:round(y(j)+dy/2)+w,round(x(j)+dx/2)-w:round(x(j)+dx/2)+w,:);
        
        intLobe1 = sum(lobe1,1:2); intensitiesLobe1(:,j) = intLobe1(:);
        intLobe2 = sum(lobe2,1:2); intensitiesLobe2(:,j) = intLobe2(:);
    end

    results(i).intLobe1 = intensitiesLobe1;
    results(i).intLobe2 = intensitiesLobe2;
    results(i).ratio12 = intensitiesLobe1./intensitiesLobe2;
    results(i).n = length(intensitiesLobe1(:));
    
end

save('results.mat','results')




%% Plot of median spectra for each quantum dot

numLocs = numel(coordinates.x);
wavelenghts = [405 488 561 638];
factorPower = [1 1 1 1];
correctionFactor = factorPower.*(405./wavelenghts);

intensityThreshold = 5e4;

meanRatio405 = [];
meanRatio488 = [];
meanRatio561 = [];
meanRatio638 = [];

for i=1:numLocs
    
    keep1 = logical(results(1).intLobe1(:,i) > intensityThreshold);
    ratio1 = results(1).ratio12(:,i);
    ratio1 = nanmean(ratio1(keep1));
    meanRatio405 = [meanRatio405 ratio1];
    
    keep2 = logical(results(2).intLobe1(:,i) > intensityThreshold);
    ratio2 = results(2).ratio12(:,i);
    ratio2 = nanmean(ratio2(keep2));
    meanRatio488 = [meanRatio488 ratio2];
    
    keep3 = logical(results(3).intLobe1(:,i) > intensityThreshold);
    ratio3 = results(3).ratio12(:,i);
    ratio3 = nanmean(ratio3(keep3));
    meanRatio561 = [meanRatio561 ratio3];
    
    keep4 = logical(results(4).intLobe1(:,i) > intensityThreshold);
    ratio4 = results(4).ratio12(:,i);
    ratio4 = nanmean(ratio4(keep4));
    meanRatio638 = [meanRatio638 ratio4];
    
    ratios = [ratio1 ratio2 ratio3 ratio4];
    ratios = correctionFactor./ratios;
    plot(wavelenghts,ratios,'-o','color',0.8*[1 1 1])
    hold on
end
xlim([395 650])
hold on
meanRatios = [nanmean(meanRatio405) nanmean(meanRatio488) nanmean(meanRatio561) nanmean(meanRatio638)];
meanRatios = correctionFactor./meanRatios;
plot(wavelenghts,meanRatios,'b-o')

grid on
xlabel('Wavelength (nm)');
ylabel('Ratio intensity lobe2/lobe1');
set(gca,'FontSize',10)
savefig(fullfile(outputdir,'estimatedExcitationSpectra.fig'))

%% Plot ratio histograms for individual QDots

close all

binEdges = 0:0.05:1.5;

numLocs = numel(coordinates.x);
wavelenghts = [405 488 561 638];
factorPower = [1 1 1 1];
correctionFactor = factorPower.*(405./wavelenghts);

intensityThreshold = 10e4;

for i=1:numLocs
    
    figure('position',[50 50 1600 200]);
    
    keep1 = logical(results(1).intLobe1(:,i) > intensityThreshold);
    ratio1 = results(1).ratio12(:,i);
    ratio1 = ratio1(keep1); numel(ratio1)
    subplot(1,4,1); histogram(1./ratio1(:),binEdges)
    
    keep2 = logical(results(2).intLobe1(:,i) > intensityThreshold);
    ratio2 = results(2).ratio12(:,i);
    ratio2 = ratio2(keep2); numel(ratio2)
    subplot(1,4,2); histogram(1./ratio2(:),binEdges)
    
    keep3 = logical(results(3).intLobe1(:,i) > intensityThreshold);
    ratio3 = results(3).ratio12(:,i);
    ratio3 = ratio3(keep3); numel(ratio3)
    subplot(1,4,3); histogram(1./ratio3(:),binEdges)
    
    keep4 = logical(results(4).intLobe1(:,i) > intensityThreshold);
    ratio4 = results(4).ratio12(:,i);
    ratio4 = ratio4(keep4); numel(ratio4)
    subplot(1,4,4); histogram(1./ratio4(:),binEdges)
    
end

