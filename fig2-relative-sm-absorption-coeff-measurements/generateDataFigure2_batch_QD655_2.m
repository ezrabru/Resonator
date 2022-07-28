%% Batch semi-automated analysis of 2-lobe resonator QDot data

clear all; close all; clc
addpath('lib2')

% PSF displacement
dx = 0;
dy = 37;
lobeDist = sqrt(dx^2 + dy^2);
wiggle = 1;

numRois = 15;
%directory_tif = 'D:\manuscripts\resonator\new data\20201024_resonator quantitative\QD585'; % path to folder with tif stacks
directory_tif = fullfile(pwd,'testdata_qd585'); % path to folder with tif stacks
directory_out = directory_tif;

flagCalib = 'B405_T405';

flag405 = 'B405_T405';
flag488 = 'B405_T488';
flag561 = 'B405_T561';
flag638 = 'B405_T638';

NA = 1.49;
magnification = 1;
nImmersion = 1.518;
camOffset = 0;
camGain = 1;
camQE = 1;
camPixSize = 1;

bkgndEstimation.method = 'min of avg';

localisationParams.DoGsigmaSmall = 2;
localisationParams.DoGsigmaLarge = 3;
localisationParams.DoGminThreshold = 500;
localisationParams.method = 'asymmetric gaussian';
localisationParams.w = 7;

%% Prep

% Create new output folder
outputdir = fullfile(directory_out,'results'); if ~exist(outputdir, 'dir'); mkdir(outputdir); end

% get filelist and organize
roiPaths = {};
filelist = dir(fullfile(directory_tif,'*.tif'));
roi = 1;
for i=1:numRois
    flagRoi = strcat('roi',num2str(i));
    for j=1:length(filelist)
        if strfind(filelist(j).name,flagRoi) % get roi
            if strfind(filelist(j).name,flag405)
                roiPaths(roi).filename405 = filelist(j).name;
            elseif strfind(filelist(j).name,flag488)
                roiPaths(roi).filename488 = filelist(j).name;
            elseif strfind(filelist(j).name,flag561)
                roiPaths(roi).filename561 = filelist(j).name;
            elseif strfind(filelist(j).name,flag638)
                roiPaths(roi).filename638 = filelist(j).name;
            end
            if strfind(filelist(j).name,flagCalib)
                roiPaths(roi).filenameCalib = filelist(j).name;
            end
        end
    end
    roi = roi + 1;
end


%% get data

id_roi = 1;

results(id_roi).roiID = id_roi;

% perform a regular SMLM localisation on the calibration dataset
pathCalibStack = fullfile(directory_tif,roiPaths(id_roi).filenameCalib);

Scope = Microscope(NA,magnification,nImmersion,camOffset,camGain,camQE,camPixSize);
Proc = Process(Scope,bkgndEstimation,localisationParams);

img = double(imread(pathCalibStack));
bkgnd = Proc.estimateBackground(pathCalibStack);
locs = Proc.processFrame(img,bkgnd);

%% Group localisations

x = [locs.x]';
y = [locs.y]';
numLocs = numel(x);
id = 1:numLocs;
group_id = nan(numLocs,1);

D = squareform(pdist([x y]/1e9,'euclidean'));

% % remove localisations that are too close together or overlapping
% minDist = localisationParams.w;
% D(D < minDist) = nan;

% remove localisations that can't be paired with another one considering
% the PSF spatial arrangement
D((D < lobeDist - wiggle)) = 0;
D((D > lobeDist + wiggle)) = 0;

% remove half of the pdist matrix (because of symmetry)
id_lower_lobe = tril(D);
id_upper_lobe = triu(D);

[row,col] = find(triu(D));
coord_upper_lobe = [x(row) y(col)]/1e9;

[row,col] = find(tril(D));
coord_lower_lobe = [x(row) y(col)]/1e9;

group_id = 1:numel(x(row));

figure
% imshow(img,[]); hold on
scatter(coord_upper_lobe(:,2),coord_upper_lobe(:,1),100,group_id,'s'); hold on
scatter(coord_lower_lobe(:,2),coord_lower_lobe(:,1),100,group_id,'o');
colormap(lines)

%% Functions

function [psf,MIP] = getMIPandPSFcutOut(filepath)

% Generate maximum intensity projection
stack = readTifStack(filepath);
MIP = max(stack,[],3);

% Let user pick rois
imshow(MIP,[]); set(gcf,'Position',get(0,'Screensize'));
roi = drawrectangle;
x = floor(roi.Position(1)); y = floor(roi.Position(2));
w = ceil(roi.Position(3));  h = ceil(roi.Position(4));
roi = MIP(y:y+h,x:x+w);
imshow(roi,[]); pause(1);
psf = roi;

end

function coordinates = localizeCrossCorr(img,psf,threshold)

c = normxcorr2(psf,img);

% offset found by correlation
c(c < threshold) = 0;
c = medfilt2(c,[3,3]);

% get peaks areas and centroids
stats = regionprops(logical(c),c,'Area','WeightedCentroid');
locs = [stats.WeightedCentroid]';
locs = reshape(locs',2,length(locs)/2)';
x = locs(:,1) - size(psf,2)/2;
y = locs(:,2) - size(psf,1)/2;

keep = ones(size(x));
imshow(img,[]); hold on
plot(x,y,'+r')
pause(1)
answer = input('Manually accept/reject localisations (y/n)?','s');
if strcmp(answer,'y')
    for j=1:length(x)
        imshow(img,[]); hold on
        plot(x(j),y(j),'+r')
        pause(1)
        answer = input('Accept localisation (y/n)?','s');
        if strcmp(answer,'n')
            keep(j) = 0;
        end
    end
end
coordinates.x = x(logical(keep));
coordinates.y = y(logical(keep));

% x_rejected = x(~logical(keep));
% y_rejected = y(~logical(keep));
%
% % Annotate MIP with chosen rois
% x = coordinates.x;
% y = coordinates.y;
% imshow(MIP,[]); hold on
% for j = 1:length(x)
%     rectangle('Position',[x(j)-size(psf,2)/2,y(j)-size(psf,1)/2,size(psf,2),size(psf,1)],'EdgeColor','y');
%     text(x(j)-size(psf,2)/2,y(j),string(j),'Color','y');
% end
end

function [dx,dy] = getRelativeCoordinatesLobes(psf)

fig = figure; imshow(psf,[]); pause(1);
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
end

function bkgnd = estimateBackground(filepath,medfiltKernel)

stack = readTifStack(filepath);
bkgnd = medfilt2(median(stack,3),[medfiltKernel medfiltKernel],'symmetric');

end

function [intensitiesLobe1,intensitiesLobe2] = measureIntensitiesLobes(filepath,w,dx,dy,coordinates)

% Estimate background
medfiltKernel = 30;
bkgnd = estimateBackground(filepath,medfiltKernel);

x = coordinates.x;
y = coordinates.y;
numLocs = length(x);

stack = readTifStack(filepath);
stack = stack - bkgnd;

intensitiesLobe1 = zeros(size(stack,3),numLocs);
intensitiesLobe2 = zeros(size(stack,3),numLocs);
for j=1:numLocs
    lobe1 = stack(round(y(j)-dy/2)-w:round(y(j)-dy/2)+w,round(x(j)-dx/2)-w:round(x(j)-dx/2)+w,:);
    lobe2 = stack(round(y(j)+dy/2)-w:round(y(j)+dy/2)+w,round(x(j)+dx/2)-w:round(x(j)+dx/2)+w,:);
    
    intLobe1 = sum(lobe1,1:2); intensitiesLobe1(:,j) = intLobe1(:);
    intLobe2 = sum(lobe2,1:2); intensitiesLobe2(:,j) = intLobe2(:);
end
end