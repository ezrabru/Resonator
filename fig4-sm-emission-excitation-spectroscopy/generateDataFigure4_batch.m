%% Batch analysis of 3-lobe resonator QDot data

clear all; close all; clc
addpath('lib')

% PSF displacement
%%%%%%%%%%%%%%%%
dx = 0;
dy = 37;

dx3 = 28;
dy3 = 16.5;

directory_tif = fullfile(pwd,'testdata'); % path to folder with tif stacks
directory_out = directory_tif;

NA = 1.49;
magnification = 1;
nImmersion = 1.518;
camOffset = 0;
camGain = 1;
camQE = 1;
camPixSize = 1;

bkgndEstimation.method = 'medfilt med';
bkgndEstimation.params.k = 50;
bkgndEstimation.params.frames = 10;

localisationParams.DoGsigmaSmall = 2;
localisationParams.DoGsigmaLarge = 3;
localisationParams.DoGminThreshold = 300;
localisationParams.method = 'asymmetric gaussian';
localisationParams.w = 7;
localisationParams.wiggle = 2;
localisationParams.wiggleDeg = 5;
localisationParams.lobeDist = sqrt(dx^2 + dy^2);

%% Prep

% Create new output folder
directory_out = fullfile(directory_out,'results');
if ~exist(directory_out,'dir'); mkdir(directory_out); end

% get filelist and organize
filelist = dir(fullfile(directory_tif,'*.tif'));


%% get data

Scope = Microscope(NA,magnification,nImmersion,camOffset,camGain,camQE,camPixSize);
Proc = Process(Scope,bkgndEstimation,localisationParams);
warning('off')
w = 3;
results = {};


% perform localisation on the MIP of the 405/405 nm dataset
filepath = fullfile(filelist.folder,filelist.name);
mip = Proc.getMaximumIntensityProjectionFromPath(filepath);
locs = Proc.processFrame(mip,0);
locs = struct2table(locs);

%%
clc

% split based on sigma into resonator and third lobes
threshold_sigmay = 1.7;
locs_resonator = locs(locs.sigmay > threshold_sigmay,:);
locs_dichroic = locs(locs.sigmay <= threshold_sigmay,:); locs_dichroic = table2struct(locs_dichroic);

% group the resonator lobes based on distance
locs_resonator = table2struct(locs_resonator);
[locs_upper,locs_lower] = Proc.pairLocalisationsByDistance(locs_resonator);
[locs_upper,locs_lower,dAngle,keep] = Proc.filterBasedOnAngle(locs_upper,locs_lower);
[locs_upper,locs_lower] = Proc.pairLocalisationsByDistance([locs_upper; locs_lower]);

% loop over each upper lobe and check if there is a matching lower (and
% that it is not accidentally a lower lobe...)
groupedLocs = struct;
numUpperLocs = numel(locs_upper);

%%
clc
figure;
imshow(mip,[]); hold on
for i = 1:numUpperLocs
    x_i = locs_upper(i).x/1e9;
    y_i = locs_upper(i).y/1e9;

    % expected position lower lobe
    x_i_ll = x_i + dy - 1;
    y_i_ll = y_i + dx;
    D = sqrt(([locs_lower.x]/1e9 - x_i_ll).^2 + ([locs_lower.y]/1e9 - y_i_ll).^2);
    [D_min_ll,id_D_min_ll] = min(D);
    if D_min_ll <= localisationParams.wiggle % lower lobe was found
        % expected position thir lobe
        x_i_dc = x_i + dy3;
        y_i_dc = y_i + dx3;
        D = sqrt(([locs_dichroic.x]/1e9 - x_i_dc).^2 + ([locs_dichroic.y]/1e9 - y_i_dc).^2);
        [D_min_3l,id_D_min_3l] = min(D);
        if D_min_3l <= 4*localisationParams.wiggle % third lobe was also found
            locsGrouped(i).upperLobe.x = x_i;
            locsGrouped(i).upperLobe.y = y_i;
            locsGrouped(i).lowerLobe.x = locs_lower(id_D_min_ll).x;
            locsGrouped(i).lowerLobe.y = locs_lower(id_D_min_ll).y;
            locsGrouped(i).thirdLobe.x = locs_dichroic(id_D_min_3l).x;
            locsGrouped(i).thirdLobe.y = locs_dichroic(id_D_min_3l).y;

            scatter(y_i,x_i,'r'); hold on
            scatter(locs_lower(id_D_min_ll).y/1e9,locs_lower(id_D_min_ll).x/1e9,'b'); hold on
            scatter(locs_dichroic(id_D_min_3l).y/1e9,locs_dichroic(id_D_min_3l).x/1e9,'y'); hold on
            pause(0.5)
        end
    end

end
legend('upper lobe','lower lobe','third lobe')


%%

x_lower = [locsGrouped.lowerLobe.x]/1e9;
y_lower = [locsGrouped.lowerLobe.y]/1e9;
x_upper = [locsGrouped.upperLobe.x];
y_upper = [locsGrouped.upperLobe.y];
x_third = [locsGrouped.thirdLobe.x]/1e9;
y_third = [locsGrouped.thirdLobe.y]/1e9;

% read stack
filepath = fullfile(filelist.folder,filelist.name);
bkgnd = Proc.estimateBackground(filepath);
figure;
subplot(1,3,1); histogram(bkgnd(:));
stack = Utils.readTiffStack(filepath);
subplot(1,3,2); histogram(stack(:));
stack = stack - bkgnd;
subplot(1,3,3); histogram(stack(:));
%%
clc
% measure intensities
[I_upperLobe,I_lowerLobe,I_thirdLobe] = measureIntensities(x_lower,y_lower,x_upper,y_upper,x_third,y_third,stack,w);

figure
plot(I_upperLobe); hold on
plot(I_lowerLobe); hold on
plot(I_thirdLobe); hold on


%%
figure;
imshow(mip,[]); hold on
scatter([locs_upper.y]/1e9,[locs_upper.x]/1e9); hold on
scatter([locs_lower.y]/1e9,[locs_lower.x]/1e9); hold on
scatter([locs_dichroic.y]/1e9,[locs_dichroic.x]/1e9);
legend('upper lobe','lower lobe','third lobe')


warning('on')

%% Functions

function [I_upperLobe,I_lowerLobe,I_thirdLobe] = measureIntensities(x_lower,y_lower,x_upper,y_upper,x_third,y_third,stack,w)
numLocs = numel(x_lower);
I_upperLobe = zeros(size(stack,3),numLocs);
I_lowerLobe = zeros(size(stack,3),numLocs);
I_thirdLobe = zeros(size(stack,3),numLocs);
for n=1:numLocs
    lobe1 = stack(round(x_lower(n))-w:round(x_lower(n))+w,round(y_lower(n))-w:round(y_lower(n))+w,:);
    lobe2 = stack(round(x_upper(n))-w:round(x_upper(n))+w,round(y_upper(n))-w:round(y_upper(n))+w,:);
    lobe3 = stack(round(x_third(n))-w:round(x_third(n))+w,round(y_third(n))-w:round(y_third(n))+w,:);
    intLobe1 = sum(lobe1,1:2); I_upperLobe(:,n) = intLobe1(:);
    intLobe2 = sum(lobe2,1:2); I_lowerLobe(:,n) = intLobe2(:);
    intLobe3 = sum(lobe3,1:2); I_thirdLobe(:,n) = intLobe3(:);
end
end
