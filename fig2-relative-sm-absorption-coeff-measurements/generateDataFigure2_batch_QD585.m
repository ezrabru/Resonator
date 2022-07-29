%% Batch semi-automated analysis of 2-lobe resonator QDot data

clear all; close all; clc
addpath('lib2')

% PSF displacement
dx = 0;
dy = 37;

numRois = 14;
directory_tif = 'D:\manuscripts\resonator\new data\20201024_resonator quantitative\QD585'; % path to folder with tif stacks
% directory_tif = fullfile(pwd,'testdata_qd585'); % path to folder with tif stacks
directory_out = directory_tif;

pathBulkSpectrum = 'D:\manuscripts\resonator\figures\figure 2\quantitative resonator\data spectra\dyes\Qdot585.txt';

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

bkgndEstimation.method = 'medfilt med';
bkgndEstimation.params.k = 50;
bkgndEstimation.params.frames = 10;

localisationParams.DoGsigmaSmall = 2;
localisationParams.DoGsigmaLarge = 3;
localisationParams.DoGminThreshold = 1000;
localisationParams.method = 'asymmetric gaussian';
localisationParams.w = 7;
localisationParams.wiggle = 1;
localisationParams.wiggleDeg = 5;
localisationParams.lobeDist = sqrt(dx^2 + dy^2);

%% Prep

directory_out = fullfile(directory_out,'results');
if ~exist(directory_out,'dir'); mkdir(directory_out); end

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

Scope = Microscope(NA,magnification,nImmersion,camOffset,camGain,camQE,camPixSize);
Proc = Process(Scope,bkgndEstimation,localisationParams);
warning('off')
w = 3;
results = {};
for id_roi=1:numRois
    fprintf('Processing roi %d/%d\n',id_roi,numRois)

    results(id_roi).roiID = id_roi;
    
    % perform localisation on the MIP of the 405/405 nm dataset
    filepath = fullfile(directory_tif,roiPaths(id_roi).filenameCalib);
    mip = Proc.getMaximumIntensityProjectionFromPath(filepath);
    locs = Proc.processFrame(mip,0);
    numLocs = length(locs.upperLobe);
    x_upper = [locs.upperLobe.x]/1e9;
    y_upper = [locs.upperLobe.y]/1e9;
    x_lower = [locs.lowerLobe.x]/1e9;
    y_lower = [locs.lowerLobe.y]/1e9;

    % 405 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename405);
    bkgnd = Proc.estimateBackground(filepath);
    figure;
    subplot(1,3,1); histogram(bkgnd(:));
    stack = Utils.readTiffStack(filepath);
    subplot(1,3,2); histogram(stack(:));
    stack = stack - bkgnd;
    subplot(1,3,3); histogram(stack(:));
    
    [intLobe1,intLobe2] = measureIntensities(x_lower,y_lower,x_upper,y_upper,stack,numLocs,w);
    results(id_roi).ratio405 = intLobe2./intLobe1;
    results(id_roi).intLobe1_405 = intLobe1;
    results(id_roi).intLobe2_405 = intLobe2;

    % 488 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename488);
    bkgnd = Proc.estimateBackground(filepath);
    stack = Utils.readTiffStack(filepath);
    stack = stack - bkgnd;
    [intLobe1,intLobe2] = measureIntensities(x_lower,y_lower,x_upper,y_upper,stack,numLocs,w);
    results(id_roi).ratio488 = intLobe2./intLobe1;
    results(id_roi).intLobe1_488 = intLobe1;
    results(id_roi).intLobe2_488 = intLobe2;

    % 561 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename561);
    bkgnd = Proc.estimateBackground(filepath);
    stack = Utils.readTiffStack(filepath);
    stack = stack - bkgnd;
    [intLobe1,intLobe2] = measureIntensities(x_lower,y_lower,x_upper,y_upper,stack,numLocs,w);
    results(id_roi).ratio561 = intLobe2./intLobe1;
    results(id_roi).intLobe1_561 = intLobe1;
    results(id_roi).intLobe2_561 = intLobe2;

    % 638 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename638);
    bkgnd = Proc.estimateBackground(filepath);
    stack = Utils.readTiffStack(filepath);
    stack = stack - bkgnd;
    [intLobe1,intLobe2] = measureIntensities(x_lower,y_lower,x_upper,y_upper,stack,numLocs,w);
    results(id_roi).ratio638 = intLobe2./intLobe1;
    results(id_roi).intLobe1_638 = intLobe1;
    results(id_roi).intLobe2_638 = intLobe2;
    
    pause(0.001)
end
warning('on')


%% Plot of median spectra for each quantum dot

intensityThreshold = 20e4;

wavelenghts = [405 488 561 638];
factorPower = [1 1 1 1];
correctionFactor = factorPower.*(405./wavelenghts);

meanRatio405 = [];
meanRatio488 = [];
meanRatio561 = [];
meanRatio638 = [];

numQDots = 0;
fig = figure('position',[50 50 350 300]);
set(0,'DefaultAxesTitleFontWeight','normal');

for id_roi = 1:numRois
    numLocs = size(results(id_roi).ratio405,2);
    for i=1:numLocs
        
        keep = logical(results(id_roi).intLobe2_405(:,i) > intensityThreshold);
        if sum(keep(:))
            ratio1 = results(id_roi).ratio405(:,i);
            ratio1 = nanmedian(ratio1(keep));
            meanRatio405 = [meanRatio405 ratio1];
    
            keep = logical(results(id_roi).intLobe2_488(:,i) > intensityThreshold);
            if sum(keep(:))
                ratio2 = results(id_roi).ratio488(:,i);
                ratio2 = nanmedian(ratio2(keep));
                meanRatio488 = [meanRatio488 ratio2];
        
                keep = logical(results(id_roi).intLobe2_561(:,i) > intensityThreshold);
                if sum(keep(:))
                    ratio3 = results(id_roi).ratio561(:,i);
                    ratio3 = nanmedian(ratio3(keep));
                    meanRatio561 = [meanRatio561 ratio3];
            
                    keep = logical(results(id_roi).intLobe2_638(:,i) > intensityThreshold);
                    if sum(keep(:))
                        ratio4 = results(id_roi).ratio638(:,i);
                        ratio4 = nanmedian(ratio4(keep));
                        meanRatio638 = [meanRatio638 ratio4];
                
                        ratios = [ratio1 ratio2 ratio3 ratio4];
                        ratios = correctionFactor./ratios;
                        plot(wavelenghts,ratios,'-','color',0.8*[1 1 1 0.4])
                        hold on
                        
                        numQDots = numQDots + 1;
                    end
                end
            end
        end
    end
end

spectrumBulk = readmatrix(pathBulkSpectrum);
x_lambda = spectrumBulk(:,1);
y_ex = spectrumBulk(:,2);
y_em = spectrumBulk(:,3);

absAt405 = y_ex(x_lambda == 405);

% plot(x_lambda,y_em/max(y_em(:))); hold on

% col = [0.1333,0.5289,0.8000]; % blue
col = [1.0000,0.5644,0.0622]; % orange
% col = [0.1725,0.6275,0.1725]; % green
% col = 1.05*[0.8392,0.1529,0.1569]; % red
% col = [0.5804,0.4039,0.7412];
% col = [0.5490,0.3373,0.2941];
% col = [0.8902,0.4667,0.7608];
plot(x_lambda,y_ex/absAt405,'-','Color',col,'linewidth',2); hold on

xlim([350 700])
ylim([0,1.2])
xticks(200:50:800)
hold on
meanRatios = [nanmedian(meanRatio405) nanmedian(meanRatio488) nanmedian(meanRatio561) nanmedian(meanRatio638)];
meanRatios = correctionFactor./meanRatios;
plot(wavelenghts,meanRatios,'k-o','linewidth',1)

grid on
xlabel('Excitation wavelength (nm)');
ylabel('Ratio intensity lobe2/lobe1');
title([num2str(numQDots) ' quantum dots, ' num2str(numRois) ' rois'])
set(gca,'FontSize',9)

savefig(fig,fullfile(directory_out,'estimatedExcitationSpectra.fig'))
exportgraphics(fig,fullfile(directory_out,'estimatedExcitationSpectra.png'),'Resolution',400)
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,'estimatedExcitationSpectra.eps'))

%% Functions

function [I_upperLobe,I_lowerLobe] = measureIntensities(x_lower,y_lower,x_upper,y_upper,stack,numLocs,w)
I_upperLobe = zeros(size(stack,3),numLocs);
I_lowerLobe = zeros(size(stack,3),numLocs);
for n=1:numLocs
    lobe1 = stack(round(x_lower(n))-w:round(x_lower(n))+w,round(y_lower(n))-w:round(y_lower(n))+w,:);
    lobe2 = stack(round(x_upper(n))-w:round(x_upper(n))+w,round(y_upper(n))-w:round(y_upper(n))+w,:);
    intLobe1 = sum(lobe1,1:2); I_upperLobe(:,n) = intLobe1(:);
    intLobe2 = sum(lobe2,1:2); I_lowerLobe(:,n) = intLobe2(:);
end
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

