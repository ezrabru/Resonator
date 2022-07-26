%% Batch semi-automated analysis of 2-lobe resonator QDot data

clear all
close all
clc

numRois = 15;
directory_tif = 'D:\manuscripts\resonator\new data\20201024_resonator quantitative\QD585'; % path to folder with tif stacks
directory_out = directory_tif;

flagCalib = 'B405_T405';

flag405 = 'B405_T405';
flag488 = 'B405_T488';
flag561 = 'B405_T561';
flag638 = 'B405_T638';

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

results = {};

for id_roi=1:numRois
    
    results(id_roi).roiID = id_roi;
    
    %% THIS STEP REQUIRES USER INPUT: selection of templates for localisation
    pathCalibStack = fullfile(directory_tif,roiPaths(id_roi).filenameCalib);
    [psf,MIP] = getMIPandPSFcutOut(pathCalibStack);
    
    %% THIS STEP REQUIRES USER INPUT: Accepting/rejecting localizations
    % Localizes using cross-correlation with selected experimental PSFs.
    % Reopen the files one by one and let user manually accept or reject
    % localisations by typing 'y' or 'n' in command window.
    threshold = 0.7;
    coordinates = localizeCrossCorr(MIP,psf,threshold);
    
    %% THIS STEP REQUIRES USER INPUT: determine distance between lobes
    [dx,dy] = getRelativeCoordinatesLobes(psf);
    
    %% Measure intensities
    windowIntensity = 7; w = ceil(windowIntensity/2);
    medfiltKernel = 30;
    
    % 405 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename405);
    [intLobe1,intLobe2] = measureIntensitiesLobes(filepath,w,dx,dy,coordinates);
    results(id_roi).ratio405 = intLobe1./intLobe2;
    results(id_roi).intLobe1_405 = intLobe1;
    results(id_roi).intLobe2_405 = intLobe2;
    
    % 488 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename488);
    [intLobe1,intLobe2] = measureIntensitiesLobes(filepath,w,dx,dy,coordinates);
    results(id_roi).ratio488 = intLobe1./intLobe2;
    results(id_roi).intLobe1_488 = intLobe1;
    results(id_roi).intLobe2_488 = intLobe2;
    
    % 561 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename561);
    [intLobe1,intLobe2] = measureIntensitiesLobes(filepath,w,dx,dy,coordinates);
    results(id_roi).ratio561 = intLobe1./intLobe2;
    results(id_roi).intLobe1_561 = intLobe1;
    results(id_roi).intLobe2_561 = intLobe2;
    
    % 638 nm
    filepath = fullfile(directory_tif,roiPaths(id_roi).filename638);
    [intLobe1,intLobe2] = measureIntensitiesLobes(filepath,w,dx,dy,coordinates);
    results(id_roi).ratio638 = intLobe1./intLobe2;
    results(id_roi).intLobe1_638 = intLobe1;
    results(id_roi).intLobe2_638 = intLobe2;
    
    close all
end

%% Plot of median spectra for each quantum dot

clear all
load('D:\manuscripts\resonator\figures\figure2\quantitative resonator\data QD585\workspace.mat');
outputdirfig = 'D:\manuscripts\resonator\figures\figure2\quantitative resonator\data QD585';

numRois = 14;

intensityThreshold = 10e4;

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
        
        numQDots = numQDots + 1;
        
        keep = logical(results(id_roi).intLobe1_405(:,i) > intensityThreshold);
        ratio1 = results(id_roi).ratio405(:,i);
        ratio1 = nanmedian(ratio1(keep));
        meanRatio405 = [meanRatio405 ratio1];

        keep = logical(results(id_roi).intLobe1_488(:,i) > intensityThreshold);
        ratio2 = results(id_roi).ratio488(:,i);
        ratio2 = nanmedian(ratio2(keep));
        meanRatio488 = [meanRatio488 ratio2];

        keep = logical(results(id_roi).intLobe1_561(:,i) > intensityThreshold);
        ratio3 = results(id_roi).ratio561(:,i);
        ratio3 = nanmedian(ratio3(keep));
        meanRatio561 = [meanRatio561 ratio3];

        keep = logical(results(id_roi).intLobe1_638(:,i) > intensityThreshold);
        ratio4 = results(id_roi).ratio638(:,i);
        ratio4 = nanmedian(ratio4(keep));
        meanRatio638 = [meanRatio638 ratio4];

        ratios = [ratio1 ratio2 ratio3 ratio4];
        ratios = correctionFactor./ratios;
        plot(wavelenghts,ratios,'-','color',0.8*[1 1 1 0.4])
        hold on
    end
end

pathBulkSpectrum = 'D:\manuscripts\resonator\figures\figure2\quantitative resonator\data spectra\dyes\Qdot585';
spectrumBulk = readmatrix(pathBulkSpectrum);
x_lambda = spectrumBulk(:,1);
y_ex = spectrumBulk(:,2);
y_em = spectrumBulk(:,3);

absAt405 = y_ex(x_lambda == 405);

% plot(x_lambda,y_em); hold on

% col = [0.1333,0.5289,0.8000]; % blue
col = [1.0000,0.5644,0.0622]; % orange
% col = [0.1725,0.6275,0.1725]; % green
% col = 1.05*[0.8392,0.1529,0.1569]; % red
% col = [0.5804,0.4039,0.7412];
% col = [0.5490,0.3373,0.2941];
% col = [0.8902,0.4667,0.7608];
plot(x_lambda,y_ex/absAt405,'-','Color',col,'linewidth',2); hold on

xlim([375 670])
ylim([-0.1,1.5])
hold on
meanRatios = [nanmean(meanRatio405) nanmean(meanRatio488) nanmean(meanRatio561) nanmean(meanRatio638)];
meanRatios = correctionFactor./meanRatios;
plot(wavelenghts,meanRatios,'k-o','linewidth',1)

grid on
xlabel('Wavelength (nm)');
ylabel('Ratio intensity lobe2/lobe1');
title([num2str(numQDots) ' quantum dots, ' num2str(numRois) ' rois'])
set(gca,'FontSize',11)
savefig(fig,fullfile(outputdirfig,'estimatedExcitationSpectra.fig'))
saveas(fig,fullfile(outputdirfig,'estimatedExcitationSpectra.svg'))
saveas(fig,fullfile(outputdirfig,'estimatedExcitationSpectra.png'))

%% Comparison width distribution ratios

clc

correctionPulse = (1/1.03);

ensembleStdRatio405 = [];
ensembleStdRatio488 = [];
ensembleStdRatio561 = [];
ensembleStdRatio638 = [];

singleStdRatio405 = [];
singleStdRatio488 = [];
singleStdRatio561 = [];
singleStdRatio638 = [];

stdRatio405 = [];
stdRatio488 = [];
stdRatio561 = [];
stdRatio638 = [];

numQDots = 0;

intensityThreshold = 2e5;

for id_roi = 1:numRois

    numLocs = size(results(id_roi).ratio405,2);

    for i=1:numLocs
        
        numQDots = numQDots + 1;
        
        keep = logical(results(id_roi).intLobe1_405(:,i) > intensityThreshold);
        ratio1 = correctionPulse*(405/405)./results(id_roi).ratio405(:,i); ratio1 = ratio1(keep);
        stdRatio405 = [stdRatio405 std(ratio1(:))];
        ensembleStdRatio405 = [ensembleStdRatio405 ratio1(:)'];
        singleStdRatio405 = [singleStdRatio405 ratio1(:)'-median(ratio1(:))];
        
        keep = logical(results(id_roi).intLobe1_488(:,i) > intensityThreshold);
        ratio2 = correctionPulse*(405/488)./results(id_roi).ratio488(:,i); ratio2 = ratio2(keep);
        stdRatio488 = [stdRatio488 std(ratio2(:))];
        ensembleStdRatio488 = [ensembleStdRatio488 ratio2(:)'];
        singleStdRatio488 = [singleStdRatio488 ratio2(:)'-median(ratio2(:))];

        keep = logical(results(id_roi).intLobe1_561(:,i) > intensityThreshold);
        ratio3 = correctionPulse*(405/561)./results(id_roi).ratio561(:,i); ratio3 = ratio3(keep);
        stdRatio561 = [stdRatio561 std(ratio3(:))];
        ensembleStdRatio561 = [ensembleStdRatio561 ratio3(:)'];
        singleStdRatio561 = [singleStdRatio561 ratio3(:)'-median(ratio3(:))];

        keep = logical(results(id_roi).intLobe1_638(:,i) > intensityThreshold);
        ratio4 = correctionPulse*(405/638)./results(id_roi).ratio638(:,i); ratio4 = ratio4(keep);
        stdRatio638 = [stdRatio638 std(ratio4(:))];
        ensembleStdRatio638 = [ensembleStdRatio638 ratio4(:)'];
        singleStdRatio638 = [singleStdRatio638 ratio4(:)'-median(ratio4(:))];
        
    end
end


fig = figure;
bins = 0:0.01:2;
histogram(ensembleStdRatio405(:),bins,'edgecolor','none'); hold on
histogram(ensembleStdRatio488(:),bins,'edgecolor','none')
histogram(ensembleStdRatio561(:),bins,'edgecolor','none')
histogram(ensembleStdRatio638(:),bins,'edgecolor','none')
xlabel('Ratio'); ylabel('Occurence'); title('Ensemble')
legend('405/405 nm','488/405 nm','561/405 nm','638/405 nm')

%%

fig = figure('position',[50 50 800 200]);
bins = -0.25:0.01:0.25;
subplot(1,4,1);
histogram(ensembleStdRatio405(:)-median(ensembleStdRatio405(:)),bins,'edgecolor','none','facecolor',[0 0 0],'facealpha',0.3); hold on
histogram(singleStdRatio405(:)-median(singleStdRatio405(:)),bins,'edgecolor','none','facecolor',[0.5804,0.4039,0.7412]);
set(gca,'fontsize',fontsize)
title('405/405 nm'); xlabel('ratio - \langle ratio \rangle'); ylabel('Occurence')

subplot(1,4,2);
histogram(ensembleStdRatio488(:)-median(ensembleStdRatio488(:)),bins,'edgecolor','none','facecolor',[0 0 0],'facealpha',0.3); hold on
histogram(singleStdRatio488(:)-median(singleStdRatio488(:)),bins,'edgecolor','none','facecolor',[0.1333,0.5289,0.8000])
set(gca,'fontsize',fontsize)
title('488/405 nm'); xlabel('ratio - \langle ratio \rangle'); ylabel('Occurence')

subplot(1,4,3);
histogram(ensembleStdRatio561(:)-median(ensembleStdRatio561(:)),bins,'edgecolor','none','facecolor',[0 0 0],'facealpha',0.3); hold on
histogram(singleStdRatio561(:)-median(singleStdRatio561(:)),bins,'edgecolor','none','facecolor',[0.7,0.6,0.0622])
set(gca,'fontsize',fontsize)
title('561/405 nm'); xlabel('ratio - \langle ratio \rangle'); ylabel('Occurence')

subplot(1,4,4);
histogram(ensembleStdRatio638(:)-median(ensembleStdRatio638(:)),bins,'edgecolor','none','facecolor',[0 0 0],'facealpha',0.3); hold on
histogram(singleStdRatio638(:)-median(singleStdRatio638(:)),bins,'edgecolor','none','facecolor',[0.8392,0.1529,0.1569])
set(gca,'fontsize',fontsize)
title('638/405 nm'); xlabel('ratio - \langle ratio \rangle'); ylabel('Occurence')

%% Comparison width distribution ratios

fontsize = 11;

fig = figure('position',[100 400 500 220]);
set(0,'DefaultAxesTitleFontWeight','normal');

x = [1 2];
ylims = [0 0.25];

subplot(1,4,1);
y = [iqr(ensembleStdRatio405(:)) nanmedian(stdRatio405)];
err = [0 iqr(stdRatio405,'all')];
bar(x,y); hold on; errorbar(x(2),y(2),err(2),'k'); ylim(ylims)
ylabel('\sigma_{ratio}');
title('405/405 nm')
set(gca,'xticklabel',{'ensemble','single'}); xtickangle(45)
set(gca,'fontsize',fontsize)

subplot(1,4,2);
y = [iqr(ensembleStdRatio488(:)) nanmedian(stdRatio488)];
err = [0 iqr(stdRatio488,'all')];
bar(x,y); hold on; errorbar(x(2),y(2),err(2),'k'); ylim(ylims)
% ylabel('\sigma_{ratio}');
title('488/405 nm')
set(gca,'xticklabel',{'ensemble','single'}); xtickangle(45)
set(gca,'fontsize',fontsize)

subplot(1,4,3);
y = [iqr(ensembleStdRatio561(:)) nanmedian(stdRatio561)];
err = [0 iqr(stdRatio561,'all')];
bar(x,y); hold on; errorbar(x(2),y(2),err(2),'k'); ylim(ylims)
% ylabel('\sigma_{ratio}');
title('561/405 nm')
set(gca,'xticklabel',{'ensemble','single'}); xtickangle(45)
set(gca,'fontsize',fontsize)

subplot(1,4,4);
y = [iqr(ensembleStdRatio638(:)) nanmedian(stdRatio638)];
err = [0 iqr(stdRatio638,'all')];
bar(x,y); hold on; errorbar(x(2),y(2),err(2),'k'); ylim(ylims)
% ylabel('\sigma_{ratio}');
title('638/405 nm')
set(gca,'xticklabel',{'ensemble','single'}); xtickangle(45)
set(gca,'fontsize',fontsize)

savefig(fig,fullfile(outputdirfig,'comparisonWidthDistributionRatios.fig'))
saveas(fig,fullfile(outputdirfig,'comparisonWidthDistributionRatios.svg'))
saveas(fig,fullfile(outputdirfig,'comparisonWidthDistributionRatios.png'))

%%

% numQDots = 0;
% fig = figure('position',[50 50 350 300]);
% set(0,'DefaultAxesTitleFontWeight','normal');
% 
% for id_roi = 1:numRois
% 
%     numLocs = size(results(id_roi).ratio405,2);
% 
%     for i=1:numLocs
%         
%         numQDots = numQDots + 1;
%         
%         keep = logical(results(id_roi).intLobe1_405(:,i) > intensityThreshold);
%         ratio1 = results(id_roi).ratio405(:,i);
%         ratio1 = nanmedian(ratio1(keep));
%         meanRatio405 = [meanRatio405 ratio1];
% 
%         keep = logical(results(id_roi).intLobe1_488(:,i) > intensityThreshold);
%         ratio2 = results(id_roi).ratio488(:,i);
%         ratio2 = nanmedian(ratio2(keep));
%         meanRatio488 = [meanRatio488 ratio2];
% 
%         keep = logical(results(id_roi).intLobe1_561(:,i) > intensityThreshold);
%         ratio3 = results(id_roi).ratio561(:,i);
%         ratio3 = nanmedian(ratio3(keep));
%         meanRatio561 = [meanRatio561 ratio3];
% 
%         keep = logical(results(id_roi).intLobe1_638(:,i) > intensityThreshold);
%         ratio4 = results(id_roi).ratio638(:,i);
%         ratio4 = nanmedian(ratio4(keep));
%         meanRatio638 = [meanRatio638 ratio4];
% 
%         scatter((405/405)./ratio1,(405/488)./ratio2,30,'.'); xlim([0 1.5]); ylim([0 1])
% %         scatter((405/405)./ratio1,(405/561)./ratio3,30,'.'); xlim([0 1.5]); ylim([0 1])
% %         scatter((405/405)./ratio1,(405/638)./ratio4,30,'.'); xlim([0 1.5]); ylim([0 1])
% %         scatter((405/488)./ratio2,(405/561)./ratio3,30,'.'); xlim([0 1.5]); ylim([0 1])
% %         scatter((405/488)./ratio2,(405/638)./ratio4,30,'.'); xlim([0 1.5]); ylim([0 1])
%         hold on
%     end
% end

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