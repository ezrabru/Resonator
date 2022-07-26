%% Batch semi-automated analysis of 2-lobe resonator QDot data

clear all
close all
clc
addpath('lib');

directory_tif  = fullfile(pwd,'testdata'); % path to folder with tif stacks
directory_out  = directory_tif;

pathCalibration = fullfile(pwd,'testdata','QD655_roi1_B405_T40500000.tif');

%% Read data and get max intensity projections

% Create new output folder
outputdir = fullfile(directory_out,'results'); if ~exist(outputdir, 'dir'); mkdir(outputdir); end
dir_MIP   = fullfile(outputdir,'MIPs'); if ~exist(dir_MIP, 'dir'); mkdir(dir_MIP); end

% Get list of tif files in folder
filelist = dir(fullfile(directory_tif,'*.tif'));
disp(['Number of tif files found: ' num2str(length(filelist))])
for i = 1:length(filelist); disp("  "+filelist(i).name); end; disp(' ')

% Loop over files and generate maximum intensity projections
disp('Generating maximum intensity projections...')
for i = 1:length(filelist)
    disp(['  ' num2str(i) '/' num2str(length(filelist))])
    stack = readTifStack(fullfile(directory_tif,filelist(i).name));
    MIP = max(stack,[],3);
    writeTifImage(MIP,fullfile(dir_MIP,strcat(extractBefore(filelist(i).name,'.'),'_MIP.tif')))
end
clear variable stack

%% THIS STEP REQUIRES USER INPUT: selection of templates for localisation

% Let user select an example PSF from data for localization
results = struct;
for i = 1:length(filelist)
    % Let user pick rois
    MIP = imread(fullfile(dir_MIP,strcat(extractBefore(filelist(i).name,'.'),'_MIP.tif')));
    fig = figure('Name',['  ' num2str(i) '/' num2str(length(filelist))]);
    imshow(MIP,[]); set(gcf,'Position',get(0,'Screensize'));
    roi = drawrectangle;
    x = floor(roi.Position(1)); y = floor(roi.Position(2));
    w = ceil(roi.Position(3));  h = ceil(roi.Position(4));
    roi = MIP(y:y+h,x:x+w);
    imshow(roi,[]); pause(1);
    results(i).filename = filelist(i).name;
    results(i).psfCalibration = roi;
end

%% THIS STEP REQUIRES USER INPUT: Accepting/rejecting localizations
% Localizes using cross-correlation with selected experimental PSFs.
% Reopen the files one by one and let user manually accept or reject
% localisations by typing 'y' or 'n' in command window.

treshold = 0.7;

for i = 1:length(filelist)
    psf = results(i).psfCalibration;
    MIP = imread(fullfile(dir_MIP,strcat(extractBefore(filelist(i).name,'.'),'_MIP.tif')));
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
    results(i).x = x(logical(keep));
    results(i).y = y(logical(keep));
    
    x_rejected = x(~logical(keep));
    y_rejected = y(~logical(keep));
    
    % Annotate MIP with chosen rois
    x = results(i).x;
    y = results(i).y;
    fig = figure;
    imshow(MIP,[]); hold on
    for j = 1:length(x)
        rectangle('Position',[x(j)-size(psf,2)/2,y(j)-size(psf,1)/2,size(psf,2),size(psf,1)],'EdgeColor','y');
        text(x(j)-size(psf,2)/2,y(j),string(j),'Color','y');
    end
    saveas(fig,fullfile(dir_MIP,strcat(extractBefore(filelist(i).name,'.'),'_MIP_annotated.png')))
    close(fig);
    
end

%% THIS STEP REQUIRES USER INPUT: determine distance between lobes

calib = readTifStack(pathCalibration);
MIP = max(calib,[],3);

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
    
    x = results(i).x;
    y = results(i).y;
    numLocs = length(x);
    
    stack = readTifStack(fullfile(directory_tif,results(i).filename));
    stack = stack - results(i).bkgnd;
    stack(stack < 0) = 0;
    
    intensitiesLobe1 = zeros(size(stack,3),numLocs);
    intensitiesLobe2 = zeros(size(stack,3),numLocs);
    for j=1:numLocs
        lobe1 = stack(round(y(j)-dy/2)-w:round(y(j)-dy/2)+w,round(x(j)-dx/2)-w:round(x(j)-dx/2)+w,:);
        lobe2 = stack(round(y(j)+dy/2)-w:round(y(j)+dy/2)+w,round(x(j)+dx/2)-w:round(x(j)+dx/2)+w,:);
        
        intLobe1 = sum(lobe1,1:2); intensitiesLobe1(:,j) = intLobe1(:);
        intLobe2 = sum(lobe2,1:2); intensitiesLobe2(:,j) = intLobe2(:);
        
%         subplot(2,1,1); imshow(lobe1(:,:,1),[]);
%         subplot(2,1,2); imshow(lobe2(:,:,1),[]); pause(0.1);
    end

    results(i).intLobe1 = intensitiesLobe1;
    results(i).intLobe2 = intensitiesLobe2;
    results(i).ratio12 = intensitiesLobe1./intensitiesLobe2;
    results(i).n = length(intensitiesLobe1(:));
    
end

save('results.mat','results')

%% Plotting

figure;
for i = 1:length(filelist)
    xVal = results(i).intLobe1;
    yVal = results(i).intLobe2;
    plot(xVal(:),yVal(:),'x')
    hold on
end
set(gca,'FontSize',12)
legend(filelist.name)
xlabel('Intensity lobe 1');
ylabel('Intensity lobe 2');
set(gca,'FontSize',14)


wavelengths = [405 488 561 638];
factorPower = [1 0.88 0.88 0.88];
figure;
for i = 1:length(filelist)
    plot(wavelengths(i),factorPower(i)*(405/wavelengths(i))*median(median(1./results(i).ratio12)),'+k','MarkerSize',10)
    hold on
end
plot(spectra1,spectra2/37.85,'LineWidth',2)
plot(spectra1,spectra3/100,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('Wavelength (nm)');
ylabel('Ratio intensity lobe1/lobe2');
set(gca,'FontSize',14)



figure;
for i = 1:length(filelist)
    histogram(results(i).ratio12)
    hold on
end
set(gca,'FontSize',12)
xlabel('Ratio intensity lobe1/lobe2');
ylabel('Occurence');
set(gca,'FontSize',14)