%% Batch analysis of 2-lobe resonator QDot data

clear all; close all; clc
addpath('lib')

directories_tif = {'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi1',...
                   'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi2',...
                   'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi3',...
                   'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi4',...
                   'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi5',...
                   'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi6'}; % path to folder with tif stacks

%%

for id_roi = 1:numel(directories_tif)
    
    directory_tif = directories_tif{id_roi};
    directory_out = directory_tif;
    
    % PSF displacement
    dx = 0;
    dy = 37;
    
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
    localisationParams.DoGminThreshold = 200;
    localisationParams.method = 'asymmetric gaussian';
    localisationParams.w = 7;
    localisationParams.wiggle = 4;
    localisationParams.wiggleDeg = 5;
    localisationParams.lobeDist = sqrt(dx^2 + dy^2);
    localisationParams.dx = dx;
    localisationParams.dy = dy;
    
    
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
    
    % get localisations
    filepath = fullfile(filelist(1).folder,filelist(1).name);
    avg = Proc.getMaximumIntensityProjectionFromPath(filepath);

    locsGrouped = Proc.processFrame_2lobe(avg,0);
    locs_upper = locsGrouped.upperLobe;
    locs_lower = locsGrouped.lowerLobe;
    
    %%

    col1 = 1.05*[0.8392,0.1529,0.1569]; % red
    col2 = [0.1333,0.5289,0.8000]; % blue
    
    ms = 80;
    
    fig = figure;
    imshow(avg,[]); hold on
    counter = 0;
    for i = 1:height(locs_upper)
        counter = counter + 1;
        scatter(locs_upper.y(i,:)/1e9,locs_upper.x(i,:)/1e9,ms,col1); hold on
        scatter(locs_lower.y(i,:)/1e9,locs_lower.x(i,:)/1e9,ms,col2);
        x_text = (locs_upper.y(i,:)/1e9 + locs_lower.y(i,:)/1e9)/2; 
        y_text = (locs_upper.x(i,:)/1e9 + locs_lower.x(i,:)/1e9)/2; 
        text(x_text,y_text,sprintf('%d',counter),'Color','w')
        pause(0.01)
    end
    legend('upper lobe','lower lobe','Color','none','EdgeColor','w','TextColor','w')
    [~,filename,~] = fileparts(filepath);
    title(filename,'Interpreter','none','FontWeight','normal')
    
    % save annotated figure light mode
    savefig(fig,fullfile(directory_out,'annotated_mip_lightMode.fig'))
    exportgraphics(fig,fullfile(directory_out,'annotated_mip_lightMode.png'),'Resolution',400)
    set(gcf,'renderer','Painters')
    exportgraphics(fig,fullfile(directory_out,'annotated_mip_lightMode.eps'))
    
    % save annotated figure dark mode
    set(gcf,'Color','k')
    title(filename,'Interpreter','none','FontWeight','normal','Color','w')
    savefig(fig,fullfile(directory_out,'annotated_mip_darkMode.fig'))
    exportgraphics(fig,fullfile(directory_out,'annotated_mip_darkMode.png'),'Resolution',400,'BackgroundColor','k')
    set(gcf,'renderer','Painters')
    exportgraphics(fig,fullfile(directory_out,'annotated_mip_darkMode.eps'))
    
    
    %% Loop over all files from experiment and measure intensities of the lobes
    
    results = {};
    tic;
    
    numFiles = numel(filelist);
    for id_file=1:numFiles
        fprintf('Processing file %d/%d\n',id_file,numFiles)
    
        filepath = fullfile(filelist(id_file).folder,filelist(id_file).name);
    
        % estimate background
        t1 = toc; fprintf('  Estimating background... ')
        bkgnd = Proc.estimateBackground(filepath);
        fprintf('done! (%.2f s)\n',toc-t1)
        
        % read stack
        t1 = toc; fprintf('  Reading stack........... ')
        stack = Utils.readTiffStack(filepath);
        stack = stack - bkgnd;
        fprintf('done! (%.2f s)\n',toc-t1)
    
        % measure intensities of the lobes
        x_lower = [locsGrouped.lowerLobe.x]/1e9;
        y_lower = [locsGrouped.lowerLobe.y]/1e9;
        x_upper = [locsGrouped.upperLobe.x]/1e9;
        y_upper = [locsGrouped.upperLobe.y]/1e9;
        t1 = toc; fprintf('  Measuring intensities... ')
        [I_upperLobe,I_lowerLobe] = measureIntensities(x_lower,y_lower,x_upper,y_upper,stack,w);
        fprintf('done! (%.2f s)\n',toc-t1)
    
        % append these to the results section
        results(id_file).roi_id = id_file;
        results(id_file).filepath = filepath;
        results(id_file).I_upperLobe = I_upperLobe;
        results(id_file).I_lowerLobe = I_lowerLobe;
        
        clear stack
    end
    
    save(fullfile(directory_out,'results.mat'),'results')
    
    %% Reshape results to be organised per quantum dot and traces merged
    
    results_curated = struct;
    results_curated.mip = avg;
    results_curated.upperLobe = table2struct(locsGrouped.upperLobe);
    results_curated.lowerLobe = table2struct(locsGrouped.lowerLobe);
    
    numQdots = height(locsGrouped.upperLobe);
    for id_qdot=1:numQdots
        concatenatedTrace_upper = [];
        concatenatedTrace_lower = [];
        for id_file=1:numFiles
            concatenatedTrace_upper = [concatenatedTrace_upper; results(id_file).I_upperLobe(:,id_qdot)];
            concatenatedTrace_lower = [concatenatedTrace_lower; results(id_file).I_lowerLobe(:,id_qdot)];
        end
        results_curated.upperLobe(id_qdot).intensityTrace = concatenatedTrace_upper;
        results_curated.lowerLobe(id_qdot).intensityTrace = concatenatedTrace_lower;
        
        totalIntensity = concatenatedTrace_upper + concatenatedTrace_lower;
        results_curated.ratios(id_qdot).totalIntensity = totalIntensity;
        results_curated.ratios(id_qdot).fractionUpperTotalInt = concatenatedTrace_upper./totalIntensity;
        results_curated.ratios(id_qdot).fractionLowerTotalInt = concatenatedTrace_lower./totalIntensity;
    end
    
    save(fullfile(directory_out,'results_curated.mat'),'results_curated')
    
    close all
end

%% Functions

function [I_upperLobe,I_lowerLobe] = measureIntensities(x_lower,y_lower,x_upper,y_upper,stack,w)
numLocs = numel(x_lower);
I_upperLobe = zeros(size(stack,3),numLocs);
I_lowerLobe = zeros(size(stack,3),numLocs);
for n=1:numLocs
    lobe1 = stack(round(x_lower(n))-w:round(x_lower(n))+w,round(y_lower(n))-w:round(y_lower(n))+w,:);
    lobe2 = stack(round(x_upper(n))-w:round(x_upper(n))+w,round(y_upper(n))-w:round(y_upper(n))+w,:);
    intLobe1 = sum(lobe1,1:2); I_upperLobe(:,n) = intLobe1(:);
    intLobe2 = sum(lobe2,1:2); I_lowerLobe(:,n) = intLobe2(:);
end
end
