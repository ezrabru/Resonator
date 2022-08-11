%% Batch analysis of 3-lobe resonator QDot data

clear all; close all; clc
addpath('lib')

% PSF displacement
dx = 0;
dy = 37;

dx3 = 28;
dy3 = 16.5;

directory_tif = 'D:\manuscripts\resonator\new data\20201013_resonator 3 lobe\split\roi8'; % path to folder with tif stacks
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
localisationParams.DoGminThreshold = 200;
localisationParams.method = 'asymmetric gaussian';
localisationParams.w = 7;
localisationParams.wiggle = 4;
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

% get localisations
filepath = fullfile(filelist(1).folder,filelist(1).name);
avg = Proc.getMaximumIntensityProjectionFromPath(filepath);
locs = Proc.processFrame(avg,0);
locs = struct2table(locs);

% split based on sigma into resonator and third lobes
threshold_sigmay = 1.7;
locs_resonator = locs(locs.sigmay > threshold_sigmay,:);
locs_dichroic = locs(locs.sigmay <= threshold_sigmay,:);

% group the resonator lobes based on distance
locs_resonator = table2struct(locs_resonator);
[locs_upper,locs_lower] = Proc.pairLocalisationsByDistance(locs_resonator);
% convert to tables
locs_upper = struct2table(locs_upper);
locs_lower = struct2table(locs_lower);

% loop over each upper lobe and check if there is a matching lower (and
% that it is not accidentally a lower lobe...)
locsGrouped.upperLobe = table;
locsGrouped.lowerLobe = table;
locsGrouped.thirdLobe = table;

%%

% remove any duplicate rows
locs_upper = unique(locs_upper,'stable');
numUpperLocs = height(locs_upper);

cols = lines(3);
col1 = cols(1,:);
col2 = cols(2,:);
col3 = cols(3,:);

ms = 80;

fig = figure;
imshow(avg,[]); hold on
counter = 0;
for i = 1:numUpperLocs
    x_i = locs_upper.x(i)/1e9;
    y_i = locs_upper.y(i)/1e9;

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
        if D_min_3l <= localisationParams.wiggle % third lobe was also found
            if isempty(locsGrouped.upperLobe)
                locsGrouped.upperLobe = locs_upper(i,:); 
                locsGrouped.lowerLobe = locs_lower(id_D_min_ll,:); 
                locsGrouped.thirdLobe = locs_dichroic(id_D_min_3l,:); 
            else
                locsGrouped.upperLobe = [locsGrouped.upperLobe; locs_upper(i,:)]; 
                locsGrouped.lowerLobe = [locsGrouped.lowerLobe; locs_lower(id_D_min_ll,:)]; 
                locsGrouped.thirdLobe = [locsGrouped.thirdLobe; locs_dichroic(id_D_min_3l,:)]; 
            end

            counter = counter + 1;
            
            scatter(y_i,x_i,ms,col1); hold on
            scatter(locs_lower.y(id_D_min_ll,:)/1e9,locs_lower.x(id_D_min_ll,:)/1e9,ms,col2);
            scatter(locs_dichroic.y(id_D_min_3l,:)/1e9,locs_dichroic.x(id_D_min_3l,:)/1e9,ms,col3);
            x_text = (y_i + locs_lower.y(id_D_min_ll,:)/1e9)/2; 
            y_text = (x_i + locs_lower.x(id_D_min_ll,:)/1e9)/2; 
            text(x_text,y_text,sprintf('%d',counter),'Color','w')
            pause(0.01)
        end
    end

end
legend('upper lobe','lower lobe','third lobe','Color','none','EdgeColor','w','TextColor','w')
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
    x_third = [locsGrouped.thirdLobe.x]/1e9;
    y_third = [locsGrouped.thirdLobe.y]/1e9;
    t1 = toc; fprintf('  Measuring intensities... ')
    [I_upperLobe,I_lowerLobe,I_thirdLobe] = measureIntensities(x_lower,y_lower,x_upper,y_upper,x_third,y_third,stack,w);
    fprintf('done! (%.2f s)\n',toc-t1)

    % append these to the results section
    results(id_file).roi_id = id_file;
    results(id_file).filepath = filepath;
    results(id_file).I_upperLobe = I_upperLobe;
    results(id_file).I_lowerLobe = I_lowerLobe;
    results(id_file).I_thirdLobe = I_thirdLobe;
    
    clear stack
end

save(fullfile(directory_out,'results.m'),'results')

%% Reshape results to be organised per quantum dot and traces merged

results_curated = struct;
results_curated.mip = avg;
results_curated.upperLobe = table2struct(locsGrouped.upperLobe);
results_curated.lowerLobe = table2struct(locsGrouped.lowerLobe);
results_curated.thirdLobe = table2struct(locsGrouped.thirdLobe);

numQdots = height(locsGrouped.upperLobe);
for id_qdot=1:numQdots
    concatenatedTrace_upper = [];
    concatenatedTrace_lower = [];
    concatenatedTrace_third = [];
    for id_file=1:numFiles
        concatenatedTrace_upper = [concatenatedTrace_upper; results(id_file).I_upperLobe(:,id_qdot)];
        concatenatedTrace_lower = [concatenatedTrace_lower; results(id_file).I_lowerLobe(:,id_qdot)];
        concatenatedTrace_third = [concatenatedTrace_third; results(id_file).I_thirdLobe(:,id_qdot)];
    end
    results_curated.upperLobe(id_qdot).intensityTrace = concatenatedTrace_upper;
    results_curated.lowerLobe(id_qdot).intensityTrace = concatenatedTrace_lower;
    results_curated.thirdLobe(id_qdot).intensityTrace = concatenatedTrace_third;
    
    totalIntensity = concatenatedTrace_upper + concatenatedTrace_lower + concatenatedTrace_third;
    results_curated.ratios(id_qdot).totalIntensity = totalIntensity;
    results_curated.ratios(id_qdot).fractionUpperTotalInt = concatenatedTrace_upper./totalIntensity;
    results_curated.ratios(id_qdot).fractionLowerTotalInt = concatenatedTrace_lower./totalIntensity;
    results_curated.ratios(id_qdot).fractionThirdTotalInt = concatenatedTrace_third./totalIntensity;

    totalIntensityResonator = concatenatedTrace_upper + concatenatedTrace_lower;
    results_curated.ratios(id_qdot).totalIntensityResonator = totalIntensityResonator;
    results_curated.ratios(id_qdot).fractionUpperOverTotalIntResonator = concatenatedTrace_upper./totalIntensityResonator;
    results_curated.ratios(id_qdot).fractionLowerOverTotalIntResonator = concatenatedTrace_lower./totalIntensityResonator;
    results_curated.ratios(id_qdot).ratioThirdOverResonator = concatenatedTrace_third./totalIntensityResonator;
    results_curated.ratios(id_qdot).ratioResonatorOverThird = totalIntensityResonator./concatenatedTrace_third;
end

save(fullfile(directory_out,'results_curated.m'),'results_curated')



%% plot intensity traces

% cols = lines(3);
% alpha = 0.8;
% 
% for i=1:size(I_upperLobe,2)
%     fig = figure('Position',[200,500,1080,200]);
%     plot(results_curated.upperLobe(i).intensityTrace + 20e5,'Color',[cols(1,:) alpha]); hold on
%     plot(results_curated.lowerLobe(i).intensityTrace + 10e5,'Color',[cols(2,:) alpha]);
%     plot(results_curated.thirdLobe(i).intensityTrace,'Color',[cols(3,:) alpha]); hold on
%     legend('upper','lower','third')
%     xlabel('Time (frames)'); ylabel('Intensity (adu)');
% 
%     % save figure light mode
%     savefig(fig,fullfile(directory_out,sprintf('traces_lightMode_id%d.fig',i)))
%     exportgraphics(fig,fullfile(directory_out,sprintf('traces_lightMode_id%d.png',i)),'Resolution',400,'BackgroundColor','w')
%     set(gcf,'renderer','Painters')
%     exportgraphics(fig,fullfile(directory_out,sprintf('traces_lightMode_id%d.eps',i)))
%     
%     % save annotated figure dark mode
%     set(gcf,'Color','k')
%     set(gca,'Color','k','xColor','w','yColor','w')
%     legend('Color','k','TextColor','w','EdgeColor','w')
%     title(filename,'Interpreter','none','FontWeight','normal','Color','w')
%     savefig(fig,fullfile(directory_out,sprintf('traces_darkMode_id%d.fig',i)))
%     exportgraphics(fig,fullfile(directory_out,sprintf('traces_darkMode_id%d.png',i)),'Resolution',400,'BackgroundColor','k')
%     set(gcf,'renderer','Painters')
%     exportgraphics(fig,fullfile(directory_out,sprintf('traces_darkMode_id%d.eps',i)),'BackgroundColor','k')
%     
%     close all
% end

%% plot scatter plots

% alpha = 0.2;
% min_total = -1e5; % minimum value used for axis limits
% 
% for i=1:size(I_upperLobe,2)
%     fig = figure('Position',[300,300,1300,350]);
%     
%     x = results_curated.upperLobe(i).intensityTrace;
%     y = results_curated.lowerLobe(i).intensityTrace;
%     z = results_curated.thirdLobe(i).intensityTrace;
%     c = 1:numel(x);
% 
%     max_res = max([x;y]);
%     max_third = max(z);
%     max_total = max([max_res max_third]);
% 
%     subplot(1,3,1)
%     scatter(y,x,10,c,'filled','MarkerFaceAlpha',alpha)
%     xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
%     xlabel('Intensity lower lobe'); ylabel('Intensity upper lobe');
%     colormap turbo; c1 = colorbar; c1.Label.String = 'Frame';
% 
%     subplot(1,3,2)
%     scatter(z,x,10,c,'filled','MarkerFaceAlpha',alpha)
%     xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
%     xlabel('Intensity third lobe'); ylabel('Intensity upper lobe');
%     colormap turbo; c2 = colorbar; c2.Label.String = 'Frame';
%     
%     subplot(1,3,3)
%     scatter(z,y,10,c,'filled','MarkerFaceAlpha',alpha)
%     xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
%     xlabel('Intensity third lobe'); ylabel('Intensity lower lobe');
%     colormap turbo; c3 = colorbar; c3.Label.String = 'Frame';
% 
%     % save figure light mode
%     sgtitle(filename,'Interpreter','none','FontWeight','normal')
%     savefig(fig,fullfile(directory_out,sprintf('scatter_lightMode_id%d.fig',i)))
%     exportgraphics(fig,fullfile(directory_out,sprintf('scatter_lightMode_id%d.png',i)),'Resolution',400,'BackgroundColor','w')
%     set(gcf,'renderer','Painters')
%     exportgraphics(fig,fullfile(directory_out,sprintf('scatter_lightMode_id%d.eps',i)))
%     
%     % save annotated figure dark mode
%     set(gcf,'Color','k')
%     subplot(1,3,1); set(gca,'Color','k','xColor','w','yColor','w'); c1.Color = 'w';
%     subplot(1,3,2); set(gca,'Color','k','xColor','w','yColor','w'); c2.Color = 'w';
%     subplot(1,3,3); set(gca,'Color','k','xColor','w','yColor','w'); c3.Color = 'w';
%     sgtitle(filename,'Interpreter','none','FontWeight','normal','Color','w')
%     savefig(fig,fullfile(directory_out,sprintf('scatter_darkMode_id%d.fig',i)))
%     exportgraphics(fig,fullfile(directory_out,sprintf('scatter_darkMode_id%d.png',i)),'Resolution',400,'BackgroundColor','k')
%     set(gcf,'renderer','Painters')
%     exportgraphics(fig,fullfile(directory_out,sprintf('scatter_darkMode_id%d.eps',i)),'BackgroundColor','k')
%     
%     close all
% end

%% Overview plot lightMode

alpha = 0.2;
min_total = -1e5; % minimum value used for axis limits

for i=1:size(I_upperLobe,2)
% for i=1:1
    fig = figure('Position',[50,50,1660,650]);
    
    t = tiledlayout(2,5);

    x = results_curated.upperLobe(i).intensityTrace;
    y = results_curated.lowerLobe(i).intensityTrace;
    z = results_curated.thirdLobe(i).intensityTrace;
    c = 1:numel(x);

    max_res = max([x;y]);
    max_third = max(z);
    max_total = max([max_res max_third]);

    nexttile([1 3])
    plot(results_curated.upperLobe(i).intensityTrace + 20e5,'Color',[cols(1,:) 1]); hold on
    plot(results_curated.lowerLobe(i).intensityTrace + 10e5,'Color',[cols(2,:) 1]);
    plot(results_curated.thirdLobe(i).intensityTrace,'Color',[cols(3,:) 1]); hold on
    box off
    legend('upper','lower','third')
    xlabel('Time (frames)'); ylabel('Intensity (adu)');

    nexttile(6)
    scatter(y,x,10,c,'filled','MarkerFaceAlpha',alpha)
    xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
    xlabel('Intensity lower lobe'); ylabel('Intensity upper lobe');
    colormap turbo; c1 = colorbar; c1.Label.String = 'Frame';
    
    nexttile(7)
    scatter(z,x,10,c,'filled','MarkerFaceAlpha',alpha)
    xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
    xlabel('Intensity third lobe'); ylabel('Intensity upper lobe');
    colormap turbo; c2 = colorbar; c2.Label.String = 'Frame';
    
    nexttile(8)
    scatter(z,y,10,c,'filled','MarkerFaceAlpha',alpha)
    xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
    xlabel('Intensity third lobe'); ylabel('Intensity lower lobe');
    colormap turbo; c3 = colorbar; c3.Label.String = 'Frame';

    nexttile([2 2])
    mip_rgb = uint8(255*results_curated.mip./max(results_curated.mip(:)));
    mip_rgb = cat(3,mip_rgb,mip_rgb,mip_rgb);
    imagesc(mip_rgb); axis image; axis off
    hold on
    scatter(results_curated.upperLobe(i).y/1e9,results_curated.upperLobe(i).x/1e9,ms,col1); hold on
    scatter(results_curated.lowerLobe(i).y/1e9,results_curated.lowerLobe(i).x/1e9,ms,col2);
    scatter(results_curated.thirdLobe(i).y/1e9,results_curated.thirdLobe(i).x/1e9,ms,col3);
    x_text = (results_curated.upperLobe(i).y/1e9 + results_curated.lowerLobe(i).y/1e9)/2; 
    y_text = (results_curated.upperLobe(i).x/1e9 + results_curated.lowerLobe(i).x/1e9)/2; 
    text(x_text,y_text,sprintf('%d',i),'Color','w')
    pause(0.01)

    % save figure light mode
    sgtitle(filename,'Interpreter','none','FontWeight','normal')
    savefig(fig,fullfile(directory_out,sprintf('overview_lightMode_id%d.fig',i)))
    exportgraphics(fig,fullfile(directory_out,sprintf('overview_lightMode_id%d.png',i)),'Resolution',400,'BackgroundColor','w')
    set(gcf,'renderer','Painters')
    exportgraphics(fig,fullfile(directory_out,sprintf('overview_lightMode_id%d.eps',i)))
    
    close all
end

%% Overview plot dark mode

alpha = 0.2;
min_total = -1e5; % minimum value used for axis limits

for i=1:size(I_upperLobe,2)
% for i=1:1
    fig = figure('Position',[50,50,1660,650]);
    
    t = tiledlayout(2,5);
    set(gcf,'Color','k')

    x = results_curated.upperLobe(i).intensityTrace;
    y = results_curated.lowerLobe(i).intensityTrace;
    z = results_curated.thirdLobe(i).intensityTrace;
    c = 1:numel(x);

    max_res = max([x;y]);
    max_third = max(z);
    max_total = max([max_res max_third]);

    nexttile([1 3])
    plot(results_curated.upperLobe(i).intensityTrace + 20e5,'Color',[cols(1,:) 1]); hold on
    plot(results_curated.lowerLobe(i).intensityTrace + 10e5,'Color',[cols(2,:) 1]);
    plot(results_curated.thirdLobe(i).intensityTrace,'Color',[cols(3,:) 1]); hold on
    box off
    legend('upper','lower','third','TextColor','w','EdgeColor','w','Color','k')
    xlabel('Time (frames)'); ylabel('Intensity (adu)');
    set(gca,'Color','k','xColor','w','yColor','w');

    nexttile(6)
    scatter(y,x,10,c,'filled','MarkerFaceAlpha',alpha)
    xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
    xlabel('Intensity lower lobe'); ylabel('Intensity upper lobe');
    colormap turbo; c1 = colorbar; c1.Label.String = 'Frame'; c1.Color = 'w';
    set(gca,'Color','k','xColor','w','yColor','w');
    
    nexttile(7)
    scatter(z,x,10,c,'filled','MarkerFaceAlpha',alpha)
    xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
    xlabel('Intensity third lobe'); ylabel('Intensity upper lobe');
    colormap turbo; c2 = colorbar; c2.Label.String = 'Frame'; c2.Color = 'w';
    set(gca,'Color','k','xColor','w','yColor','w');
    
    nexttile(8)
    scatter(z,y,10,c,'filled','MarkerFaceAlpha',alpha)
    xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
    xlabel('Intensity third lobe'); ylabel('Intensity lower lobe');
    colormap turbo; c3 = colorbar; c3.Label.String = 'Frame'; c3.Color = 'w';
    set(gca,'Color','k','xColor','w','yColor','w');

    nexttile([2 2])
    mip_rgb = uint8(255*results_curated.mip./max(results_curated.mip(:)));
    mip_rgb = cat(3,mip_rgb,mip_rgb,mip_rgb);
    imagesc(mip_rgb); axis image; axis off
    hold on
    scatter(results_curated.upperLobe(i).y/1e9,results_curated.upperLobe(i).x/1e9,ms,col1); hold on
    scatter(results_curated.lowerLobe(i).y/1e9,results_curated.lowerLobe(i).x/1e9,ms,col2);
    scatter(results_curated.thirdLobe(i).y/1e9,results_curated.thirdLobe(i).x/1e9,ms,col3);
    x_text = (results_curated.upperLobe(i).y/1e9 + results_curated.lowerLobe(i).y/1e9)/2; 
    y_text = (results_curated.upperLobe(i).x/1e9 + results_curated.lowerLobe(i).x/1e9)/2; 
    text(x_text,y_text,sprintf('%d',i),'Color','w')
    pause(0.01)
 
    % save annotated figure dark mode
    sgtitle(filename,'Interpreter','none','FontWeight','normal','Color','w')
    savefig(fig,fullfile(directory_out,sprintf('scatter_darkMode_id%d.fig',i)))
    exportgraphics(fig,fullfile(directory_out,sprintf('scatter_darkMode_id%d.png',i)),'Resolution',400,'BackgroundColor','k')
    set(gcf,'renderer','Painters')
    exportgraphics(fig,fullfile(directory_out,sprintf('scatter_darkMode_id%d.eps',i)),'BackgroundColor','k')
%     
    close all
end

% %% plot ratios
% 
% minIntThreshold = 0.75e5;
% 
% cols = lines(3);
% alpha = 0.5;
% ms = 1;
% 
% for i=1:size(I_upperLobe,2)
%     fig = figure('Position',[200,200,1080,650]);
%     
%     subplot(4,1,1)
%     plot(results_curated.upperLobe(i).intensityTrace + 20e5,'Color',[cols(1,:) 1]); hold on
%     plot(results_curated.lowerLobe(i).intensityTrace + 10e5,'Color',[cols(2,:) 1]);
%     plot(results_curated.thirdLobe(i).intensityTrace,'Color',[cols(3,:) 1]); hold on
%     box off
%     legend('upper','lower','third')
%     xlabel('Time (frames)'); ylabel('Intensity (adu)');
% 
%     totalInt = results_curated.upperLobe(i).intensityTrace + results_curated.lowerLobe(i).intensityTrace + results_curated.thirdLobe(i).intensityTrace;
%     t = 1:numel(totalInt);
%     keep = totalInt > minIntThreshold;
% 
%     subplot(4,1,2)
%     scatter(t(keep),results_curated.upperLobe(i).intensityTrace(keep)./totalInt(keep),ms,cols(1,:),'o','filled','MarkerFaceAlpha',alpha); hold on
%     scatter(t(keep),results_curated.lowerLobe(i).intensityTrace(keep)./totalInt(keep),ms,cols(2,:),'o','filled','MarkerFaceAlpha',alpha);
%     scatter(t(keep),results_curated.thirdLobe(i).intensityTrace(keep)./totalInt(keep),ms,cols(3,:),'o','filled','MarkerFaceAlpha',alpha);
%     ylim([-0.2 1.2]); xlim([0 numel(totalInt)]); grid on;
%     xlabel('Time (frame)'); ylabel('Intensity fraction')
%     legend('% upper','% lower','% third','location','northeast')
% 
%     subplot(4,1,3)
%     totalIntResonator = results_curated.upperLobe(i).intensityTrace + results_curated.lowerLobe(i).intensityTrace;
%     scatter(t(keep),results_curated.upperLobe(i).intensityTrace(keep)./totalIntResonator(keep),ms,cols(1,:),'o','filled','MarkerFaceAlpha',alpha);hold on
%     scatter(t(keep),results_curated.lowerLobe(i).intensityTrace(keep)./totalIntResonator(keep),ms,cols(2,:),'o','filled','MarkerFaceAlpha',alpha);
%     ylim([-0.2 1.2]); xlim([0 numel(totalInt)]); grid on
%     xlabel('Time (frame)'); ylabel('Intensity fraction')
%     legend('% upper resonator','% lower resonator','location','northeast')
% 
%     subplot(4,1,4)
%     scatter(t(keep),totalIntResonator(keep)./totalInt(keep),ms,(cols(1,:)+cols(2,:))/2,'o','filled','MarkerFaceAlpha',alpha); hold on
%     scatter(t(keep),results_curated.thirdLobe(i).intensityTrace(keep)./totalInt(keep),ms,cols(3,:),'o','filled','MarkerFaceAlpha',alpha);
%     ylim([-0.2 1.2]); xlim([0 numel(totalInt)]); grid on
%     xlabel('Time (frame)'); ylabel('Intensity fraction')
%     legend('% resonator','% lower','% third','location','northeast')
% 
%     % save figure light mode
%     sgtitle(filename,'Interpreter','none','FontWeight','normal')
%     savefig(fig,fullfile(directory_out,sprintf('ratioTraces_lightMode_id%d.fig',i)))
%     exportgraphics(fig,fullfile(directory_out,sprintf('ratioTraces_lightMode_id%d.png',i)),'Resolution',400,'BackgroundColor','w')
%     set(gcf,'renderer','Painters')
%     exportgraphics(fig,fullfile(directory_out,sprintf('ratioTraces_lightMode_id%d.eps',i)))
%     
%     % save annotated figure dark mode
%     set(gcf,'Color','k')
%     subplot(4,1,1); set(gca,'Color','k','xColor','w','yColor','w'); legend('Color','k','TextColor','w','EdgeColor','w')
%     subplot(4,1,2); set(gca,'Color','k','xColor','w','yColor','w'); legend('Color','k','TextColor','w','EdgeColor','w')
%     subplot(4,1,3); set(gca,'Color','k','xColor','w','yColor','w'); legend('Color','k','TextColor','w','EdgeColor','w')
%     subplot(4,1,4); set(gca,'Color','k','xColor','w','yColor','w'); legend('Color','k','TextColor','w','EdgeColor','w')
%     sgtitle(filename,'Interpreter','none','FontWeight','normal','Color','w')
%     savefig(fig,fullfile(directory_out,sprintf('ratioTraces_darkMode_id%d.fig',i)))
%     exportgraphics(fig,fullfile(directory_out,sprintf('ratioTraces_darkMode_id%d.png',i)),'Resolution',400,'BackgroundColor','k')
%     set(gcf,'renderer','Painters')
%     exportgraphics(fig,fullfile(directory_out,sprintf('ratioTraces_darkMode_id%d.eps',i)),'BackgroundColor','k')
%     
%     close all
% end



%% Plot the rate of blueing in emission and excitation




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
