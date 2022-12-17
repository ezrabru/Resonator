clear all
close all
clc
addpath('lib')

directories = {'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi1\results',...
               'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi2\results',...
               'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi3\results',...
               'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi4\results',...
               'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi5\results',...
               'D:\manuscripts\resonator\new data\20201013_resonator 2 lobe sparser\data\roi6\results'}; % path to folder with tif stacks

for id_dir = 3:numel(directories)


    %% Load data
    
    directory = directories{id_dir};
    load(fullfile(directory,'results_curated.mat'))
    load(fullfile(directory,'results.mat'))
    
    
    %% Colours
    
    col1 = 1.05*[0.8392,0.1529,0.1569]; % red
    col2 = [0.1333,0.5289,0.8000]; % blue
    
    % Create new output folder
    directory_out = fullfile(directory,'figures');
    if ~exist(directory_out,'dir'); mkdir(directory_out); end
    
    
    %% Overview dark mode
    
    alpha = 1;
    min_total = -1e5; % minimum value used for axis limits
    minIntThreshold = 1e5;
    ms = 1;
    
    colourcode_axislabels = 0;
    [~,filename,~] = fileparts(results(1).filepath);
    
    for i=1:height(results_curated.upperLobe)
        fig = figure('Position',[50,50,1700,700]);
    
        tiledlayout(4,9);
        set(gcf,'Color','k')
    
        x = results_curated.upperLobe(i).intensityTrace;
        y = results_curated.lowerLobe(i).intensityTrace;
        t = 1:numel(x);
    
        max_total = max([x;y]);
    
        nexttile([4 4])
        mip_rgb = uint8(255*results_curated.mip./max(results_curated.mip(:)));
        mip_rgb = cat(3,mip_rgb,mip_rgb,mip_rgb);
        imagesc(mip_rgb); axis image; axis off
        hold on
        scatter(results_curated.upperLobe(i).y/1e9,results_curated.upperLobe(i).x/1e9,70,col1); hold on
        scatter(results_curated.lowerLobe(i).y/1e9,results_curated.lowerLobe(i).x/1e9,70,col2);
        x_text = (results_curated.upperLobe(i).y/1e9 + results_curated.lowerLobe(i).y/1e9)/2;
        y_text = (results_curated.upperLobe(i).x/1e9 + results_curated.lowerLobe(i).x/1e9)/2;
        text(x_text,y_text,sprintf('%d',i),'Color','w')
        title(sprintf('%s, QDot %d',filename,i),'Fontweight','normal','Interpreter','none','Color','w')
    
    
        nexttile([1 5])
        plot(results_curated.upperLobe(i).intensityTrace + 10e5,'Color',col1); hold on
        plot(results_curated.lowerLobe(i).intensityTrace,'Color',col2);
        box off
        legend('upper','lower','Location','bestoutside','Color','k','TextColor','w','EdgeColor','none')
        xlabel('Time (frames)'); ylabel('Intensity (adu)');
        set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    
        totalInt = results_curated.upperLobe(i).intensityTrace + results_curated.lowerLobe(i).intensityTrace;
        keep = totalInt > minIntThreshold;
    
        nexttile([1 5])
        scatter(t(keep),results_curated.upperLobe(i).intensityTrace(keep)./totalInt(keep),ms,col1,'o','filled','MarkerFaceAlpha',alpha); hold on
        scatter(t(keep),results_curated.lowerLobe(i).intensityTrace(keep)./totalInt(keep),ms,col2,'o','filled','MarkerFaceAlpha',alpha);
        ylim([-0.2 1.2]); xlim([0 numel(totalInt)]); grid on;
        xlabel('Time (frame)'); ylabel('Intensity fraction')
        legend('% upper','% lower','Location','bestoutside','Color','k','TextColor','w','EdgeColor','none')
        set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    
        max_total = max([x;y]);
    
        nexttile([2 2])
        scatter(y,x,5,t,'filled','MarkerFaceAlpha',alpha)
        xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
        xlabel('Intensity lower lobe'); ylabel('Intensity upper lobe');
        colormap turbo; c = colorbar; c.Label.String = 'Frame'; c.Color = 'w';
        set(gca,'Color','k','XColor','w','YColor','w','Layer','top')

%         x = results_curated.upperLobe(i).intensityTrace(keep);
%         y = results_curated.lowerLobe(i).intensityTrace(keep);
%         normInt = sqrt(x.^2 + y.^2);
%         % slope = atan2(x./normInt,y./normInt);
%         slope = (x-y)./(x + y);

%         nexttile([2 3])
%         scatter(t(keep),slope,ms,'w','o','filled','MarkerFaceAlpha',alpha); grid on; hold on
%         scatter(t(keep),smooth(slope,1000),3*ms,'r','o','filled','MarkerFaceAlpha',alpha);
%         plot(t(keep),smooth(slope,500),'r','LineWidth',1.5); grid on
%         xlabel('Time (frame)'); ylabel('Intensity upper lobe');
%         ylim([-1.5 1.5])
%         colormap turbo; c = colorbar; c.Label.String = 'Frame'; c.Color = 'w';
%         set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    
        pause(0.001)
        
        % save annotated figure dark mode
        savefig(fig,fullfile(directory_out,sprintf('overview_darkMode_id%d.fig',i)))
        exportgraphics(fig,fullfile(directory_out,sprintf('overview_darkMode_id%d.png',i)),'Resolution',400,'BackgroundColor','k')
        set(gcf,'renderer','Painters')
        exportgraphics(fig,fullfile(directory_out,sprintf('overview_darkMode_id%d.eps',i)),'BackgroundColor','k')
        
        close all
    end
    
    
%     %% Overview light mode
%     
%     alpha = 1;
%     min_total = -1e5; % minimum value used for axis limits
%     minIntThreshold = 1e5;
%     ms = 1;
%     
%     colourcode_axislabels = 0;
%     [~,filename,~] = fileparts(results(1).filepath);
%     
%     for i=1:height(results_curated.upperLobe)
%         fig = figure('Position',[50,50,1700,700]);
%     
%         t = tiledlayout(4,9);
%         set(gcf,'Color','w')
%     
%         x = results_curated.upperLobe(i).intensityTrace;
%         y = results_curated.lowerLobe(i).intensityTrace;
%         c = 1:numel(x);
%     
%         max_total = max([x;y]);
%     
%         nexttile([4 4])
%         mip_rgb = uint8(255*results_curated.mip./max(results_curated.mip(:)));
%         mip_rgb = cat(3,mip_rgb,mip_rgb,mip_rgb);
%         imagesc(mip_rgb); axis image; axis off
%         hold on
%         scatter(results_curated.upperLobe(i).y/1e9,results_curated.upperLobe(i).x/1e9,70,col1); hold on
%         scatter(results_curated.lowerLobe(i).y/1e9,results_curated.lowerLobe(i).x/1e9,70,col2);
%         x_text = (results_curated.upperLobe(i).y/1e9 + results_curated.lowerLobe(i).y/1e9)/2;
%         y_text = (results_curated.upperLobe(i).x/1e9 + results_curated.lowerLobe(i).x/1e9)/2;
%         text(x_text,y_text,sprintf('%d',i),'Color','w')
%         title(sprintf('%s, QDot %d',filename,i),'Fontweight','normal','Interpreter','none')
%     
%     
%         nexttile([1 5])
%         plot(results_curated.upperLobe(i).intensityTrace + 10e5,'Color',col1); hold on
%         plot(results_curated.lowerLobe(i).intensityTrace,'Color',col2);
%         box off
%         legend('upper','lower','Location','bestoutside','EdgeColor','none')
%         xlabel('Time (frames)'); ylabel('Intensity (adu)');
%         set(gca,'Layer','top')
%     
%         totalInt = results_curated.upperLobe(i).intensityTrace + results_curated.lowerLobe(i).intensityTrace;
%         t = 1:numel(totalInt);
%         keep = totalInt > minIntThreshold;
%     
%         nexttile([1 5])
%         scatter(t(keep),results_curated.upperLobe(i).intensityTrace(keep)./totalInt(keep),ms,col1,'o','filled','MarkerFaceAlpha',alpha); hold on
%         scatter(t(keep),results_curated.lowerLobe(i).intensityTrace(keep)./totalInt(keep),ms,col2,'o','filled','MarkerFaceAlpha',alpha);
%         ylim([-0.2 1.2]); xlim([0 numel(totalInt)]); grid on;
%         xlabel('Time (frame)'); ylabel('Intensity fraction')
%         legend('% upper','% lower','Location','bestoutside','EdgeColor','none')
%         set(gca,'Layer','top')
%     
%         max_total = max([x;y]);
%     
%         nexttile([2 2])
%         scatter(y,x,5,c,'filled','MarkerFaceAlpha',alpha)
%         xlim([min_total max_total]); ylim([min_total max_total]); axis square; grid on
%         xlabel('Intensity lower lobe'); ylabel('Intensity upper lobe');
%         colormap turbo; c = colorbar; c.Label.String = 'Frame';
%         set(gca,'Layer','top')
%     
%         pause(0.001)
%         
%         % save annotated figure dark mode
%         savefig(fig,fullfile(directory_out,sprintf('overview_lightMode_id%d.fig',i)))
%         exportgraphics(fig,fullfile(directory_out,sprintf('overview_lightMode_id%d.png',i)),'Resolution',400)
%         set(gcf,'renderer','Painters')
%         exportgraphics(fig,fullfile(directory_out,sprintf('overview_lightMode_id%d.eps',i)))
%         
%         close all
%     end
end

