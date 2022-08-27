clear all
close all
clc
addpath('lib')

directory_out = fullfile(pwd,'results','simulated spectral blueing 3 lobe');
if ~exist(directory_out,'dir'); mkdir(directory_out); end

fontsize = 10;

col_lightgray = 0.8*[1 1 1];
col_405 = [135 90 175]/255; % 405 nm
col_488 = [2 100 173]/255; % 488 nm
col_561 = [177 212 55]/255; % 561 nm
col_638 = [212 0 0]/255; % 638 nm

col1 = col_488;
col2 = col_638;

wavelength1 = 488; % nm
wavelength2 = 638; % nm
dichroicCutoff = 630; % resonator dichroic edge (nm)

% qdotName = 'QDot585';
% peakEmissionWavelength = 585;
% qdot = readmatrix(fullfile(pwd,'spectra dyes and filters','dyes','Qdot585.txt'));
% T_qdot.wavelength = qdot(:,1);
% T_qdot.ex = qdot(:,2)/100;
% T_qdot.em = qdot(:,3)/100;

qdotName = 'QDot655';
peakEmissionWavelength = 655;
qdot = readmatrix(fullfile(pwd,'spectra dyes and filters','dyes','Qdot655.txt'));
T_qdot.wavelength = qdot(:,1);
T_qdot.ex = qdot(:,2)/100;
T_qdot.em = qdot(:,3)/100;

emfilter = readmatrix(fullfile(pwd,'spectra dyes and filters','filters','ZET405-488-561-640mv2.txt'));
T_emfilter.wavelength = emfilter(:,1);
T_emfilter.transmission = emfilter(:,2);

dichroic = readmatrix(fullfile(pwd,'spectra dyes and filters','filters','ZT405-488-561-640rpcv2.txt'));
T_dichroic.wavelength = dichroic(:,1);
T_dichroic.transmission = dichroic(:,2);

qe = readmatrix(fullfile(pwd,'spectra dyes and filters','filters','imagEM_C9100-13_QE.txt'));
T_qe.wavelength = qe(:,1);
T_qe.qe = qe(:,2)/100;

dichroiclobe3 = readmatrix(fullfile(pwd,'spectra dyes and filters','filters','DMLP638.txt'));
T_dichroiclobe3.wavelength = dichroiclobe3(:,1);
T_dichroiclobe3.transmission = dichroiclobe3(:,2)/100;

wavelengthRange   = [350 800]; % nm
wavelengthBleuing = 200; % nm

theme = 'light'; % 'light' or 'dark'
writeGif = 1; % write gif away (1), or only display (0)


%% Generate blueing animation

wavelengthStep = 1; % nm (do not change value)
wavelength = wavelengthRange(1):wavelengthStep:(wavelengthRange(2)+wavelengthBleuing);

[~,emission]      = Utils.resampleSpectrum(T_qdot.wavelength,T_qdot.em,wavelength);
[~,excitation]    = Utils.resampleSpectrum(T_qdot.wavelength,T_qdot.ex,wavelength);
[~,emfilter]      = Utils.resampleSpectrum(T_emfilter.wavelength,T_emfilter.transmission,wavelength);
[~,dichroic]      = Utils.resampleSpectrum(T_dichroic.wavelength,T_dichroic.transmission,wavelength);
[~,dichroiclobe3] = Utils.resampleSpectrum(T_dichroiclobe3.wavelength,T_dichroiclobe3.transmission,wavelength);
[~,qe]            = Utils.resampleSpectrum(T_qe.wavelength,T_qe.qe,wavelength);

% normalise emission to integrate to 1
areaEmission = sum(emission(:));
emission = emission/areaEmission;

% get the powers of the two lasers such that the intensity of both
% resonator PSF lobes are equal. This is true when:
% power1*epsilon1*lambda1 = power2*epsilon2*lambda2
[~,idx1] = min(abs(wavelength - wavelength1));
[~,idx2] = min(abs(wavelength - wavelength2));
epsilon1 = excitation(idx1);
epsilon2 = excitation(idx2);
power1 = 1;
power2 = power1*(epsilon1/epsilon2)*(wavelength1/wavelength2);

T_sys = emfilter.*dichroic.*qe;

% initialise
lobe1 = nan(numel(wavelengthBleuing),1);
lobe2 = nan(numel(wavelengthBleuing),1);
lobe3 = nan(numel(wavelengthBleuing),1);

% get index of two wavelengths (for accessing excitation spectrum values)
[~,idx1] = min(abs(wavelength - wavelength1));
[~,idx2] = min(abs(wavelength - wavelength2));

% prepare figure
figure('Position',[50 200 1000 500]);
if strcmp(theme,'light')
    col_emfilter = 'k';
    col_dcscope = col_lightgray;
    col_dcres = 'b';
    col_qe = '--k';
elseif strcmp(theme,'dark')
    col_emfilter = 'w';
    col_dcscope = 0.5*col_lightgray;
    col_dcres = 'g';
    col_qe = '--w';
    set(gcf,'Color','k');
end

for i=1:wavelengthBleuing
    
    % get colour triplet matching wavelentgh current emission peak
    rgbTriplet = Utils.getRGBtripletFromWavelength(peakEmissionWavelength - i); 
    
    % shift emission and excitation spectra 1 nm to the blue
    [emission,excitation] = Utils.blueSpectra(emission,excitation);
    
    % calculate photon flux at 2 wavelengths
    F1 = power1*excitation(idx1)*wavelength1;
    F2 = power2*excitation(idx2)*wavelength2;

    % calculate measured intensity in 3 PSF lobes
    lobe1(i) = sum(T_sys.*dichroiclobe3.*emission*F1,'all');
    lobe2(i) = sum(T_sys.*dichroiclobe3.*emission*F2,'all');
    lobe3(i) = sum(T_sys.*(1 - dichroiclobe3).*emission*(F1 + F2),'all');
    
    % update plots
    subplot(2,2,1)
    plot(wavelength,excitation,'Color',rgbTriplet,'LineWidth',1.5); hold on
    plot(wavelength,dichroic,'Color',col_dcscope)
    line([wavelength1 wavelength1],[0 1],'LineWidth',2,'Color',col1);
    line([wavelength2 wavelength2],[0 1],'LineWidth',2,'Color',col2);
    xlim(wavelengthRange); ylim([0 1]); box off
    xlabel('Wavelength (nm)'); ylabel('Spectrum')
    if strcmp(theme,'light'); set(gca,'Layer','top')
    elseif strcmp(theme,'dark')
        set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    else; disp('Unsupported value for parameter "theme".'); return
    end

    subplot(2,2,2)
    plot(wavelength,areaEmission*emission,'Color',rgbTriplet,'LineWidth',1.5); hold on
    plot(wavelength,emfilter,col_emfilter)
    plot(wavelength,dichroic,'Color',col_dcscope)
    plot(wavelength,qe,col_qe)
    plot(wavelength,dichroiclobe3,col_dcres)
    fill(wavelength,areaEmission*T_sys.*emission,rgbTriplet,'FaceAlpha',0.4,'EdgeColor','none')
    xlim(wavelengthRange); ylim([0 1]); box off
    xlabel('Wavelength (nm)'), ylabel('Spectrum')
    if strcmp(theme,'light'); set(gca,'Layer','top')
    elseif strcmp(theme,'dark')
        set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    else; disp('Unsupported value for parameter "theme".'); return
    end

    subplot(2,2,3)
    plot(1:i,lobe3,'Color',col_lightgray,'LineWidth',1.5); hold on
    plot(1:i,lobe1,'Color',col1,'LineWidth',1.5);
    plot(1:i,lobe2,'Color',col2,'LineWidth',1.5);
    xlim([0 wavelengthBleuing]); box off
    xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)'); ylabel('Lobe intensity')
    if strcmp(theme,'light'); set(gca,'Layer','top')
    elseif strcmp(theme,'dark')
        set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    else; disp('Unsupported value for parameter "theme".'); return
    end

    subplot(2,2,4)
    plot(1:i,lobe3./(lobe1+lobe2+lobe3),'Color',col_lightgray,'LineWidth',1.5); hold on
    plot(1:i,lobe1./(lobe1+lobe2+lobe3),'Color',col1,'LineWidth',1.5);
    plot(1:i,lobe2./(lobe1+lobe2+lobe3),'Color',col2,'LineWidth',1.5);
    xlim([0 wavelengthBleuing]); box off
    ylim([0 1])
    xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)'); box off
    ylabel('Lobe intensity/Total intensity')
    if strcmp(theme,'light'); set(gca,'Layer','top')
    elseif strcmp(theme,'dark')
        set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    else; disp('Unsupported value for parameter "theme".'); return
    end
    pause(0.00001)

    if writeGif
        if strcmp(theme,'light')
            if i == 1; exportgraphics(gcf,fullfile(directory_out,sprintf('blueingGif_%s_lightmode.gif',qdotName)),'Resolution',400);
            else; exportgraphics(gcf,fullfile(directory_out,sprintf('blueingGif_%s_lightmode.gif',qdotName)),'Append',true,'Resolution',400);
            end
        elseif strcmp(theme,'dark')
            if i == 1; exportgraphics(gcf,fullfile(directory_out,sprintf('blueingGif_%s_darkmode.gif',qdotName)),'BackgroundColor','k','Resolution',400);
            else; exportgraphics(gcf,fullfile(directory_out,sprintf('blueingGif_%s_darkmode.gif',qdotName)),'Append',true,'BackgroundColor','k','Resolution',400);
            end
        end
    end

    if i ~= wavelengthBleuing; clf; end
end


%%

lw = 1.5;
fontsize = 10;
legendNames{1} = strcat('\lambda_{em} >',sprintf(' %d ',dichroicCutoff),' nm, \lambda_{ex} = ',sprintf(' %d ',wavelength1),' nm');
legendNames{2} = strcat('\lambda_{em} >',sprintf(' %d ',dichroicCutoff),' nm, \lambda_{ex} = ',sprintf(' %d ',wavelength2),' nm');
legendNames{3} = strcat('\lambda_{em} <',sprintf(' %d ',dichroicCutoff),' nm, \lambda_{ex} = ',sprintf(' %d ',wavelength1),' and ',sprintf(' %d ',wavelength2),' nm');


fig = figure('Position',[50 200 1000 320]);

subplot(1,2,1)
plot(1:wavelengthBleuing,lobe1,'Color',col1,'LineWidth',lw); hold on
plot(1:wavelengthBleuing,lobe2,'Color',col2,'LineWidth',lw);
plot(1:wavelengthBleuing,lobe3,'Color',col_lightgray,'LineWidth',lw);
legend(legendNames,'Location','northoutside','EdgeColor','none')
ylim([0 inf]); box off
xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)')
ylabel('Lobe intensity')
set(gca,'FontSize',fontsize)

subplot(1,2,2)
plot(1:wavelengthBleuing,lobe1./(lobe1+lobe2+lobe3),'Color',col1,'LineWidth',lw); hold on
plot(1:wavelengthBleuing,lobe2./(lobe1+lobe2+lobe3),'Color',col2,'LineWidth',lw);
plot(1:wavelengthBleuing,lobe3./(lobe1+lobe2+lobe3),'Color',col_lightgray,'LineWidth',lw);
legend(legendNames,'Location','northoutside','EdgeColor','none')
xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)'); box off
ylabel('Lobe intensity/Total intensity')
set(gca,'FontSize',fontsize)
set(gcf,'Color','w')

figureBasename = sprintf('experimentSpectra_3lobe_%s',qdotName);

% save light mode version
savefig(fig,fullfile(directory_out,sprintf('%s_lightMode.fig',figureBasename)))
exportgraphics(fig,fullfile(directory_out,sprintf('%s_lightMode.png',figureBasename)),'Resolution',400)
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,sprintf('%s_lightMode.eps',figureBasename)))

% convert to dark mode
set(gcf,'Color','k')
subplot(1,2,1);
set(gca,'Color','k','XColor','w','YColor','w','Layer','top');
legend(legendNames,'Location','northoutside','Color','k','TextColor','w','EdgeColor','none')
subplot(1,2,2); set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
legend(legendNames,'Location','northoutside','Color','k','TextColor','w','EdgeColor','none')

% save dark mode version
savefig(fig,fullfile(directory_out,sprintf('%s_darkMode.fig',figureBasename)))
exportgraphics(fig,fullfile(directory_out,sprintf('%s_darkMode.png',figureBasename)),'Resolution',400,'BackgroundColor','k')
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,sprintf('%s_darkMode.eps',figureBasename)),'BackgroundColor','k')


