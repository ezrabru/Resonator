clear all
close all
clc
addpath('lib')

directory_out = fullfile(pwd,'simulated spectral blueing 3 lobe');
if ~exist(directory_out,'dir'); mkdir(directory_out); end

fontsize = 10;

col_lightgray = 0.8*[1 1 1];
col_405 = [135 90 175]/255; % 405 nm
col_488 = [2 100 173]/255; % 488 nm
col_561 = [177 212 55]/255; % 561 nm
col_638 = [212 0 0]/255; % 638 nm

% qdot585 = readmatrix(fullfile(pwd,'spectra dyes and filters','dyes','Qdot585.txt'));
% T_qdot585.wavelength = qdot585(:,1);
% T_qdot585.ex = qdot585(:,2)/100;
% T_qdot585.em = qdot585(:,3)/100;

qdot655 = readmatrix(fullfile(pwd,'spectra dyes and filters','dyes','Qdot655.txt'));
T_qdot655.wavelength = qdot655(:,1);
T_qdot655.ex = qdot655(:,2)/100;
T_qdot655.em = qdot655(:,3)/100;

% exfilter = readmatrix(fullfile(pwd,'spectra dyes and filters','filters','ZET405-488-561-640xv2.txt'));
% T_exfilter.wavelength = exfilter(:,1);
% T_exfilter.transmission = exfilter(:,2);

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

writeLightModeGif = 1;
writeDarkModeGif = 1;


%% Generate light mode blueing animation

wavelengthStep = 1; % nm
wavelength = wavelengthRange(1):wavelengthStep:(wavelengthRange(2)+wavelengthBleuing);

[~,emission]      = resampleSpectrum(T_qdot655.wavelength,T_qdot655.em,wavelength);
[~,excitation]    = resampleSpectrum(T_qdot655.wavelength,T_qdot655.ex,wavelength);
[~,emfilter]      = resampleSpectrum(T_emfilter.wavelength,T_emfilter.transmission,wavelength);
[~,dichroic]      = resampleSpectrum(T_dichroic.wavelength,T_dichroic.transmission,wavelength);
[~,dichroiclobe3] = resampleSpectrum(T_dichroiclobe3.wavelength,T_dichroiclobe3.transmission,wavelength);
[~,qe]            = resampleSpectrum(T_qe.wavelength,T_qe.qe,wavelength);

% normalise emission to integrate to 1
areaEmission = sum(emission(:));
emission = emission/areaEmission;

% get the powers of the two lasers such that the intensity of both
% resonator PSF lobes are equal. This is true when:
% P488*epsilon488*lambda488 = P638*epsilon638*lambda638
[~,idx_488] = min(abs(wavelength - 488));
[~,idx_638] = min(abs(wavelength - 638));
epsilon488 = excitation(idx_488);
epsilon638 = excitation(idx_638);
power488 = 1;
power638 = power488*(epsilon488/epsilon638)*(488/638);

T_sys = emfilter.*dichroic.*qe;

% initialise
totalFluorescence = nan(numel(wavelengthBleuing),1);
absorption488 = nan(numel(wavelengthBleuing),1);
absorption638 = nan(numel(wavelengthBleuing),1);
lobe488 = nan(numel(wavelengthBleuing),1);
lobe638 = nan(numel(wavelengthBleuing),1);
lobe3 = nan(numel(wavelengthBleuing),1);

figure('Position',[50 200 1000 500])
for i=1:wavelengthBleuing

    rgbTriplet = Utils.getRGBtripletFromWavelength(655 - i); 

    [emission,excitation] = blueSpectra(emission,excitation);
    
    [~,idx_488] = min(abs(wavelength - 488));
    [~,idx_638] = min(abs(wavelength - 638));
    absorption488(i) = excitation(idx_488);
    absorption638(i) = excitation(idx_638);
    F488 = power488*absorption488(i)*488;
    F638 = power638*absorption638(i)*638;

    lobe3(i)   = sum(T_sys.*(1 - dichroiclobe3).*emission*(F488 + F638),'all');
    lobe488(i) = sum(T_sys.*dichroiclobe3.*emission*F488,'all');
    lobe638(i) = sum(T_sys.*dichroiclobe3.*emission*F638,'all');

    subplot(2,2,1)
    plot(wavelength,excitation,'Color',rgbTriplet,'LineWidth',1.5); hold on
    plot(wavelength,dichroic,'Color',col_lightgray)
    line([488 488],[0 1],'LineWidth',2,'Color',col_488);
    line([638 638],[0 1],'LineWidth',2,'Color',col_638);
    xlim(wavelengthRange); ylim([0 1]); box off
    set(gca,'Layer','top')
    xlabel('Wavelength (nm)')
    ylabel('Spectrum')

    subplot(2,2,2)
    plot(wavelength,areaEmission*emission,'Color',rgbTriplet,'LineWidth',1.5); hold on
    plot(wavelength,emfilter,'k')
    plot(wavelength,dichroic,'Color',0.8*[1 1 1])
    plot(wavelength,qe,'--k')
    plot(wavelength,dichroiclobe3,'b')
    fill(wavelength,areaEmission*T_sys.*emission,rgbTriplet,'FaceAlpha',0.4,'EdgeColor','none')
    xlim(wavelengthRange); ylim([0 1]); box off
    set(gca,'Layer','top')
    xlabel('Wavelength (nm)')
    ylabel('Spectrum')

    subplot(2,2,3)
    plot(1:i,lobe3,'Color',col_lightgray,'LineWidth',1.5); hold on
    plot(1:i,lobe488,'Color',col_488,'LineWidth',1.5);
    plot(1:i,lobe638,'Color',col_638,'LineWidth',1.5);
    xlim([0 wavelengthBleuing]); box off
    xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)')
    ylabel('Lobe intensity')

    subplot(2,2,4)
    plot(1:i,lobe3./(lobe488+lobe638+lobe3),'Color',col_lightgray,'LineWidth',1.5); hold on
    plot(1:i,lobe488./(lobe488+lobe638+lobe3),'Color',col_488,'LineWidth',1.5);
    plot(1:i,lobe638./(lobe488+lobe638+lobe3),'Color',col_638,'LineWidth',1.5);
    xlim([0 wavelengthBleuing]); box off
    ylim([0 1])
    xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)'); box off
    ylabel('Lobe intensity/Total intensity')

    pause(0.00001)

    if writeLightModeGif
        if i == 1
            exportgraphics(gcf,fullfile(directory_out,'blueingGif_lightmode.gif'),'Resolution',400);
        else
            exportgraphics(gcf,fullfile(directory_out,'blueingGif_lightmode.gif'),'Append',true,'Resolution',400);
        end
    end

    if i ~= wavelengthBleuing; clf; end
end

%% Generate light mode blueing animation

wavelengthStep = 1; % nm
wavelength = wavelengthRange(1):wavelengthStep:(wavelengthRange(2)+wavelengthBleuing);

[~,emission]      = resampleSpectrum(T_qdot655.wavelength,T_qdot655.em,wavelength);
[~,excitation]    = resampleSpectrum(T_qdot655.wavelength,T_qdot655.ex,wavelength);
[~,emfilter]      = resampleSpectrum(T_emfilter.wavelength,T_emfilter.transmission,wavelength);
[~,dichroic]      = resampleSpectrum(T_dichroic.wavelength,T_dichroic.transmission,wavelength);
[~,dichroiclobe3] = resampleSpectrum(T_dichroiclobe3.wavelength,T_dichroiclobe3.transmission,wavelength);
[~,qe]            = resampleSpectrum(T_qe.wavelength,T_qe.qe,wavelength);

% normalise emission to integrate to 1
areaEmission = sum(emission(:));
emission = emission/areaEmission;

% get the powers of the two lasers such that the intensity of both
% resonator PSF lobes are equal. This is true when:
% P488*epsilon488*lambda488 = P638*epsilon638*lambda638
[~,idx_488] = min(abs(wavelength - 488));
[~,idx_638] = min(abs(wavelength - 638));
epsilon488 = excitation(idx_488);
epsilon638 = excitation(idx_638);
power488 = 1;
power638 = power488*(epsilon488/epsilon638)*(488/638);

T_sys = emfilter.*dichroic.*qe;

% initialise
totalFluorescence = nan(numel(wavelengthBleuing),1);
absorption488 = nan(numel(wavelengthBleuing),1);
absorption638 = nan(numel(wavelengthBleuing),1);
lobe488 = nan(numel(wavelengthBleuing),1);
lobe638 = nan(numel(wavelengthBleuing),1);
lobe3 = nan(numel(wavelengthBleuing),1);

figure('Position',[50 200 1000 500])
for i=1:wavelengthBleuing

    rgbTriplet = Utils.getRGBtripletFromWavelength(655 - i); 

    [emission,excitation] = blueSpectra(emission,excitation);
    
    [~,idx_488] = min(abs(wavelength - 488));
    [~,idx_638] = min(abs(wavelength - 638));
    absorption488(i) = excitation(idx_488);
    absorption638(i) = excitation(idx_638);
    F488 = power488*absorption488(i)*488;
    F638 = power638*absorption638(i)*638;

    lobe3(i)   = sum(T_sys.*(1 - dichroiclobe3).*emission*(F488 + F638),'all');
    lobe488(i) = sum(T_sys.*dichroiclobe3.*emission*F488,'all');
    lobe638(i) = sum(T_sys.*dichroiclobe3.*emission*F638,'all');

    subplot(2,2,1)
    plot(wavelength,excitation,'Color',rgbTriplet,'LineWidth',1.5); hold on
    plot(wavelength,dichroic,'Color',0.4*col_lightgray)
    line([488 488],[0 1],'LineWidth',2,'Color',col_488);
    line([638 638],[0 1],'LineWidth',2,'Color',col_638);
    xlim(wavelengthRange); ylim([0 1]); box off
    set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    xlabel('Wavelength (nm)')
    ylabel('Spectrum')

    subplot(2,2,2)
    plot(wavelength,areaEmission*emission,'Color',rgbTriplet,'LineWidth',1.5); hold on
    plot(wavelength,emfilter,'w')
    plot(wavelength,dichroic,'Color',0.4*col_lightgray)
    plot(wavelength,qe,'--w')
    plot(wavelength,dichroiclobe3,'g')
    fill(wavelength,areaEmission*T_sys.*emission,rgbTriplet,'FaceAlpha',0.4,'EdgeColor','none')
    xlim(wavelengthRange); ylim([0 1]); box off
    set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
    xlabel('Wavelength (nm)')
    ylabel('Spectrum')

    subplot(2,2,3)
    plot(1:i,lobe3,'Color',col_lightgray,'LineWidth',1.5); hold on
    plot(1:i,lobe488,'Color',col_488,'LineWidth',1.5);
    plot(1:i,lobe638,'Color',col_638,'LineWidth',1.5);
    xlim([0 wavelengthBleuing]); box off
    xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)')
    ylabel('Lobe intensity')
    set(gca,'Color','k','XColor','w','YColor','w','Layer','top')

    subplot(2,2,4)
    plot(1:i,lobe3./(lobe488+lobe638+lobe3),'Color',col_lightgray,'LineWidth',1.5); hold on
    plot(1:i,lobe488./(lobe488+lobe638+lobe3),'Color',col_488,'LineWidth',1.5);
    plot(1:i,lobe638./(lobe488+lobe638+lobe3),'Color',col_638,'LineWidth',1.5);
    xlim([0 wavelengthBleuing]); box off
    ylim([0 1])
    xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)'); box off
    ylabel('Lobe intensity/Total intensity')
    set(gca,'Color','k','XColor','w','YColor','w','Layer','top')

    set(gcf,'Color','k')
    pause(0.00001)

    if writeDarkModeGif
        if i == 1
            exportgraphics(gcf,fullfile(directory_out,'blueingGif_darkmode.gif'),'BackgroundColor','k','Resolution',400);
        else
            exportgraphics(gcf,fullfile(directory_out,'blueingGif_darkmode.gif'),'Append',true,'BackgroundColor','k','Resolution',400);
        end
    end
    
    if i ~= wavelengthBleuing; clf; end
end


%%

lw = 1.5;
fontsize = 10;
legendNames = {'\lambda_{em} > 630 nm, \lambda_{ex} = 488 nm',...
               '\lambda_{em} > 630 nm, \lambda_{ex} = 638 nm',...
               '\lambda_{em} < 630 nm, \lambda_{ex} = 488 and 638 nm'};


fig = figure('Position',[50 200 1000 320]);

subplot(1,2,1)
plot(1:wavelengthBleuing,lobe488,'Color',col_488,'LineWidth',lw); hold on
plot(1:wavelengthBleuing,lobe638,'Color',col_638,'LineWidth',lw);
plot(1:wavelengthBleuing,lobe3,'Color',col_lightgray,'LineWidth',lw);
legend(legendNames,'Location','northoutside','EdgeColor','none')
ylim([0 inf]); box off
xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)')
ylabel('Lobe intensity')
set(gca,'FontSize',fontsize)

subplot(1,2,2)
plot(1:wavelengthBleuing,lobe488./(lobe488+lobe638+lobe3),'Color',col_488,'LineWidth',lw); hold on
plot(1:wavelengthBleuing,lobe638./(lobe488+lobe638+lobe3),'Color',col_638,'LineWidth',lw);
plot(1:wavelengthBleuing,lobe3./(lobe488+lobe638+lobe3),'Color',col_lightgray,'LineWidth',lw);
legend(legendNames,'Location','northoutside','EdgeColor','none')
xlabel('Spectral shift, \Delta\lambda_{blueing} (nm)'); box off
ylabel('Lobe intensity/Total intensity')
set(gca,'FontSize',fontsize)

set(gcf,'Color','w')


% save light mode version
savefig(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_lightMode.fig'))
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_lightMode.png'),'Resolution',400)
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_lightMode.eps'))

% convert to dark mode
set(gcf,'Color','k')
subplot(1,2,1);
set(gca,'Color','k','XColor','w','YColor','w','Layer','top');
legend(legendNames,'Location','northoutside','Color','k','TextColor','w','EdgeColor','none')
subplot(1,2,2); set(gca,'Color','k','XColor','w','YColor','w','Layer','top')
legend(legendNames,'Location','northoutside','Color','k','TextColor','w','EdgeColor','none')

% save dark mode version
savefig(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_darkMode.fig'))
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_darkMode.png'),'Resolution',400,'BackgroundColor','k')
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_darkMode.eps'),'BackgroundColor','k')


%% Function

function [emission,excitation] = blueSpectra(emission,excitation)
emission = blueArrayStep(emission);
excitation = blueArrayStep(excitation);
end

function X_new = blueArrayStep(X)
X_new = X(:);
X_new = X_new(2:end);
X_new = [X_new; X_new(end)];
end

function [xq,vq] = resampleSpectrum(x,v,xq)
vq = interp1(x(:),v(:),xq(:),'linear');
vq(xq < x(1)) = v(1);
vq(xq > x(end)) = v(end);
vq = vq';
vq = vq(:);
end