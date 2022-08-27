clear all
close all
clc
addpath('lib')

directory_out = fullfile(pwd,'results','simulated spectral blueing illustration');
if ~exist(directory_out,'dir'); mkdir(directory_out); end

fontsize = 10;

col_lightgray = 0.8*[1 1 1];
col_405 = [135 90 175]/255; % 405 nm
col_488 = [2 100 173]/255; % 488 nm
col_561 = [177 212 55]/255; % 561 nm
col_638 = [212 0 0]/255; % 638 nm

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

wavelengthRange   = [300 700]; % nm
wavelengthBleuing = 200; % nm

writeLightModeGif = 1;
writeDarkModeGif = 1;


%% Generate light mode blueing animation

wavelengthStep = 1; % nm
wavelength = wavelengthRange(1):wavelengthStep:(wavelengthRange(2)+wavelengthBleuing);

[~,emission]   = resampleSpectrum(T_qdot.wavelength,T_qdot.em,wavelength);
[~,excitation] = resampleSpectrum(T_qdot.wavelength,T_qdot.ex,wavelength);
[~,emfilter]   = resampleSpectrum(T_emfilter.wavelength,T_emfilter.transmission,wavelength);
[~,dichroic]   = resampleSpectrum(T_dichroic.wavelength,T_dichroic.transmission,wavelength);
[~,qe]         = resampleSpectrum(T_qe.wavelength,T_qe.qe,wavelength);

% normalise emission to integrate to 1
areaEmission = sum(emission(:));
emission = emission/areaEmission;
excitation = excitation/excitation(1);

rgbTriplet = Utils.getRGBtripletFromWavelength(655); 

fig = figure('Position',[50 200 1000 250]);
subplot(2,2,1)
plot(wavelength,excitation,'Color',rgbTriplet,'LineWidth',1.5); hold on
% plot(wavelength,dichroic,'Color',col_lightgray)
line([488 488],[0 1],'LineWidth',2,'Color',col_488);
line([638 638],[0 1],'LineWidth',2,'Color',col_638);
xlim(wavelengthRange); ylim([0 1]); box off
set(gca,'Layer','top')
xlabel('Wavelength (nm)')
ylabel('Excitation spectrum')

subplot(2,2,2)
plot(wavelength,emission,'Color',rgbTriplet,'LineWidth',1.5); hold on
% plot(wavelength,emfilter,'k')
% plot(wavelength,dichroic,'Color',0.8*[1 1 1])
% plot(wavelength,qe,'--k')
xlim(wavelengthRange); ylim([0 1]); box off
set(gca,'Layer','top')
xlabel('Wavelength (nm)')
ylabel('Emission spectrum')

for i=1:wavelengthBleuing
    
    rgbTriplet = Utils.getRGBtripletFromWavelength(peakEmissionWavelength - i); 
    [emission,excitation] = blueSpectra(emission,excitation);

    if ~mod(i,20)
        
        subplot(1,2,1)
        plot(wavelength,excitation,'Color',rgbTriplet,'LineWidth',1.5); hold on
        xlim(wavelengthRange); ylim([0 1]); box off
        set(gca,'Layer','top')
        xlabel('Wavelength (nm)')
        ylabel('Excitation spectrum')
    
        subplot(1,2,2)
        plot(wavelength,areaEmission*emission,'Color',rgbTriplet,'LineWidth',1.5); hold on
        xlim(wavelengthRange); ylim([0 1]); box off
        set(gca,'Layer','top')
        xlabel('Wavelength (nm)')
        ylabel('Emission spectrum')
    
        pause(0.00001)
    end

end

figureBasename = sprintf('experimentSpectra_illustration_%s',qdotName);

% save light mode version
savefig(fig,fullfile(directory_out,strcat(figureBasename,'_lightMode.fig')))
exportgraphics(fig,fullfile(directory_out,strcat(figureBasename,'_lightMode.png')),'Resolution',400)
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,strcat(figureBasename,'_lightMode.eps')))

% convert to dark mode
set(gcf,'Color','k')
subplot(1,2,1);
set(gca,'Color','k','XColor','w','YColor','w','Layer','top');
subplot(1,2,2);
set(gca,'Color','k','XColor','w','YColor','w','Layer','top')

% save dark mode version
savefig(fig,fullfile(directory_out,strcat(figureBasename,'_darkMode.fig')))
exportgraphics(fig,fullfile(directory_out,strcat(figureBasename,'_darkMode.png')),'Resolution',400,'BackgroundColor','k')
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,strcat(figureBasename,'_darkMode.eps')),'BackgroundColor','k')


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