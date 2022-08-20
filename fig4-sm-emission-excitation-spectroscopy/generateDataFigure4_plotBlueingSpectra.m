clear all
close all
clc

directory_out = fullfile(pwd,'simulated spectral blueing');
if ~exist(directory_out,'dir'); mkdir(directory_out); end

fontsize = 8;

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


%% Prepare for multiplication


x = T_emfilter.wavelength;
v = T_emfilter.transmission;
xq = 200:0.1:1400;

[xq,vq] = resampleSpectrum(x,v,xq);

plot(x,v,'o',xq,vq,':.');


%%

wavelengthRange   = [350 800]; % nm
wavelengthStep    = 1; % nm
wavelengthBleuing = 100; % nm


wavelength = wavelengthRange(1):wavelengthStep:(wavelengthRange(2)+wavelengthBleuing);

[~,emission] = resampleSpectrum(T_qdot655.wavelength,T_qdot655.em,wavelength);
[~,excitation] = resampleSpectrum(T_qdot655.wavelength,T_qdot655.ex,wavelength);
[~,emfilter] = resampleSpectrum(T_emfilter.wavelength,T_emfilter.transmission,wavelength);
[~,dichroic] = resampleSpectrum(T_dichroic.wavelength,T_dichroic.transmission,wavelength);
[~,dichroicLobe3] = resampleSpectrum(T_dichroiclobe3.wavelength,T_dichroiclobe3.transmission,wavelength);
[~,qe] = resampleSpectrum(T_qe.wavelength,T_qe.qe,wavelength);

systemTransmission = emfilter.*dichroic.*qe;
systemTransmission = systemTransmission';

totalFluorescence = nan(numel(wavelengthBleuing),1);
absorption488 = nan(numel(wavelengthBleuing),1);
absorption638 = nan(numel(wavelengthBleuing),1);

lobe488 = nan(numel(wavelengthBleuing),1);
lobe638 = nan(numel(wavelengthBleuing),1);
lobe3   = nan(numel(wavelengthBleuing),1);

figure('Position',[200 200 1000 200])
for i=1:wavelengthBleuing
   
    emission = blueArrayStep(emission);
    excitation = blueArrayStep(excitation);
    
    subplot(1,2,1)
    plot(wavelength,excitation,'Color',col_638,'LineWidth',1.5); hold on
    plot(wavelength,dichroic,'Color',0.8*[1 1 1])
    
    [~,idx_488] = min(abs(wavelength - 488));
    [~,idx_638] = min(abs(wavelength - 638));
    absorption488(i) = excitation(idx_488);
    absorption638(i) = excitation(idx_638);
    line([488 488],[0 1],'LineWidth',2,'Color',col_488); hold on
    line([638 638],[0 1],'LineWidth',2,'Color',col_638); hold on

    xlim(wavelengthRange); ylim([0 1]); box off
    set(gca,'Layer','top')

    subplot(1,2,2)
    plot(wavelength,emission,'Color',col_638,'LineWidth',1.5); hold on
    plot(wavelength,emfilter,'k')
    plot(wavelength,dichroic,'Color',0.8*[1 1 1])
    plot(wavelength,qe,'--k')
    % plot(wavelength,systemTransmission,'b')
    fill(wavelength,systemTransmission.*emission,'r','FaceAlpha',0.4,'EdgeColor','none')

    xlim(wavelengthRange); ylim([0 1]); box off
    set(gca,'Layer','top')
    pause(0.001)
    

    totalFluorescence(i) = sum(systemTransmission.*emission,'all');
    
    lobe3(i) = sum(systemTransmission.*(1-dichroicLobe3).*emission,'all');
    
    lobe488(i) = abs(totalFluorescence(i) - lobe3(i))*absorption488(i)/(absorption488(i) + absorption638(i));
    lobe638(i) = abs(totalFluorescence(i) - lobe3(i))*absorption638(i)/(absorption488(i) + absorption638(i));
    
    if i ~= wavelengthBleuing; clf; end
end

%%


figure
plot(1:wavelengthBleuing,lobe488); hold on
plot(1:wavelengthBleuing,lobe638);
plot(1:wavelengthBleuing,lobe3);
ylim([0 inf])

%% Function

function X_new = blueArrayStep(X)
X_new = X(:);
X_new = X_new(2:end);
X_new = [X_new; X_new(end)];

end

function [xq,vq] = resampleSpectrum(x,v,xq)
vq = interp1(x,v,xq,'linear');
vq(xq < x(1)) = v(1);
vq(xq > x(end)) = v(end);
end