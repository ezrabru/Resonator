clear all
close all
clc

directory_out = fullfile(pwd,'figures experiment spectra');
if ~exist(directory_out,'dir'); mkdir(directory_out); end

fontsize = 8;

col_lightgray = 0.8*[1 1 1];
col_405 = [135 90 175]/255; % 405 nm
col_488 = [2 100 173]/255; % 488 nm
col_561 = [177 212 55]/255; % 561 nm
col_638 = [212 0 0]/255; % 638 nm

qdot585 = readmatrix(fullfile(pwd,'spectra dyes and filters','dyes','Qdot585.txt'));
T_qdot585.wavelength = qdot585(:,1);
T_qdot585.ex = qdot585(:,2)/100;
T_qdot585.em = qdot585(:,3)/100;

qdot655 = readmatrix(fullfile(pwd,'spectra dyes and filters','dyes','Qdot655.txt'));
T_qdot655.wavelength = qdot655(:,1);
T_qdot655.ex = qdot655(:,2)/100;
T_qdot655.em = qdot655(:,3)/100;

exfilter = readmatrix(fullfile(pwd,'spectra dyes and filters','filters','ZET405-488-561-640xv2.txt'));
T_exfilter.wavelength = exfilter(:,1);
T_exfilter.transmission = exfilter(:,2);

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

%% Plot data in light mode

fig = figure('Position',[600 300 850 200]);
subplot(1,2,1)
% line([405 405],[0 1],'LineWidth',2,'Color',col_405); hold on
line([488 488],[0 1],'LineWidth',2,'Color',col_488); hold on
% line([561 561],[0 1],'LineWidth',2,'Color',col_561); hold on
line([638 638],[0 1],'LineWidth',2,'Color',col_638); hold on
plot(T_qdot655.wavelength,T_qdot655.ex/0.59,'Color',col_638,'LineWidth',2);
plot(T_exfilter.wavelength,T_exfilter.transmission,'k');
plot(T_dichroic.wavelength,T_dichroic.transmission,'Color',col_lightgray);
xlim([350 750]); ylim([0 1]); xticks(300:50:800)
xlabel('Wavelength (nm)'); box off
set(gca,'FontSize',fontsize,'Layer','top')

subplot(1,2,2)
plot(T_qdot655.wavelength,T_qdot655.em,'Color',col_638,'LineWidth',2); hold on
plot(T_emfilter.wavelength,T_emfilter.transmission,'k');
plot(T_dichroic.wavelength,T_dichroic.transmission,'Color',col_lightgray);
plot(T_dichroiclobe3.wavelength,T_dichroiclobe3.transmission,'b');
plot(T_qe.wavelength,T_qe.qe,'--k');
xlim([350 750]); ylim([0 1]); xticks(300:50:800)
xlabel('Wavelength (nm)'); box off
set(gca,'FontSize',fontsize,'Layer','top')

savefig(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_lightMode.fig'))
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_lightMode.png'),'Resolution',400)
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_lightMode.eps'))

%% Plot data in dark mode

fig = figure('Position',[600 300 850 200]);
set(gcf,'Color','k')

subplot(1,2,1)
% line([405 405],[0 1],'LineWidth',2,'Color',col_405); hold on
line([488 488],[0 1],'LineWidth',2,'Color',col_488); hold on
% line([561 561],[0 1],'LineWidth',2,'Color',col_561); hold on
line([638 638],[0 1],'LineWidth',2,'Color',col_638); hold on
plot(T_qdot655.wavelength,T_qdot655.ex/0.59,'Color',col_638,'LineWidth',2);
plot(T_exfilter.wavelength,T_exfilter.transmission,'w');
plot(T_dichroic.wavelength,T_dichroic.transmission,'Color',0.4*col_lightgray);
xlim([350 750]); ylim([0 1]); xticks(300:50:800)
xlabel('Wavelength (nm)'); box off
set(gca,'FontSize',fontsize,'XColor','w','YColor','w','Color','k','Layer','top')

subplot(1,2,2)
plot(T_qdot655.wavelength,T_qdot655.em,'Color',col_638,'LineWidth',2); hold on
plot(T_emfilter.wavelength,T_emfilter.transmission,'w');
plot(T_dichroic.wavelength,T_dichroic.transmission,'Color',0.4*col_lightgray);
plot(T_dichroiclobe3.wavelength,T_dichroiclobe3.transmission,'g');
plot(T_qe.wavelength,T_qe.qe,'--w');
xlim([350 750]); ylim([0 1]); xticks(300:50:800)
xlabel('Wavelength (nm)'); box off
set(gca,'FontSize',fontsize,'XColor','w','YColor','w','Color','k','Layer','top')

savefig(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_darkMode.fig'))
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_darkMode.png'),'Resolution',400,'BackgroundColor','k')
set(gcf,'renderer','Painters')
exportgraphics(fig,fullfile(directory_out,'experimentSpectra_3lobe_QDot655_darkMode.eps'),'BackgroundColor','k')