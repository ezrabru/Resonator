clear all
close all
clc
addpath('lib')

% microscope/simulations parameters
detectionType = 'standard';
NA            = 1.45;
fTubeLens     = 200e-3;
magnification = (80/50)*100;
fObj          = fTubeLens/magnification;
nImmersion    = 1.518;
nMedium       = 1.33;
wavelength    = 650e-9;
numPixBFP     = 100;
typeNormalisation = 'individual';
fov = 101;

% camera parameters
cameraType          = 'nonoise'; 
camPixSize          = 16e-6;
noiseParam.offset   = 100;
noiseParam.QE       = 0.72;
noiseParam.sigmaReadNoise = 2.1;
noiseParam.eADU     = 1;
bit                 = 12;

% resonator parameters
resonatorParams.scanAmplitude = 37; % in pixels
resonatorParams.freqMirror    = 1000; % in Hz
resonatorParams.timeRes       = 100; % timepoint resolution per oscillation
resonatorParams.pulseLength1  = 104*1e-6; % duration of a pulse laser 1
resonatorParams.pulseLength2  = 104*1e-6; % duration of a pulse laser 2

wavelengthLaser1 = 638*1e-9; % m, wavelength laser 1
wavelengthLaser2 = 638*1e-9; % m, wavelength laser 2
powerLaser1      = 200; % W/cm2, laser power density lobe 1
powerLaser2      = 200; % W/cm2, laser power density lobe 2

QY1 = 0.7;
QY2 = 0.7;
extCoeffLaser1 = 300000;
extCoeffLaser2 = 300000;

bkgnd = 0;

%% Start with stationary mirror

resonatorParams.scanAmplitude = 0; % in pixels
resonatorParams.timeRes       = 1; % timepoint resolution per oscillation

Scope = Microscope(detectionType,NA,fObj,fTubeLens,nImmersion,nMedium,wavelength,numPixBFP);
Cam = Camera(cameraType,camPixSize,noiseParam,bit);
Sim = Simulation(Scope,Cam,fov,typeNormalisation,resonatorParams);

[psf_lobe1,psf_lobe2,stackLobe12] = Sim.getPSF_2lobes;

delaytime = 50; % ms
numFrames = 30;
photons = 3000;
bkgnd = 2;

% camera parameters
cameraType          = 'sCMOS'; 
camPixSize          = 16e-6;
noiseParam.offset   = 100;
noiseParam.QE       = 0.95;
noiseParam.sigmaReadNoise = 1;
noiseParam.eADU     = 1;
bit                 = 12;
Cam = Camera(cameraType,camPixSize,noiseParam,bit);

filename = 'stationary_emitter.gif';
% normalise 
stackLobe12 = photons*stackLobe12/sum(stackLobe12(:));

stackLobe12 = stackLobe12 + bkgnd;
for n = 1:numFrames

    % add noise to image
    newFrame = stackLobe12;
    newFrame = Cam.cameraNoiseModel(newFrame);
    newFrame = kron(newFrame,ones(10));
    
    % rescale to 8-bit for writing away as a gif
    cMin = 90;
    cMax = 250;
    newFrame(newFrame < cMin) = cMin;
    newFrame(newFrame > cMax) = cMax;
    newFrame = 255*(newFrame - cMin)/(cMax - cMin);

    if n == 1
        imwrite(newFrame,filename,'gif', 'Loopcount',inf,'DelayTime',delaytime/1000);
    else
        imwrite(newFrame,filename,'gif','WriteMode','append','DelayTime',delaytime/1000);
    end
end