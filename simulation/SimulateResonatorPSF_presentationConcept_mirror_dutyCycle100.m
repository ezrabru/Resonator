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
CamNoNoise = Camera(cameraType,camPixSize,noiseParam,bit);
cameraType          = 'sCMOS'; 
camPixSize          = 16e-6;
noiseParam.offset   = 100;
noiseParam.QE       = 0.95;
noiseParam.sigmaReadNoise = 1;
noiseParam.eADU     = 1;
bit                 = 12;
CamNoise = Camera(cameraType,camPixSize,noiseParam,bit);

% resonator parameters
resonatorParams.scanAmplitude = 15; % in pixels
resonatorParams.freqMirror    = 1000; % in Hz
resonatorParams.timeRes       = 30; % timepoint resolution per oscillation
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

photons = 10000;
bkgnd = 2;

cMin = 90;
cMax = 250;

%% 0.5 Hz oscillation

Scope = Microscope(detectionType,NA,fObj,fTubeLens,nImmersion,nMedium,wavelength,numPixBFP);
Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

[~,~,stackLobe12] = Sim.getPSF_2lobes;

animationFrequency = 0.5; % Hz

filename = 'oscillating_emitter_0.5Hz.gif';

% normalise 
normalisation = sum(stackLobe12(:,:,1),'all');
stackLobe12 = photons*stackLobe12/normalisation;
stackLobe12 = stackLobe12 + bkgnd;

stackLobe12reverse = flip(stackLobe12,3);
stackLobe12 = cat(3,stackLobe12,stackLobe12reverse);

numFrames = size(stackLobe12,3);
for n = 1:numFrames

    % add noise to image
    newFrame = CamNoise.cameraNoiseModel(stackLobe12(:,:,n));
    newFrame = kron(newFrame,ones(10));
    
    % rescale to 8-bit for writing away as a gif
    newFrame(newFrame < cMin) = cMin;
    newFrame(newFrame > cMax) = cMax;
    newFrame = 255*(newFrame - cMin)/(cMax - cMin);

    if n == 1
        imwrite(newFrame,filename,'gif', 'Loopcount',inf,'DelayTime',(1/animationFrequency)/numFrames);
    else
        imwrite(newFrame,filename,'gif','WriteMode','append','DelayTime',(1/animationFrequency)/numFrames);
    end
end

%% 1 Hz oscillation

[~,~,stackLobe12] = Sim.getPSF_2lobes;

animationFrequency = 1; % Hz

filename = 'oscillating_emitter_1Hz.gif';

% normalise 
normalisation = sum(stackLobe12(:,:,1),'all');
stackLobe12 = photons*stackLobe12/normalisation;
stackLobe12 = stackLobe12 + bkgnd;

stackLobe12reverse = flip(stackLobe12,3);
stackLobe12 = cat(3,stackLobe12,stackLobe12reverse);

numFrames = size(stackLobe12,3);
for n = 1:numFrames

    % add noise to image
    newFrame = CamNoise.cameraNoiseModel(stackLobe12(:,:,n));
    newFrame = kron(newFrame,ones(10));
    
    % rescale to 8-bit for writing away as a gif
    newFrame(newFrame < cMin) = cMin;
    newFrame(newFrame > cMax) = cMax;
    newFrame = 255*(newFrame - cMin)/(cMax - cMin);

    if n == 1
        imwrite(newFrame,filename,'gif', 'Loopcount',inf,'DelayTime',(1/animationFrequency)/numFrames);
    else
        imwrite(newFrame,filename,'gif','WriteMode','append','DelayTime',(1/animationFrequency)/numFrames);
    end
end

%% 1000 Hz oscillation

[psf_lobe1,psf_lobe2,stackLobe12] = Sim.getPSF_2lobes;

animationFrequency = 0.5; % Hz

filename = 'oscillating_emitter_1000Hz.gif';

% normalise 
stackLobe12 = sum(stackLobe12,3);
normalisation = sum(stackLobe12(:));
stackLobe12 = photons*stackLobe12/normalisation;
stackLobe12 = stackLobe12 + bkgnd;

numFrames = 2*resonatorParams.timeRes;
for n = 1:numFrames

    % add noise to image
    newFrame = CamNoise.cameraNoiseModel(stackLobe12);
    newFrame = kron(newFrame,ones(10));
    
    % rescale to 8-bit for writing away as a gif
    newFrame(newFrame < cMin) = cMin;
    newFrame(newFrame > cMax) = cMax;
    newFrame = 255*(newFrame - cMin)/(cMax - cMin);

    if n == 1
        imwrite(newFrame,filename,'gif', 'Loopcount',inf,'DelayTime',(1/animationFrequency)/numFrames);
    else
        imwrite(newFrame,filename,'gif','WriteMode','append','DelayTime',(1/animationFrequency)/numFrames);
    end
end


%% Stationary mirror

cameraType          = 'nonoise'; 
camPixSize          = 16e-6;
noiseParam.offset   = 100;
noiseParam.QE       = 0.95;
noiseParam.sigmaReadNoise = 1;
noiseParam.eADU     = 1;
bit                 = 12;
Cam = Camera(cameraType,camPixSize,noiseParam,bit);

resonatorParams.scanAmplitude = 0; % in pixels
Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

[~,~,stackLobe12] = Sim.getPSF_2lobes;
stackLobe12 = stackLobe12(:,:,1);

numFrames = 2*resonatorParams.timeRes;

filename = 'stationary_emitter.gif';
% normalise 
stackLobe12 = photons*stackLobe12/sum(stackLobe12(:));

stackLobe12 = stackLobe12 + bkgnd;
for n = 1:numFrames

    % add noise to image
    newFrame = stackLobe12;
    newFrame = CamNoise.cameraNoiseModel(newFrame);
    newFrame = kron(newFrame,ones(10));
    
    % rescale to 8-bit for writing away as a gif
    newFrame(newFrame < cMin) = cMin;
    newFrame(newFrame > cMax) = cMax;
    newFrame = 255*(newFrame - cMin)/(cMax - cMin);

    if n == 1
        imwrite(newFrame,filename,'gif', 'Loopcount',inf,'DelayTime',(1/animationFrequency)/numFrames);
    else
        imwrite(newFrame,filename,'gif','WriteMode','append','DelayTime',(1/animationFrequency)/numFrames);
    end
end
