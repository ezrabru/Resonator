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
scanAmplitude = 0:0.5:38; % in pixels
resonatorParams.scanAmplitude = 0; % in pixels
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

photons = 5000;
bkgnd = 2;

cMin = 90;
cMax = 250;

% normalise intensity
normaliseIntensity = 0;

numFrames = 5;


%% Simulate 

Scope = Microscope(detectionType,NA,fObj,fTubeLens,nImmersion,nMedium,wavelength,numPixBFP);
Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

[~,~,stackLobe12] = Sim.getPSF_2lobes;
stackLobe12reverse = flip(stackLobe12,3);
stackLobe12 = cat(3,stackLobe12,stackLobe12reverse);
stackLobe12 = rot90(stackLobe12);

% normalisation
normalisationVol = sum(stackLobe12(:,:,1),'all');
normalisedFrame = stackLobe12(:,:,1)/normalisationVol;
normalisationMax = max(normalisedFrame(:));


%% Simulate with temporal binning

fullPSF = sum(stackLobe12,3);

if normaliseIntensity
    filename = sprintf('varyingScanAngle_oscillation_normalised.gif');
else
    filename = sprintf('varyingScanAngle_oscillation.gif');
end

for id_angle = 1:numel(scanAmplitude)
    
    resonatorParams.scanAmplitude = scanAmplitude(id_angle);
    Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

    [~,~,stackLobe12] = Sim.getPSF_2lobes;
    frame = sum(stackLobe12,3);

    if normaliseIntensity
        frame = frame/normalisationVol;
        frame = 255*frame/normalisationMax;
    else
        frame = 255*frame/max(frame(:));
    end

    for i=1:numFrames
        
        % (optional add noise)

        if (i == 1) && (id_angle == 1)
            imwrite(kron(frame,ones(10)),filename,'gif', 'Loopcount',inf,'DelayTime',1/(2*resonatorParams.timeRes));
        else
            imwrite(kron(frame,ones(10)),filename,'gif','WriteMode','append','DelayTime',1/(2*resonatorParams.timeRes));
        end
        
    end
end