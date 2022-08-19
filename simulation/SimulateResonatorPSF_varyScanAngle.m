clear all; close all; clc
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
cameraType          = 'sCMOS'; 
camPixSize          = 16e-6;
noiseParam.offset   = 100;
noiseParam.QE       = 0.95;
noiseParam.sigmaReadNoise = 0.9;
noiseParam.eADU     = 1;
bit                 = 16;
CamNoise = Camera(cameraType,camPixSize,noiseParam,bit);
CamNoNoise = Camera('nonoise',camPixSize,noiseParam,bit);

% resonator parameters
scanAmplitude = [0:1:38, 38:-1:0]; % in pixels
resonatorParams.freqMirror   = 1000; % in Hz
resonatorParams.timeRes      = 100; % timepoint resolution per oscillation
resonatorParams.pulseLength1 = 104*1e-6; % duration of a pulse laser 1
resonatorParams.pulseLength2 = 104*1e-6; % duration of a pulse laser 2

% normalise intensity
normaliseIntensity = 1;

numFrames = 1;


%% Simulate 

Scope = Microscope(detectionType,NA,fObj,fTubeLens,nImmersion,nMedium,wavelength,numPixBFP);
resonatorParams.scanAmplitude = 0; % in pixels
Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

[~,~,stackLobe12] = Sim.getPSF_2lobes;
stackLobe12reverse = flip(stackLobe12,3);
stackLobe12 = cat(3,stackLobe12,stackLobe12reverse);
stackLobe12 = rot90(stackLobe12);

% normalisation
frame = sum(stackLobe12,3);
frame = frame/sum(frame(:));
normalisationMax = max(frame(:));


%% Simulate with temporal binning

if normaliseIntensity
    filename = 'varyingScanAngle_normalised.gif';
else
    filename = 'varyingScanAngle.gif';
end
outputpath = fullfile(pwd,'results',filename);


for id_angle = 1:numel(scanAmplitude)
    fprintf('Simulating scan amplitude %d/%d\n',id_angle,numel(scanAmplitude))
        
    resonatorParams.scanAmplitude = scanAmplitude(id_angle);
    Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

    [~,~,stackLobe12] = Sim.getPSF_2lobes;
    frame = sum(stackLobe12,3);

    if normaliseIntensity
        frame = frame/sum(frame(:));
        frame = 255*frame/normalisationMax;
    else
        frame = 255*frame/max(frame(:));
    end
    
    for i=1:numFrames
        frame_i = kron(frame,ones(10));
        imagesc(frame_i); pause(0.1)
        
        if (i == 1) && (id_angle == 1)
            imwrite(frame_i,outputpath,'gif', 'Loopcount',inf,'DelayTime',1/30);
        else
            imwrite(frame_i,outputpath,'gif','WriteMode','append','DelayTime',1/30);
        end

    end

end