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
resonatorParams.scanAmplitude = 37; % in pixels
resonatorParams.freqMirror    = 1000; % in Hz
resonatorParams.timeRes       = 100; % timepoint resolution per oscillation
resonatorParams.pulseLength1  = 104*1e-6; % duration of a pulse laser 1
resonatorParams.pulseLength2  = 104*1e-6; % duration of a pulse laser 2

% binning = 1 means 1 Hz oscillations
binning = [1,2,3,4,5,10,20,50,100,1000];

% normalise intensity
normaliseIntensity = 0;

%% Simulate 

Scope = Microscope(detectionType,NA,fObj,fTubeLens,nImmersion,nMedium,wavelength,numPixBFP);
Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

[~,~,stackLobe12] = Sim.getPSF_2lobes;
stackLobe12reverse = flip(stackLobe12,3);
stackLobe12 = cat(3,stackLobe12,stackLobe12reverse);
stackLobe12 = rot90(stackLobe12);

% normalisation
frame = stackLobe12(:,:,1);
frame = frame/sum(frame(:));
normalisationMax = max(frame(:));

%% Simulate with temporal binning

fullPSF = sum(stackLobe12,3);

if normaliseIntensity
    filename = sprintf('acceleratingOscillation_aascanAngle%d_normalised.gif',round(resonatorParams.scanAmplitude));
else
    filename = sprintf('acceleratingOscillation_aascanAngle%d.gif',round(resonatorParams.scanAmplitude));
end
outputpath = fullfile(pwd,'results',filename);

currentFrame = 0;
for id_bin = 1:numel(binning)
    
    %numFrames = 2*resonatorParams.timeRes;
    if  id_bin == numel(binning)
        numFrames = 500;
    else
        numFrames = 100;
    end

    for i=1:numFrames

        idx_start = mod(currentFrame,2*resonatorParams.timeRes)+1;
        idx_end = mod(currentFrame+binning(id_bin),2*resonatorParams.timeRes)+1;

        if idx_start > idx_end
            substack1 = stackLobe12(:,:,idx_start:end);
            substack2 = stackLobe12(:,:,1:idx_end);
            frame = sum(substack1,3) + sum(substack2,3);
        else
            substack = stackLobe12(:,:,idx_start:idx_end);
            frame = sum(substack,3);
        end

        if binning(id_bin) > 2*resonatorParams.timeRes
            frame = frame + floor(binning(id_bin)/(2*resonatorParams.timeRes))*fullPSF;
        end
        
        if normaliseIntensity
            frame = frame/sum(frame(:));
            frame = 255*frame/normalisationMax;
        else
            frame = 255*frame/max(frame(:));
        end
        frame = kron(frame,ones(10));
    
        if (i == 1) && (id_bin == 1)
            imwrite(frame,outputpath,'gif', 'Loopcount',inf,'DelayTime',1/(2*resonatorParams.timeRes));
        else
            imwrite(frame,outputpath,'gif','WriteMode','append','DelayTime',1/(2*resonatorParams.timeRes));
        end

        currentFrame = currentFrame + binning(id_bin);
    end
end