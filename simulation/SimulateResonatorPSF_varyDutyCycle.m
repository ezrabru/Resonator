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
resonatorParams.freqMirror   = 1000; % in Hz
resonatorParams.timeRes      = 100; % timepoint resolution per oscillation
dutyCycle = [0:5:100, 100:-5:0];

% normalise intensity
normaliseIntensity = 1;

numFrames = 1;


%% Simulate 

% Get normalisation from a 100 duty cycle PSF
Scope = Microscope(detectionType,NA,fObj,fTubeLens,nImmersion,nMedium,wavelength,numPixBFP);
resonatorParams.pulseLength1 = 0.5*(1/resonatorParams.freqMirror); % duration of a pulse laser 1
resonatorParams.pulseLength2 = 0.5*(1/resonatorParams.freqMirror); % duration of a pulse laser 2
Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);
[psf_lobe1,psf_lobe2,~] = Sim.getPSF_2lobes;
full_psf = psf_lobe1 + psf_lobe2;
normalisationMax = max(full_psf(:));

%% Simulate with temporal binning

if normaliseIntensity
    filename = 'varyingDutyCycle_normalised.gif';
else
    filename = 'varyingDutyCycle.gif';
end
outputpath = fullfile(pwd,'results',filename);


for id_dc = 1:numel(dutyCycle)
    fprintf('Simulating duty cycle %d/%d\n',id_dc,numel(dutyCycle))
    
    pulseLength = dutyCycle(id_dc)/(2*resonatorParams.freqMirror*100);

    resonatorParams.pulseLength1 = pulseLength; % duration of a pulse laser 1
    resonatorParams.pulseLength2 = pulseLength; % duration of a pulse laser 2
    Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

    [psf_lobe1,psf_lobe2,~] = Sim.getPSF_2lobes;
    frame = psf_lobe1 + psf_lobe2;

    if normaliseIntensity
        frame = 255*frame/normalisationMax;
    else
        frame = 255*frame/max(frame(:));
    end
    
    for i=1:numFrames
        frame_i = kron(frame,ones(10));
        imagesc(frame_i); pause(0.1)
        
        if (i == 1) && (id_dc == 1)
            imwrite(frame_i,outputpath,'gif', 'Loopcount',inf,'DelayTime',1/30);
        else
            imwrite(frame_i,outputpath,'gif','WriteMode','append','DelayTime',1/30);
        end

    end

end