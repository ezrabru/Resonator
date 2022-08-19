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

photons = 5000;
bkgnd = 2;

cMin = 90;
cMax = 250;

% binning = 1 means 1 Hz oscillations
binning = [1,2,5,10,20,50,100];


%% Simulate 

Scope = Microscope(detectionType,NA,fObj,fTubeLens,nImmersion,nMedium,wavelength,numPixBFP);
Sim = Simulation(Scope,CamNoNoise,fov,typeNormalisation,resonatorParams);

[~,~,stackLobe12] = Sim.getPSF_2lobes;
stackLobe12reverse = flip(stackLobe12,3);
stackLobe12 = cat(3,stackLobe12,stackLobe12reverse);
stackLobe12 = rot90(stackLobe12);


%% Simulate with temporal binning

for id_bin = 1:numel(binning)
    filename = sprintf('oscillating_emitter_%dHz.gif',binning(id_bin));
    
    numFrames = 2*resonatorParams.timeRes;
    currentFrame = 0;
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
        
        currentFrame = currentFrame + binning(id_bin);
        
        frame = bkgnd + photons*frame/sum(frame(:));
        frame = CamNoise.cameraNoiseModel(frame);

        % rescale to 8-bit for writing away as a gif
        cMin = min(frame(:));
        cMax = max(frame(:));
        frame(frame < cMin) = cMin;
        frame(frame > cMax) = cMax;
        frame = 255*(frame - cMin)/(cMax - cMin);
        frame = round(frame);
    
        frame = kron(frame,ones(10));
    
        if i == 1
            imwrite(frame,filename,'gif', 'Loopcount',inf,'DelayTime',1/(2*resonatorParams.timeRes));
        else
            imwrite(frame,filename,'gif','WriteMode','append','DelayTime',1/(2*resonatorParams.timeRes));
        end
    end
end