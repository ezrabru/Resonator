classdef Camera
    %CAMERA object contains camera properties and models.
    
    properties
        cameraType % camera type
        noiseParam % noise parameters (struct)
        camPixSize % size camera pixels
        bit % dynamic range (e.g. 12 for CS505MUP, Thorlabs)
    end
    
    methods
        
        function obj = Camera(cameraType,camPixSize,noiseParam,bit)
            %CAMERA Construct an instance of this class
            obj.cameraType = cameraType;
            obj.camPixSize = camPixSize;
            obj.noiseParam = noiseParam;
            obj.bit = bit;
        end
        
        function img = cameraNoiseModel(obj,img)
            % CAMERANOISEMODEL
            switch obj.cameraType
                case 'nonoise'
                    % do nothing
                case 'ideal'
                    img = poissrnd(img);
                case 'sCMOS'
                    img = noiseModelsCMOS(obj,img);
                case 'EMCCD'
                    img = noiseModelEMCCD(obj,img);
                otherwise
                    fprintf('Unexpted value for property "cameraType" in class "Camera".\n')
            end
        end
        
        function img = noiseModelsCMOS(obj,img)
            % NOISEMODELSCMOS Model detection process by a sCMOS camera
            % img: image in units of photons
            % QE: quantum efficiency of the sensor at average wavelength of the signal (electrons/photon)
            % offset: camera offset, bias or baseline (ADU, analog-to-digital unit)
            % e_adu: analog-to-digital conversion factor (electrons per analog-to-digital unit)
            % sigmaReadNoise: read noise in electrons (root mean square read noise)
            % bit: resolution of the analog-to-digital converter (saturation = 2^bit - 1)
            % 
            % Example parameters of the Prime 95B from Photometrics
            % QE = 0.95, offset = 100, e_adu = 1, sigmaReadNoise = 1.3
            % 
            % returns: original image with detection noise added
            
            offset         = obj.noiseParam.offset; % camera offset
            QE             = obj.noiseParam.QE; % quantum efficiency
            sigmaReadNoise = obj.noiseParam.sigmaReadNoise; % the standard deviation (or RMS) of the read noise
            eADU           = obj.noiseParam.eADU; % analog to digital conversion factor
            
            photoElectrons = poissrnd(QE*img) + normrnd(0,sigmaReadNoise,size(img));
            img = round(photoElectrons/eADU) + offset;
            saturation = 2^obj.bit - 1;
            img(img > saturation) = saturation;
        end

        function img = noiseModelEMCCD(obj,img)
            % NOISEMODELEMCCD Model detection process by an EMCCD camera
            % img: image in units of photons
            % QE: quantum efficiency of the sensor at average wavelength of the signal (electrons/photon)
            % offset: camera offset, bias or baseline (ADU, analog-to-digital unit)
            % e_adu: analog-to-digital conversion factor (electrons per analog-to-digital unit)
            % sigmaReadNoise: read noise in electrons (root mean square read noise)
            % emgain: electron-multiplying gain
            % c: spurious charge (clock-induced charge only, dark counts negligible)
            % bit: resolution of the analog-to-digital converter (saturation = 2^bit - 1)
            %
            % Example parameters of the Evolve Delta 512 from Photometrics
            % QE = 0.9, offset = 100, e_adu = 45, sigmaReadNoise = 74.4,emgain = 300, c = 0.002
            %
            % returns: original image with detection noise added

            offset         = obj.noiseParam.offset; % camera offset
            QE             = obj.noiseParam.QE;
            sigmaReadNoise = obj.noiseParam.sigmaReadNoise;
            eADU           = obj.noiseParam.eADU;
            emgain         = obj.noiseParam.emgain;
            c              = obj.noiseParam.c;

            photoElectrons = poissrnd(QE*img + c);
            photoElectrons = gamrnd(photoElectrons,emgain,size(img,1),size(img,2)) + normrnd(0,sigmaReadNoise,size(img));
            img = floor(photoElectrons/eADU) + offset;
            saturation = 2^obj.bit - 1;
            img(img > saturation) = saturation;
        end
        
    end
end
