classdef Simulation
    %SIMULATION object contains properties and methods to simulate a
    %dataset given a localisation list.
    
    properties
        Microscope % instance of the Microscope class
        Camera % instance of the Camera class
        fov % field of view (in pixels)
        typeNormalisation % 'relativeToStandardPSF' or 'individual'
        resonatorParams % struct with resonator parameters
        
        % contents of object 'Microscope'
        detectionType % 'standard'
        NA % numerical aperture of the objective
        fObj % focal length of the objective lens (fObj = M/fTubeLens)
        fTubeLens % focal length of the tube lens
        magnification % system magnification
        nImmersion % refractive index of the cover glass
        nMedium % refractive index of the medium in which the molecule has been embedded
        wavelength % wavelength of the emission
        numPixBFP % number of pixels in the BFP (for simulation purposes)
        rhoMax % maximum collection angle of the objective lens
        bfpRadiusNative % radius of the native back focal plane
        bfpDiameterNative % diameter of the native back focal plane
        BFPphi
        BFPrho
        BFPaperture

        pad % calculated padding before fft2 in paraxial image formation to achieve correct virtual pixel size
        
        k0 % wavenumber
        GreensTensorZ0_bfp   % Z0 = for a molecule in focus (i.e. on the refractive index interface if there is one)
        GreensTensorZ0_img   % Z0 = for a molecule in focus (i.e. on the refractive index interface if there is one)
        BasisFunctionsZ0_bfp % Z0 = for a molecule in focus (i.e. on the refractive index interface if there is one)
        BasisFunctionsZ0_img % Z0 = for a molecule in focus (i.e. on the refractive index interface if there is one)
        
        fontsize % fontsize used for any plotting methods in this class
    end
    
    methods
        
        function obj = Simulation(Microscope,Camera,fov,typeNormalisation,resonatorParams)
            %SIMULATION Construct an instance of this class
            obj.Microscope        = Microscope;
            obj.Camera            = Camera;
            obj.fov               = fov;
            obj.typeNormalisation = typeNormalisation;
            obj.resonatorParams   = resonatorParams;

            % copy contents of object 'Microscope' for more readable code
            obj.detectionType     = obj.Microscope.detectionType;
            obj.NA                = obj.Microscope.NA;
            obj.fObj              = obj.Microscope.fObj;
            obj.fTubeLens         = obj.Microscope.fTubeLens;
            obj.magnification     = obj.fTubeLens/obj.fObj;
            obj.nImmersion        = obj.Microscope.nImmersion;
            obj.nMedium           = obj.Microscope.nMedium;
            obj.wavelength        = obj.Microscope.wavelength;
            obj.numPixBFP         = obj.Microscope.numPixBFP;
            obj.rhoMax            = obj.Microscope.rhoMax;
            obj.bfpRadiusNative   = obj.Microscope.bfpRadiusNative;
            obj.bfpDiameterNative = obj.Microscope.bfpDiameterNative;
            obj.BFPphi            = obj.Microscope.BFPphi;
            obj.BFPrho            = obj.Microscope.BFPrho;
            obj.BFPaperture       = obj.Microscope.BFPapertureNative;
            obj.k0                = obj.Microscope.k0;
            
            obj.pad = obj.calculatePadding;
            
            % at (x,y,z,Delta) = (0,0,0,0);
            G = obj.getGreensTensor_bfp(0);
            g = obj.getGreensTensor_img(0,0,0,0);
            B_bfp = obj.getBasisFunctions_bfp(0);
            B_img = obj.getBasisFunctions_img(0,0,0,0);
            obj.GreensTensorZ0_bfp = G;
            obj.GreensTensorZ0_img = g;
            obj.BasisFunctionsZ0_bfp = B_bfp;
            obj.BasisFunctionsZ0_img = B_img;
            
            obj.fontsize = 9;
        end
        
        function [psf_lobe1,psf_lobe2,stackLobe12] = getPSF_2lobes(obj)
            %GETPSF_2LOBES

            timeRes       = obj.resonatorParams.timeRes;
            scanAmplitude = obj.resonatorParams.scanAmplitude;
            freqMirror    = obj.resonatorParams.freqMirror;
            pulseLength1  = obj.resonatorParams.pulseLength1;
            pulseLength2  = obj.resonatorParams.pulseLength2;

            pixelsize = obj.Camera.camPixSize/obj.magnification;
                        
            % Get path image movement
            t = linspace(-pi/2,pi/2,timeRes);
            signal = sin(t)';
            tPoints = numel(signal);
            locsTable = table;
            locsTable.x = zeros(tPoints,1);
            locsTable.y = pixelsize*(scanAmplitude/2)*signal;
            locsTable.z = zeros(tPoints,1);
            locsTable.frame = (1:tPoints)';
            locsTable.gamma = zeros(tPoints,1);
            locsTable.photons = ones(tPoints,1);
            locsTable.delta = zeros(tPoints,1);
            locsTable.phi = zeros(tPoints,1);
            locsTable.theta = zeros(tPoints,1);
            [stackLobe12,~] = obj.simulateSample(locsTable,0,0,1);
            
            % Generate blinking trace
            oscillationsPerExposure = 1;
            timePoints = timeRes*oscillationsPerExposure;
            blinkingTrace = ones(timePoints,1); % continuously "on", can be changed
            
            % Get resonant scanner output voltage to trigger laser
            signal = cos(pi/2 + linspace(0,2*pi*oscillationsPerExposure,timePoints));
            t = 1e3*linspace(0,pi,timePoints);
            
            % Generate laser pulse trains
            [~,laserTrace1,laserTrace2] = obj.getLaserPulseTrace(t,signal,pulseLength1*freqMirror,pulseLength2*freqMirror,0);
            
            % Multiply laser pulse trains with blinking and excitation efficiency
            responseLobe1 = blinkingTrace'.*laserTrace1;
            responseLobe2 = blinkingTrace'.*laserTrace2;
            
            psf_lobe1 = zeros(size(stackLobe12(:,:,1)));
            psf_lobe2 = zeros(size(stackLobe12(:,:,1)));
            for i=1:timePoints-1
                psf_lobe1 = psf_lobe1 + responseLobe1(i)*stackLobe12(:,:,i);
                psf_lobe2 = psf_lobe2 + responseLobe2(i)*stackLobe12(:,:,i);
            end

            % normalise lobes
            psf_lobe1 = psf_lobe1/sum(psf_lobe1(:));
            psf_lobe2 = psf_lobe2/sum(psf_lobe2(:));
            
        end

        function [psf_lobe1,psf_lobe2,psf_lobe3] = getPSF_3lobes(obj)
            %GETPSF_3LOBES

            timeRes       = obj.resonatorParams.timeRes;
            scanAmplitude = obj.resonatorParams.scanAmplitude;
            freqMirror    = obj.resonatorParams.freqMirror;
            pulseLength1  = obj.resonatorParams.pulseLength1;
            pulseLength2  = obj.resonatorParams.pulseLength2;
            
            pixelsize = obj.Camera.camPixSize/obj.magnification;
            dispLobe3 = pixelsize*sqrt(2*(scanAmplitude/2)^2);
            
            % get resonator lobes -----------------------------------------
            
            % Get path image movement
            t = linspace(-pi/2,pi/2,timeRes);
            signal = sin(t)';
            tPoints = numel(signal);
            locsTable = table;
            locsTable.x = zeros(tPoints,1) - dispLobe3/2;
            locsTable.y = pixelsize*(scanAmplitude/2)*signal;
            locsTable.z = zeros(tPoints,1);
            locsTable.frame = (1:tPoints)';
            locsTable.gamma = zeros(tPoints,1);
            locsTable.photons = ones(tPoints,1);
            locsTable.delta = zeros(tPoints,1);
            locsTable.phi = zeros(tPoints,1);
            locsTable.theta = zeros(tPoints,1);
            [stackLobe12,~] = obj.simulateSample(locsTable,0,0,1);
            
            % Generate blinking trace
            oscillationsPerExposure = 1;
            timePoints = timeRes*oscillationsPerExposure;
            blinkingTrace = ones(timePoints,1); % continuously "on", can be changed
            
            % Get resonant scanner output voltage to trigger laser on
            signal = cos(pi/2 + linspace(0,2*pi*oscillationsPerExposure,timePoints));
            t = 1e3*linspace(0,pi,timePoints);
            
            % Generate laser pulse trains
            [~,laserTrace1,laserTrace2] = obj.getLaserPulseTrace(t,signal,pulseLength1*freqMirror,pulseLength2*freqMirror,0);
            
            % Multiply laser pulse trains with blinking and excitation efficiency
            responseLobe1 = blinkingTrace'.*laserTrace1;
            responseLobe2 = blinkingTrace'.*laserTrace2;
            
            psf_lobe1 = zeros(size(stackLobe12(:,:,1)));
            psf_lobe2 = zeros(size(stackLobe12(:,:,1)));
            for i=1:timePoints-1
                psf_lobe1 = psf_lobe1 + responseLobe1(i)*stackLobe12(:,:,i);
                psf_lobe2 = psf_lobe2 + responseLobe2(i)*stackLobe12(:,:,i);
            end

            % get lobe 3 --------------------------------------------------
            locsTable = table;
            locsTable.x = dispLobe3/2;
            locsTable.y = 0;
            locsTable.z = 0;
            locsTable.frame = 1;
            locsTable.gamma = 0;
            locsTable.photons = 1;
            locsTable.delta = 0;
            locsTable.phi = 0;
            locsTable.theta = 0;
            [psf_lobe3,~] = obj.simulateSample(locsTable,0,0,0);
            
            % normalise lobes
            psf_lobe1 = psf_lobe1/sum(psf_lobe1(:));
            psf_lobe2 = psf_lobe2/sum(psf_lobe2(:));
            psf_lobe3 = psf_lobe3/sum(psf_lobe3(:));
            
        end

        function [simulations,filepath] = simulateSample(obj,locsTable,filepath,bkgnd,verbose)
            %SIMULATESAMPLE
            
            % get image dimensions
            nx = obj.fov;
            ny = obj.fov;
            
            % check if filepath is a string
            if ~ischar(filepath) % if it's not a string, don't write away results
                flag_write = 0;
                
                % initialise stack to contain simulated images
                simulations = nan(nx,ny,numel(unique(locsTable.frame)));
                
            else
                flag_write = 1; % write stack away
                
                % make a folder to write results to
                [parentFolder,filename,~] = fileparts(filepath);
                outputdir = fullfile(parentFolder,filename);
                if ~exist(outputdir,'dir'); mkdir(outputdir); end
                
                % write localisation file
                writetable(locsTable,fullfile(outputdir,'groundTruthLocalisations.txt'),'Delimiter','\t','WriteRowNames',true);

                % set up tif for writing tif stack
                t = Tiff(fullfile(outputdir,strcat(filename,'.tif')),'w'); % set up tif for writing
                tagstruct.ImageLength = nx; % image height
                tagstruct.ImageWidth  = ny; % image width
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Software = 'MATLAB';
                if strcmp(obj.Camera.cameraType,'nonoise')
                    tagstruct.Compression = Tiff.Compression.None;
                    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
                    tagstruct.Photometric = Tiff.Photometric.LinearRaw;
                    tagstruct.BitsPerSample = 32;
                    tagstruct.SamplesPerPixel = 1;
                else
                    tagstruct.Photometric = Tiff.Photometric.MinIsBlack; % https://de.mathworks.com/help/matlab/ref/tiff.html
                    tagstruct.BitsPerSample = 16;
                end
                
                simulations = 0;
            end
            
            F     = locsTable.frame;
            X     = locsTable.x;
            Y     = locsTable.y;
            Z     = locsTable.z;
            Phi   = locsTable.phi;
            Theta = locsTable.theta;
            Gamma = locsTable.gamma;
            I     = locsTable.photons;
            Delta = locsTable.delta;
            
            switch obj.typeNormalisation
                case 'relativeToStandardPSF'
                    Delta_normalisation = 0; % THIS SHOULD MAYBE BE A LOOKUP TABLE TO BE MORE CORRECT (separate normalisation for each Delta)
                    normalisation = obj.calculateNormalisation(Delta_normalisation);
                case 'individual'
                    % do nothing yet, the normalisation will be calculated as
                    % sum(img(:)) in the loop
                otherwise
                    disp("The entered value for the property typeNormalisation of class Simulation is unvalid.");
                    disp("It should be 'individual' or 'relativeToStandardPSF'.");
                    return
            end
            
            numFrames = max(F);
            for id_frame = 1:numFrames
                if verbose; fprintf('simulating frame %d/%d\n',id_frame,numFrames); end

                if logical(sum(ismember(id_frame,F))) % if there are localisations in this frame

                    % get list of localisations that occur in this frame
                    keep = logical(F == id_frame);
                    X_i     = X(keep);
                    Y_i     = Y(keep);
                    Z_i     = Z(keep);
                    Phi_i   = Phi(keep);
                    Theta_i = Theta(keep);
                    Gamma_i = Gamma(keep);
                    I_i     = I(keep);
                    Delta_i = Delta(keep);
                    
                    [mu_x_i,mu_y_i,mu_z_i] = obj.anglesToDipoleComponents(Phi_i,Theta_i);
                    
                    for id_emitter = 1:numel(X_i)
                        
                        % Generate basis functions
                        b = obj.getBasisFunctions_img(X_i(id_emitter),Y_i(id_emitter),Z_i(id_emitter),Delta_i(id_emitter));
                                
                        % Get electric field in pupil due to an emitter in point (x,y,z)
                        I_img_i = obj.getIntensityFromBasisFunctions_rotationInCone(mu_x_i(id_emitter),mu_y_i(id_emitter),mu_z_i(id_emitter),Gamma_i(id_emitter),b);
                        switch obj.typeNormalisation
                            case 'relativeToStandardPSF'
                                I_img_i = I_i(id_emitter)*I_img_i/normalisation;
                            case 'individual'
                                normalisation = sum(I_img_i(:));
                                I_img_i = I_i(id_emitter)*I_img_i/normalisation;
                        end
                        
                        % fprintf('%d\n',round(sum(I_img_i(:)))) % photon number (before bkgnd and camera model)
                        
                        if id_emitter == 1
                            I_img = I_img_i;
                        else
                            I_img = I_img + I_img_i;
                        end
                    end
                    
                    % add background
                    I_img = I_img + bkgnd;

                else % if there are no localisations in this frame
                    I_img = zeros(obj.fov,obj.fov) + bkgnd;
                end
                
                % simulate detection noise
                I_img = obj.Camera.cameraNoiseModel(I_img);
                
                if flag_write
                    % write away frame
                    if strcmp(obj.Camera.cameraType,'nonoise')
                        I_img = single(I_img);
                    else
                        I_img = uint16(I_img);
                    end
                    setTag(t,tagstruct);
                    write(t,I_img);
                    writeDirectory(t)
                else
                    simulations(:,:,id_frame) = I_img;
                end
                %imshow(kron(I_img,ones(2)),[]); title(sprintf('frame %d',id_frame)); pause(0.001)

            end
            if flag_write
                close(t); if verbose; fprintf('done!\nData saved in %s\n',outputdir); end
            else
                if verbose; fprintf('done!\n'); end
            end
        
        end

        function G = getGreensTensor_bfp(obj,Delta)
            %GETGREENSTENSOR Get the 9 components of the Green's tensor
            if abs(obj.nMedium - obj.nImmersion) < 0.0001 % if there is no refractive index mismatch
                G = getGreensTensor_bfp_RImatched(obj);
            else % if there is a refractive index mismatch
                G = getGreensTensor_bfp_planarRIinterface(obj,Delta);
            end
        end
        
        function G = getGreensTensor_bfp_RImatched(obj)
            %% GETGREENSTENSOR_BFP_RIMATCHED Generate Green's tensor for dipole
            % in an index matched environment.
            %
            % Green's tensor for a dipole in an index-matched environment (e.g.
            % immersed in oil, above a glass cover slip imaged with oil-immersion
            % objective). The Green's tensor that is returned is already multiplied
            % with a rotation matrix that mimicks an objective that satisfies the Abbe
            % sine condition. Multiplication with dipole components will give the field
            % in the back focal plane of the objective.
            %
            % Reference: Backer and Moerner, J Phys Chem B, 2014 (https://doi.org/10.1021/jp501778z)
            %
            % @param phi, rho..: array, polar coordinates of points in back focal plane
            %                    for which to calculate the rotated Green's tensor
            % @param wavelength: float, wavelength of the emitted light (in si units)
            % @param n.........: float, refractive index of the medium in which the dipole is immersed
            % @param f_obj.....: float, focal length of the objective lens (e.g. 3 mm, in si units)
            % @param f_obj.....: float, focal length of the objective lens (e.g. 3 mm, in si units)
            %
            % @output G........: struct, structure containing the Green's tensor components
            
            phi = obj.BFPphi;
            rho = obj.BFPrho;
            aperture = obj.BFPaperture;
            
            % Get the 9 components of the Green's tensor (last three are zero)
            nBFP = 1; % bfp is in air
            C = (exp(1j*obj.nMedium*obj.k0*obj.fObj)/(4*pi*obj.fObj))*sqrt(obj.nMedium./(nBFP*sqrt(1 - rho.^2)));
            
            G = struct;
            G.xx = aperture.*C.*(sin(phi).^2 + (cos(phi).^2).*sqrt(1 - rho.^2));
            G.xy = aperture.*C.*(sin(2*phi).*(sqrt(1 - rho.^2) - 1)/2);
            G.xz = aperture.*C.*(-rho.*cos(phi));
            G.yx = aperture.*C.*(sin(2*phi).*(sqrt(1 - rho.^2) - 1)/2);
            G.yy = aperture.*C.*(cos(phi).^2 + (sin(phi).^2).*sqrt(1 - rho.^2));
            G.yz = aperture.*C.*(-rho.*sin(phi));
            G.zx = zeros(size(aperture));
            G.zy = zeros(size(aperture));
            G.zz = zeros(size(aperture));
        end
        
        function G = getGreensTensor_bfp_planarRIinterface(obj,Delta)
            %% GETGREENSTENSOR_BFP_PLANARRIINTERFACE Generate rotated Green's tensor
            % for dipole in index matched environment.
            %
            % Modified Green's tensor for a dipole near a planar refractive index
            % interface (e.g. immersed in water, above a glass cover slip). The Green's
            % tensor that is returned is already multiplied with a rotation matrix that
            % mimicks an objective that satisfies the Abbe sine condition.
            % Multiplication with dipole components will give the field in the back
            % focal plane of the objective.
            %
            % Reference: Backer and Moerner, J Phys Chem B, 2014 (https://doi.org/10.1021/jp501778z)
            %
            % @param phi, rho..: array, polar coordinates of points in back focal plane
            %                    for which to calculate the rotated Green's tensor
            % @param wavelength: float, wavelength of the emitted light (in si units)
            % @param n.........: float, refractive index of the medium in which the dipole is immersed
            % @param f_obj.....: float, focal length of the objective lens (e.g. 3 mm, in si units)
            % @param Delta.....: float, distance of dipole from interface (in si units)
            %
            % @output G........: struct, structure containing the Green's tensor components

            phi = obj.BFPphi;
            rho = obj.BFPrho;
            aperture = obj.BFPaperture;
            n1 = obj.nImmersion;
            n2 = obj.nMedium;
            
            % Polar inclination of rays traveling from medium 1 to medium 2
            theta1 = asin(rho);
            theta2 = asin(rho*n1/n2);

            % Fresnel transmission coefficients for P- and S-polarized light
            [t_s, t_p] = obj.getFresnelTransmissionCoefficients(theta2,theta1,n1,n2);

            % Get some other coefficients
            exp_term = exp(1j*obj.k0*Delta*n2*sqrt(1 - (n1/n2)*rho.^2));
            c1 = ((n1/n2)^2)*(cos(theta1)./cos(theta2)).*t_p.*exp_term;
            c2 = (n1/n2).*t_p.*exp_term;
            c3 = (n1/n2)*(cos(theta1)./cos(theta2)).*t_s.*exp_term;

            % Get the 9 components of the Green's tensor (last three are zero)
            nBFP = 1; % bfp is in air
            C = (exp(1j*n1*obj.k0*obj.fObj)/(4*pi*obj.fObj))*sqrt(n1./(nBFP*sqrt(1 - rho.^2)));

            G = struct;
            G.xx = aperture.*C.*(c3.*sin(phi).^2 + c2.*(cos(phi).^2).*sqrt(1 - rho.^2));
            G.xy = aperture.*C.*(-sin(2*phi).*(c3 - c2.*sqrt(1 - rho.^2))/2);
            G.xz = aperture.*C.*(-c1.*rho.*cos(phi));
            G.yx = aperture.*C.*(-sin(2*phi).*(c3 - c2.*sqrt(1 - rho.^2))/2);
            G.yy = aperture.*C.*(c2.*sqrt(1 - rho.^2) + c3.*cos(phi).^2 - c2.*(cos(phi).^2).*sqrt(1 - rho.^2));
            G.yz = aperture.*C.*(-c1.*rho.*sin(phi));
            G.zx = zeros(size(aperture));
            G.zy = zeros(size(aperture));
            G.zz = zeros(size(aperture));
        end
        
        function g = getGreensTensor_img(obj,x,y,z,Delta)
            %GETGREENSTENSOR_IMAGE Get the Fourier transforms of the 9
            %components of the Green's tensor.
            G = obj.getGreensTensor_bfp(Delta);
            G = obj.applyTipTiltDefocusToGreensTensorElements(G,x,y,z);
            g = obj.applyImageFormationToGreensTensorElements(G);
        end
        
        function B = getBasisFunctions_bfp(obj,Delta)
            %GENERATEBFPBASISFUNCTIONS
            % Get the 9 components of the Green's tensor
            G = obj.getGreensTensor_bfp(Delta);
            % get the back focal plane basis functions
            B = obj.getBasisFunctionsFromGreensTensor(G);
        end
        
        function b = getBasisFunctions_img(obj,x,y,z,Delta)
            %GETBASISFUNCTIONS_IMG
            g = obj.getGreensTensor_img(x,y,z,Delta);
            b = obj.getBasisFunctionsFromGreensTensor(g);
        end
                
        function G = applyTipTiltDefocusToGreensTensorElements(obj,G,x,y,z)
            %APPLYTIPTILTDEFOCUSTOGREENSTENSORELEMENTS Apply tip, tilt and
            %defocus to put the Green's tensor elements for a dipole in
            %(x,y,z) in object space.
            phase_x = obj.nMedium*obj.k0*x*obj.BFPrho.*cos(obj.BFPphi);
            phase_y = obj.nMedium*obj.k0*y*obj.BFPrho.*sin(obj.BFPphi);
            phase_z = obj.nMedium*obj.k0*z*sqrt(1 - obj.BFPrho.^2);
            total_phase = phase_x + phase_y + phase_z;
            G.xx = G.xx.*exp(1j*total_phase);
            G.xy = G.xy.*exp(1j*total_phase);
            G.xz = G.xz.*exp(1j*total_phase);
            G.yx = G.yx.*exp(1j*total_phase);
            G.yy = G.yy.*exp(1j*total_phase);
            G.yz = G.yz.*exp(1j*total_phase);
        end
            
        function g = applyImageFormationToGreensTensorElements(obj,G)
            %APPLYIMAGEFORMATIONTOGREENSTENSORELEMENTS
            g = struct;
            g.xx = obj.cropToROICentered(obj.imageFormationParaxial(G.xx),obj.fov);
            g.xy = obj.cropToROICentered(obj.imageFormationParaxial(G.xy),obj.fov);
            g.xz = obj.cropToROICentered(obj.imageFormationParaxial(G.xz),obj.fov);
            g.yx = obj.cropToROICentered(obj.imageFormationParaxial(G.yx),obj.fov);
            g.yy = obj.cropToROICentered(obj.imageFormationParaxial(G.yy),obj.fov);
            g.yz = obj.cropToROICentered(obj.imageFormationParaxial(G.yz),obj.fov);
        end
        
        function [Ex,Ey] = getElectricFieldFromGreensTensor(obj,mu_x,mu_y,mu_z,greensTensor)
            %GETELECTRICFIELDFROMGREENSTENSOR Calculate the electric field due to
            %an emitter with an emission dipole moment mu, from the Green's tensor.
            [~,mu_x,mu_y,mu_z] = obj.checkIsUnitVector(mu_x,mu_y,mu_z); % check that norm of dipole vector is 1, otherwise rescale
            Ex = mu_x*greensTensor.xx + mu_y*greensTensor.xy + mu_z*greensTensor.xz; % x-component electric field in bfp
            Ey = mu_x*greensTensor.yx + mu_y*greensTensor.yy + mu_z*greensTensor.yz; % y-component electric field in bfp
        end
                
        function I = getIntensityFromBasisFunctions_rotationInCone(obj,mu_x,mu_y,mu_z,gamma,basisFunctions)
            %GETINTENSITYFROMBASISFUNCTIONS_ROTATIONINCONE
            
            % check that norm of dipole vector is 1, otherwise rescale
            [~,mu_x,mu_y,mu_z] = obj.checkIsUnitVector(mu_x,mu_y,mu_z);
            
            % get image basis function coefficients assuming rotation in a cone
            mxx = mu_x*mu_x*gamma + (1-gamma)/3;
            myy = mu_y*mu_y*gamma + (1-gamma)/3;
            mzz = mu_z*mu_z*gamma + (1-gamma)/3; 
            mxy = mu_x*mu_y*gamma;
            mxz = mu_x*mu_z*gamma;
            myz = mu_y*mu_z*gamma;

            % get image rotating molecule
            I = (basisFunctions.XXx + basisFunctions.XXy)*mxx + ...
                (basisFunctions.YYx + basisFunctions.YYy)*myy + ...
                (basisFunctions.ZZx + basisFunctions.ZZy)*mzz + ...
                (basisFunctions.XYx + basisFunctions.XYy)*mxy + ...
                (basisFunctions.XZx + basisFunctions.XZy)*mxz + ...
                (basisFunctions.YZx + basisFunctions.YZy)*myz;
        end
        
        function pad = calculatePadding(obj)
            %CALCULATEPADDING
            
            % calculate desired virtual pixel size
            virtualPixSize = obj.Camera.camPixSize/obj.Microscope.magnificationObjective;
            
            % calculate the number of 'simulation pixels' in the bfp
            bfpPixSize = obj.bfpDiameterNative/obj.numPixBFP; % pixel size back focal plane (simulation)

            totalPadding = (obj.wavelength*obj.fObj/virtualPixSize) - obj.bfpDiameterNative; % padding needed to get desired pixel size
            padding = totalPadding/2; % will pad on each side of bfp
            
            % convert to pixels
            pad = round(padding/bfpPixSize); % padding
        end
        
        function E = applyTipTiltDefocus(obj,E,x,y,z)
            % APPLYTIPTILTDEFOCUS Apply tip, tilt and defocus to put
            % dipole in (x,y,z) in object space.
            phase_x = obj.nMedium*obj.k0*x*obj.BFPrho.*cos(obj.BFPphi);
            phase_y = obj.nMedium*obj.k0*y*obj.BFPrho.*sin(obj.BFPphi);
            phase_z = obj.nMedium*obj.k0*z*sqrt(1 - obj.BFPrho.^2);
            total_phase = phase_x + phase_y + phase_z;
            E = E.*exp(1j*total_phase);
        end
        
        function E = imageFormationParaxial(obj,E)
            %IMAGEFORMATIONPARAXIAL Simulate image formation from a bfp
            %electric field using paraxial approximation (i.e. fourier
            %transform after padding).
            E = fftshift(fft2(fftshift(padarray(E,[obj.pad,obj.pad]))));
        end
        
        function normalisation = calculateNormalisation(obj,Delta)
            %CALCULATENORMALISATION
            G = obj.getGreensTensor_bfp(Delta);
            g = obj.applyImageFormationToGreensTensorElements(G);
            b = obj.getBasisFunctionsFromGreensTensor(g);
            if obj.nMedium >= obj.NA % use in in-plane oriented emitter for normalisation
                I = obj.getIntensityFromBasisFunctions_rotationInCone(1,0,0,1,b);
            else % use an out-of-plane oriented emitter for normalisation
                I = obj.getIntensityFromBasisFunctions_rotationInCone(0,0,1,1,b);
            end
            normalisation = sum(I(:));
        end
        
        function [fig1,fig2] = showGreensTensorComponents_bfp(obj)
            %SHOWGREENSTENSORCOMPONENTS_BFP
            normalisation = max(max(obj.getAmplitude(obj.GreensTensorZ0_bfp.xx)));
            fig1 = figure('Name',"Green's tensor components at back focal plane (amplitude)",'position',[100 100 1000 450]);
            subplot(2,3,1); imshow(obj.getAmplitude(obj.GreensTensorZ0_bfp.xx)/normalisation,[]); c = colorbar; title('abs(G_{xx})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,2); imshow(obj.getAmplitude(obj.GreensTensorZ0_bfp.xy)/normalisation,[]); c = colorbar; title('abs(G_{xy})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,3); imshow(obj.getAmplitude(obj.GreensTensorZ0_bfp.xz)/normalisation,[]); c = colorbar; title('abs(G_{xz})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,4); imshow(obj.getAmplitude(obj.GreensTensorZ0_bfp.yx)/normalisation,[]); c = colorbar; title('abs(G_{yx})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,5); imshow(obj.getAmplitude(obj.GreensTensorZ0_bfp.yy)/normalisation,[]); c = colorbar; title('abs(G_{yy})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,6); imshow(obj.getAmplitude(obj.GreensTensorZ0_bfp.yz)/normalisation,[]); c = colorbar; title('abs(G_{yz})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            
            fig2 = figure('Name',"Green's tensor components at back focal plane (phase)",'position',[800 100 1000 450]);
            subplot(2,3,1); imshow(obj.getPhase(obj.GreensTensorZ0_bfp.xx),[]); c = colorbar; colormap hsv; title('angle(G_{xx})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,2); imshow(obj.getPhase(obj.GreensTensorZ0_bfp.xy),[]); c = colorbar; colormap hsv; title('angle(G_{xy})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,3); imshow(obj.getPhase(obj.GreensTensorZ0_bfp.xz),[]); c = colorbar; colormap hsv; title('angle(G_{xz})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,4); imshow(obj.getPhase(obj.GreensTensorZ0_bfp.yx),[]); c = colorbar; colormap hsv; title('angle(G_{yx})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,5); imshow(obj.getPhase(obj.GreensTensorZ0_bfp.yy),[]); c = colorbar; colormap hsv; title('angle(G_{yy})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,6); imshow(obj.getPhase(obj.GreensTensorZ0_bfp.yz),[]); c = colorbar; colormap hsv; title('angle(G_{yz})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
        end
        
        function [fig1,fig2] = showGreensTensorComponents_img(obj)
            %SHOWGREENSTENSORCOMPONENTS_IMAGE
            normalisation = max(max(obj.getAmplitude(obj.GreensTensorZ0_img.xx)));
            fig1 = figure('Name',"Green's tensor components at image plane (amplitude)",'position',[100 100 1000 450]);
            subplot(2,3,1); imshow(obj.getAmplitude(obj.GreensTensorZ0_img.xx)/normalisation,[]); c = colorbar; title('abs(g_{xx})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,2); imshow(obj.getAmplitude(obj.GreensTensorZ0_img.xy)/normalisation,[]); c = colorbar; title('abs(g_{xy})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,3); imshow(obj.getAmplitude(obj.GreensTensorZ0_img.xz)/normalisation,[]); c = colorbar; title('abs(g_{xz})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,4); imshow(obj.getAmplitude(obj.GreensTensorZ0_img.yx)/normalisation,[]); c = colorbar; title('abs(g_{yx})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,5); imshow(obj.getAmplitude(obj.GreensTensorZ0_img.yy)/normalisation,[]); c = colorbar; title('abs(g_{yy})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(2,3,6); imshow(obj.getAmplitude(obj.GreensTensorZ0_img.yz)/normalisation,[]); c = colorbar; title('abs(g_{yz})','FontWeight','normal'); caxis([0 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            
            fig2 = figure('Name',"Green's tensor components at image plane (phase)",'position',[800 100 1000 450]);
            subplot(2,3,1); imshow(obj.getPhase(obj.GreensTensorZ0_img.xx),[]); c = colorbar; colormap hsv; title('angle(g_{xx})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,2); imshow(obj.getPhase(obj.GreensTensorZ0_img.xy),[]); c = colorbar; colormap hsv; title('angle(g_{xy})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,3); imshow(obj.getPhase(obj.GreensTensorZ0_img.xz),[]); c = colorbar; colormap hsv; title('angle(g_{xz})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,4); imshow(obj.getPhase(obj.GreensTensorZ0_img.yx),[]); c = colorbar; colormap hsv; title('angle(g_{yx})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,5); imshow(obj.getPhase(obj.GreensTensorZ0_img.yy),[]); c = colorbar; colormap hsv; title('angle(g_{yy})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
            subplot(2,3,6); imshow(obj.getPhase(obj.GreensTensorZ0_img.yz),[]); c = colorbar; colormap hsv; title('angle(g_{yz})','FontWeight','normal'); caxis([-pi,pi]); c.Ticks = -pi:pi:pi; c.TickLabels = {'-\pi','0','\pi'}; set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize*1.2;
        end
        
        function fig = showBasisFunctions_bfp(obj)
            %SHOWBASISFUNCTIONS_BFP
            normalisation = max(obj.BasisFunctionsZ0_bfp.XXx(:));
            fig = figure('Name','Back focal plane basis functions','position',[50 50 1000 500]);
            subplot(3,4,1); imshow(obj.BasisFunctionsZ0_bfp.XXx/normalisation,[]); c = colorbar; title('XX^x','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,2); imshow(obj.BasisFunctionsZ0_bfp.XYx/normalisation,[]); c = colorbar; title('XY^x','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,3); imshow(obj.BasisFunctionsZ0_bfp.XXy/normalisation,[]); c = colorbar; title('XX^y','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,4); imshow(obj.BasisFunctionsZ0_bfp.XYy/normalisation,[]); c = colorbar; title('XY^y','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;

            subplot(3,4,5); imshow(obj.BasisFunctionsZ0_bfp.YYx/normalisation,[]); c = colorbar; title('YY^x','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,6); imshow(obj.BasisFunctionsZ0_bfp.XZx/normalisation,[]); c = colorbar; title('XZ^x','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,7); imshow(obj.BasisFunctionsZ0_bfp.YYy/normalisation,[]); c = colorbar; title('YY^y','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,8); imshow(obj.BasisFunctionsZ0_bfp.XZy/normalisation,[]); c = colorbar; title('XZ^y','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;

            subplot(3,4,9);  imshow(obj.BasisFunctionsZ0_bfp.ZZx/normalisation,[]); c = colorbar; title('ZZ^x','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,10); imshow(obj.BasisFunctionsZ0_bfp.YZx/normalisation,[]); c = colorbar; title('YZ^x','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,11); imshow(obj.BasisFunctionsZ0_bfp.ZZy/normalisation,[]); c = colorbar; title('ZZ^y','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,12); imshow(obj.BasisFunctionsZ0_bfp.YZy/normalisation,[]); c = colorbar; title('YZ^y','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
        end
        
        function fig = showBasisFunctions_img(obj)
            %SHOWBASISFUNCTIONS_IMAGE
            normalisation = max(obj.BasisFunctionsZ0_img.XXx(:));
            fig = figure('Name','Image plane basis functions','position',[50 50 1000 500]);
            subplot(3,4,1); imshow(obj.BasisFunctionsZ0_img.XXx/normalisation,[]); c = colorbar; title('XX^x','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,2); imshow(obj.BasisFunctionsZ0_img.XYx/normalisation,[]); c = colorbar; title('XY^x','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,3); imshow(obj.BasisFunctionsZ0_img.XXy/normalisation,[]); c = colorbar; title('XX^y','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,4); imshow(obj.BasisFunctionsZ0_img.XYy/normalisation,[]); c = colorbar; title('XY^y','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            
            subplot(3,4,5); imshow(obj.BasisFunctionsZ0_img.YYx/normalisation,[]); c = colorbar; title('YY^x','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,6); imshow(obj.BasisFunctionsZ0_img.XZx/normalisation,[]); c = colorbar; title('XZ^x','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,7); imshow(obj.BasisFunctionsZ0_img.YYy/normalisation,[]); c = colorbar; title('YY^y','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,8); imshow(obj.BasisFunctionsZ0_img.XZy/normalisation,[]); c = colorbar; title('XZ^y','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;

            subplot(3,4,9);  imshow(obj.BasisFunctionsZ0_img.ZZx/normalisation,[]); c = colorbar; title('ZZ^x','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,10); imshow(obj.BasisFunctionsZ0_img.YZx/normalisation,[]); c = colorbar; title('YZ^x','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,11); imshow(obj.BasisFunctionsZ0_img.ZZy/normalisation,[]); c = colorbar; title('ZZ^y','FontWeight','normal'); caxis([0  1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
            subplot(3,4,12); imshow(obj.BasisFunctionsZ0_img.YZy/normalisation,[]); c = colorbar; title('YZ^y','FontWeight','normal'); caxis([-1 1]); set(gca,'Fontsize',obj.fontsize); c.FontSize = obj.fontsize;
        end
        
    end    
    
    methods (Static)
       
        function flux = getPhotonFlux(QY,molarExtCoeff,wavelength,laserPower)
            % GETPHOTONFLUX Calculate the flux of photon emission (emitted
            % photons/s) of a dye given descriptive parameters of the
            % molecule and laser power. Based on
            % https://www.nature.com/articles/s41592-019-0364-4#Sec19
            % 
            % QY...........: quantum yield of the molecule
            % molarExtCoeff: molar extinction coefficient or absorptivity (cm^2/molecule)
            % wavelength...: wavelength of the laser (m)
            % laserPower...: laser power density (W/cm2)
    
            N_a = 6.0221409e+23; % Avogadro's constant
            h = 6.62607004e-34; % m^2*kg/s, Planck's constant
            c = 299792458; % m/s, speed of light
            
            a = 1000*log(10)*molarExtCoeff/N_a; % cm^2, absorption cross-sectoin
            e = h*c/wavelength; % energy one photon
            
            flux = QY*laserPower*a/e;
        
        end

        function [t,pulse1,pulse2] = getLaserPulseTrace(t_angle,angle,w1,w2,flagshow)
            % GETLASERPULSETRAIN

            threshold1 = sin((pi*(1-2*w1))/2);
            threshold2 = -sin((pi*(1-2*w2))/2);
            pulse1 = zeros(size(angle));
            pulse2 = zeros(size(angle));
            pulse1(angle >= threshold1) = 1;
            pulse2(angle <= threshold2) = 1;
            t = t_angle;
            
            if flagshow
                figure;
                plot(t,angle); hold on
                plot(t,pulse1);
                plot(t,-pulse2);
                xlim([0,inf]); ylim([-1.1,1.1]);
                xlabel('Time (ms)'); ylabel('Voltage (a.u.)')
                legend('Output RS','Input laser 1','Input laser 2')
                set(gca,'FontSize',10)
                
                faceAlpha = 0.7;
                figure;
                area(t,10*pulse1,'FaceColor',[0.5804,0.4039,0.7412],'FaceAlpha',faceAlpha,'EdgeColor','none'); hold on
                area(t,-10*pulse1,'FaceColor',[0.5804,0.4039,0.7412],'FaceAlpha',faceAlpha,'EdgeColor','none');
                area(t,10*pulse2,'FaceColor',[0.1725,0.6275,0.1725],'FaceAlpha',faceAlpha,'EdgeColor','none');
                area(t,-10*pulse2,'FaceColor',[0.1725,0.6275,0.1725],'FaceAlpha',faceAlpha,'EdgeColor','none');
                plot(t,angle,'k','LineWidth',1.2);
                xlim([0,inf]); ylim([-1.1,1.1]);
                xlabel('Time (ms)'); ylabel('Position');
                set(gca,'FontSize',10)
            end
        end

        function [t_s,t_p] = getFresnelTransmissionCoefficients(theta_i,theta_t,n1,n2)
            % GETFRESNELTRANSMISSIONCOEFFICIENTS Calculate the Fresnel
            % transmission coefficients. These describe the transmission of
            % light when incident on an interface between different optical media.
            % Reference: M. Born and E. Wolf. Principles of Optics, Ch 1
            % Cambridge University Press, 7th edition, 2019
            %   theta: float or array, angle(s) at which coefficients are calculated
            %   n1: float, refractive index of the first medium (e.g. of glass, immersion oil)
            %   n2: float, refractive index of the second medium (medium of sample, e.g. water or air)
            %   t_s: float or array (depending on type input parameter theta),
            %        Fresnel transmission coefficients for s-polarized light
            %   t_p: float or array (depending on type input parameter theta)
            %        Fresnel transmission coefficients for p-polarized light
            % 
            % returns: t_s: fresnel transmission coefficient for s-polarized component
            %               of the electric field (normal to the plane of incidence)
            %          t_p: fresnel transmission coefficient for p-polarized component
            %               of the electric field (parallel to the plane of incidence)

            t_s_num = 2*n1*cos(theta_i);
            t_s_den = n1*cos(theta_i) + n2*cos(theta_t);
            t_s = t_s_num./t_s_den;

            t_p_num = 2*n1*cos(theta_i);
            t_p_den = n2*cos(theta_i) + n1*cos(theta_t);
            t_p = t_p_num./t_p_den;
        end
        
        function basisFunctions = getBasisFunctionsFromGreensTensor(g)
            %GETBASISFUNCTIONSFROMGREENSTENSOR
            basisFunctions.XXx = abs(g.xx).^2;
            basisFunctions.XXy = abs(g.yx).^2;
            basisFunctions.YYx = abs(g.xy).^2;
            basisFunctions.YYy = abs(g.yy).^2;
            basisFunctions.ZZx = abs(g.xz).^2;
            basisFunctions.ZZy = abs(g.yz).^2;
            basisFunctions.XYx = 2*real(g.xx.*conj(g.xy));
            basisFunctions.XYy = 2*real(g.yx.*conj(g.yy));
            basisFunctions.XZx = 2*real(g.xx.*conj(g.xz));
            basisFunctions.XZy = 2*real(g.yx.*conj(g.yz));
            basisFunctions.YZx = 2*real(g.xy.*conj(g.xz));
            basisFunctions.YZy = 2*real(g.yy.*conj(g.yz));
        end
        
        function [mu_x,mu_y,mu_z] = anglesToDipoleComponents(phi,theta)
            %ANGLESTODIPOLECOMPONENTS Convert angles (phi,theta) to dipole
            %components mu_x,mu_y,mu_z.
            mu_x = sin(theta).*cos(phi);
            mu_y = sin(theta).*sin(phi);
            mu_z = cos(theta);
        end
        
        function [unitVector,mu_x,mu_y,mu_z] = checkIsUnitVector(mu_x,mu_y,mu_z)
            %CHECKISUNITVECTOR Check whether dipole vector is a unit
            %vector. If not, renormalise such that the norm is 1.
            norm = sqrt(mu_x^2 + mu_y^2 + mu_z^2);
            if (norm - 1) < 0.00001
                unitVector = true;
            else
                fprintf('Norm of dipole components is %.4f. Components renormalised such that norm is 1.\n',norm);
                unitVector = false;
                mu_x = mu_x/norm;
                mu_y = mu_y/norm;
                mu_z = mu_z/norm;
            end
        end
        
        function crop = cropToROICentered(img,fov)
            %CROPTOROICENTERED
            [nx,ny] = size(img);
            xc = ceil(nx/2) + 1;
            yc = ceil(ny/2) + 1;
            w = floor(fov/2);
            crop = img(xc-w:xc-w+fov-1,...
                       yc-w:yc-w+fov-1,:);
        end
        
        function amplitude = getAmplitude(E)
            %GETAMPLITUDE Get amplitude of electric field.
            amplitude = abs(E);
        end
        
        function phase = getPhase(E)
            %GETPHASE Get phase of electric field.
            phase = angle(E);
        end
        
        function intensity = getIntensity(E)
            %GETINTENSITY Get intensity of electric field.
            intensity = abs(E).^2;
        end
        
    end
    
end
