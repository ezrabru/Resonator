classdef Characterisation
    %CHARACTERISATION object contains properties and methods to
    %characterise performance of an algorithm, setup and conditions.
    
    properties
        Scope % instance of the Microscope class
        Cam % instance of the Camera class
        Sim % instance of the Simulation class
    end
    
    methods
        
        function obj = Characterisation(Microscope,Camera,Simulation)
            %CHARACTERISATION Construct an instance of this class
            obj.Scope = Microscope;
            obj.Cam   = Camera;
            obj.Sim   = Simulation;
        end
        
        function crlbPhi = getCRLB_phi(obj,thetaRange,phiRange,gamma,photons,bkgnd)
            %GETCRLB_PHI
            
            % check if phiRange is an equally spaced array
            diffRange = diff(phiRange);
            dPhi = diffRange(1);
            equallySpaced = all(diffRange - dPhi < 1e-6);
            
            if equallySpaced
                crlb = nan(numel(thetaRange),numel(phiRange)-1);
                for i=1:numel(thetaRange)
                    fprintf('Calculating CRLB phi for theta %d/%d\n',i,numel(thetaRange));
                    locsTable = obj.fixedThetaVaryPhi(phiRange,thetaRange(i),gamma,photons);
                    [simulations_phi,~] = obj.Sim.simulateSample(locsTable,0,bkgnd,0);
                    crlb(i,:) = obj.getCRLB(simulations_phi,dPhi); % Calculate CRLB of Phi for fixed Theta
                end
                crlbPhi.crlb = crlb*180/pi;
                crlbPhi.phi = phiRange*180/pi;
                crlbPhi.theta = thetaRange*180/pi;
            else
                disp('The points in phiRange need to be equally spaced.'); return
            end
        end
        
        function crlbTheta = getCRLB_theta(obj,thetaRange,phiRange,gamma,photons,bkgnd)
            %GETCRLB_THETA
            
            % check if phiRange is an equally spaced array
            diffRange = diff(phiRange);
            dTheta = diffRange(1);
            equallySpaced = all(diffRange - dTheta < 1e-6);
            
            if equallySpaced
                crlb = nan(numel(thetaRange)-1,numel(phiRange));
                for i=1:numel(phiRange)
                   fprintf('Calculating CRLB theta for phi %d/%d\n',i,numel(phiRange));
                    locsTable = obj.fixedPhiVaryTheta(phiRange(i),thetaRange,gamma,photons);
                    [simulations_theta,~] = obj.Sim.simulateSample(locsTable,0,bkgnd,0);
                    crlb(:,i) = obj.getCRLB(simulations_theta,dTheta); % Calculate CRLB of Theta for fixed Phi
                end
                crlbTheta.crlb   = crlb*180/pi;
                crlbTheta.phi    = phiRange*180/pi;
                crlbTheta.theta  = thetaRange*180/pi;
            else
                disp('The points in thetaRange need to be equally spaced.'); return
            end
        end
        
        function [estimates,simulationsNoise,simulations,allSimulations] = ...
                calculateBiasPrecision(obj,phiRange,thetaRange,variedParam,paramRange,repeats,...
                fixedParams,bkgnd,fov,conversionmethod,localisationParam,w,returnAllSimulations)
            %CALCULATEBIASPRECISION
            % Note: the molecule will be in focus at z = 0, Delta = 0
            
            if returnAllSimulations
                allSimulations = nan(fov,fov,numel(phiRange)*numel(thetaRange)*numel(paramRange)*repeats);
            else
                allSimulations = [];
            end
            
            % initialise
            estimates = {};
            
            estimates.variedParam = variedParam;
            estimates.variedParamRange = paramRange;
            estimates.phiRange = phiRange;
            estimates.thetaRange = thetaRange;
            estimates.repeats = repeats;
            estimates.fixedParams = fixedParams;
            estimates.bkgnd = bkgnd;
            estimates.explanationDimensions.dim1 = 'phi';
            estimates.explanationDimensions.dim2 = 'theta';
            estimates.explanationDimensions.dim3 = variedParam;
            estimates.explanationDimensions.dim4 = 'repeats';
            
            estimates.phi     = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            estimates.theta   = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            estimates.avgDoLP = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            estimates.netDoLP = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            estimates.x       = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            estimates.y       = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            % estimates.photons = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            
            count = 0;
            for id_param = 1:numel(paramRange)
                
                switch variedParam
                    case 'photons'
                        photons = paramRange(id_param);
                        gamma   = fixedParams.gamma;
                    case 'gamma'
                        photons = fixedParams.photons;
                        gamma   = paramRange(id_param);
                    otherwise
                        disp("Invalid value for input parameter 'variedParam' in method 'calculateBiasPrecision' of class 'Characterisation'.");
                        disp("Supported values are 'photons' or 'gamma'."); return
                end
                
                for id_theta = 1:numel(thetaRange)
                    
                    % simulate images without noise or background, for fixed theta
                    Cam_nonoise = Camera('nonoise',obj.Cam.camPixSize,0,obj.Cam.bit,obj.Cam.polariserUnitMosaic);
                    locsTable = obj.fixedThetaVaryPhi(phiRange,thetaRange(id_theta),gamma,1);
                    Sim = Simulation(locsTable,obj.Scope,Cam_nonoise,fov,obj.typeNormalisation);
                    [simulations,~] = Sim.simulateSample(0,0,0);
                    
                    masks = Sim.masksFov;
                    Proc = Process(obj.Scope,obj.Cam,masks,conversionmethod,localisationParam,w);
                    
                    % normalise and add noise in a loop and analyse
                    simulationsNoise = nan(Sim.fov,Sim.fov,numel(phiRange),numel(thetaRange),numel(paramRange)); % save a few example images with noise
                    for id_phi = 1:numel(phiRange)
                        for id_rep = 1:repeats
                            count = count + 1;
                            if ~mod(count,100); fprintf('Processing frame %d/%d\n',count,numel(phiRange)*numel(thetaRange)*numel(paramRange)*repeats); end
                            
                            frame = abs(photons*simulations(:,:,id_phi) + bkgnd); % abs to avoid loosing whole image to 'nan' after adding poisson noise (poissrnd(negative number) gives 'nan')
                            frame = obj.Cam.cameraNoiseModel(frame);
                            if or(strcmp(obj.Cam.cameraType,'nonoise'),strcmp(obj.Cam.cameraType,'ideal'))
                                frame = frame - bkgnd;
                            elseif or(strcmp(obj.Cam.cameraType,'sCMOS'),strcmp(obj.Cam.cameraType,'EMCCD'))
                                frame = frame - bkgnd - obj.Cam.noiseParam.offset;
                            end
                            
                            if returnAllSimulations
                                allSimulations(:,:,count) = frame;
                            end
                            
                            if id_rep == 1
                                simulationsNoise(:,:,id_phi,id_theta,id_param) = frame;
                            end
                            
                            [x,y,phi,theta,avgDoLP,netDoLP,int] = Proc.analyseSimulation(frame);
                            
                            estimates.x(      id_phi,id_theta,id_param,id_rep) = x;
                            estimates.y(      id_phi,id_theta,id_param,id_rep) = y;
                            estimates.phi(    id_phi,id_theta,id_param,id_rep) = phi;
                            estimates.theta(  id_phi,id_theta,id_param,id_rep) = theta;
                            estimates.avgDoLP(id_phi,id_theta,id_param,id_rep) = avgDoLP;
                            estimates.netDoLP(id_phi,id_theta,id_param,id_rep) = netDoLP;
                            estimates.photons(id_phi,id_theta,id_param,id_rep) = int;
                        end
                    end
                end
            end
            fprintf('done!\n')
        end
        
        function [results,simulationsNoise,simulations] = ...
                calculateBiasPrecision_phi_vs_photons(obj,phiRange,photonRange,repeats,...
                Theta,gamma,bkgnd,fov,conversionmethod,w,flagshow)
            %CALCULATEBIASPRECISION_PHI_VS_PHOTONS
            % Note: the molecule will be in focus at z = 0, Delta = 0
            
            % simulate images without noise or background
            Cam_nonoise = Camera('nonoise',obj.Cam.camPixSize,0,obj.Cam.bit,obj.Cam.polariserUnitMosaic);
            locsTable = obj.fixedThetaVaryPhi(phiRange,Theta,gamma,1);
            Sim = Simulation(locsTable,obj.Scope,Cam_nonoise,fov,obj.typeNormalisation);
            [simulations,~] = Sim.simulateSample(0,0);
            
            masks = Sim.masksFov;
            Proc = Process(obj.Scope,obj.Cam,masks,conversionmethod,w);
            
            % normalise and add noise in a loop and analyse
            simulationsNoise = nan(Sim.fov,Sim.fov,numel(phiRange),numel(photonRange)); % save a few example images with noise
            phiEst = nan(numel(phiRange),numel(photonRange),repeats);
            count = 0;
            for id_phi=1:numel(phiRange)
                for id_photons = 1:numel(photonRange)
                    for id_rep=1:repeats
                        count = count + 1;
                        if ~mod(count,100); fprintf('Analysing frame %d/%d\n',count,numel(phiRange)*numel(photonRange)*repeats); end
                        
                        frame = abs(photonRange(id_photons)*simulations(:,:,id_phi) + bkgnd); % abs to avoid loosing whole image to 'nan' after adding poisson noise (poissrnd(negative number) gives 'nan')
                        frame = obj.Cam.cameraNoiseModel(frame);
                        if or(strcmp(obj.Cam.cameraType,'nonoise'),strcmp(obj.Cam.cameraType,'ideal'))
                            frame = frame - bkgnd;
                        elseif or(strcmp(obj.Cam.cameraType,'sCMOS'),strcmp(obj.Cam.cameraType,'EMCCD'))
                            frame = frame - bkgnd - obj.Cam.noiseParam.offset;
                        end
                        
                        if flagshow
                            if id_rep == 1
                                simulationsNoise(:,:,id_phi,id_photons) = frame;
                            end
                        end
                        
                        [phi,~,~,~] = Proc.analyseSimulation(frame);
                        phiEst(id_phi,id_photons,id_rep) = phi;
                    end
                end
            end
            fprintf('done!\n')
            
            % show example images with and without noise
            if flagshow
                figure('Name','Example simulations without noise');
                for i=1:floor(sqrt(numel(phiRange)))^2
                    subplot(floor(sqrt(numel(phiRange))),floor(sqrt(numel(phiRange))),i);
                    imshow(simulations(:,:,i),[])
                end
                for id_photons=1:numel(photonRange)
                    figure('Name',sprintf('Example simulations %d photons',photonRange(id_photons)));
                    for i=1:floor(sqrt(numel(phiRange)))^2
                        subplot(floor(sqrt(numel(phiRange))),floor(sqrt(numel(phiRange))),i);
                        img = simulationsNoise(:,:,i,id_photons);
                        imshow(img,[]); colorbar; title(sprintf('sum(img(:)) = %.1f',sum(img(:))))
                    end
                end
            end
            
            fprintf('Calculating bias and precision... ')
            % calculate bias
            truePhi = repmat(phiRange'*pi/180,1,numel(photonRange));
            phiErrorPhoton = nan(numel(phiRange),numel(photonRange),repeats);
            phiRemappedPhoton = nan(numel(phiRange),numel(photonRange),repeats);
            for id_reps=1:repeats
                [phiError, phiEstRemapped] = obj.remapPhiError(phiEst(:,:,id_reps),truePhi);
                phiErrorPhoton(:,:,id_reps) = phiError;
                phiRemappedPhoton(:,:,id_reps) = phiEstRemapped;
            end
            bias = mean(phiErrorPhoton,3);
            meanBias = mean(bias)*180/pi;
            stdBias = std(bias)*180/pi;
            minBias = min(bias)*180/pi;
            maxBias = max(bias)*180/pi;
            
            % calculate precision
            precision = std(phiErrorPhoton,[],3);
            meanPrecision = mean(precision)*180/pi;
            stdPrecision = std(precision)*180/pi;
            minPrecision = min(precision)*180/pi;
            maxPrecision = max(precision)*180/pi;
            
            % group results
            results.bias          = bias;
            results.meanError     = meanBias;
            results.stdError      = stdBias;
            results.minError      = minBias;
            results.maxError      = maxBias;
            results.precision     = precision;
            results.meanPrecision = meanPrecision;
            results.stdPrecision  = stdPrecision;
            results.minPrecision  = minPrecision;
            results.maxPrecision  = maxPrecision;
            fprintf('done!\n')
            
            %+++++++++++++++++++++
            % Plotting
            fontsize = 10;
            if flagshow
                fprintf('Generating plots... ')
                for id_photons = 1:numel(photonRange)
                    % check by plotting true vs estimate in scatter plot (rawest results)
                    figure('position',[100 450 1200 250],'Name',sprintf('%d photons',photonRange(id_photons)));
                    subplot(1,4,1)
                    for i=1:repeats
                        phiEst_rep_i = phiEst(:,id_photons,i);
                        scatter(phiRange,phiEst_rep_i(:)*180/pi,5,'k','filled','markerfacealpha',0.1); hold on
                    end
                    xlabel('True \phi (deg.)'); ylabel('Estimate \phi (deg.)');
                    xticks(-90:30:90); yticks(-90:30:90)
                    axis equal; axis tight; grid on; set(gca,'FontSize',fontsize)
                    
                    % check whether angle remapping worked an if error makes sense
                    subplot(1,4,2)
                    for i=1:repeats
                        phiEstRemapped_rep_i = phiRemappedPhoton(:,id_photons,i);
                        scatter(phiRange,phiEstRemapped_rep_i(:)*180/pi,5,'k','filled','markerfacealpha',0.1); hold on
                    end
                    xlabel('True \phi (deg.)'); ylabel('Remapped estimate \phi (deg.)');
                    xticks(-90:30:90); yticks(-90:30:90)
                    axis equal; axis tight; grid on; set(gca,'FontSize',fontsize)
                    
                    % check whether error makes sense
                    subplot(1,4,3)
                    for i=1:repeats
                        phiErrorPhoton_rep_i = phiErrorPhoton(:,id_photons,i);
                        scatter(phiRange,phiErrorPhoton_rep_i(:)*180/pi,5,'k','filled','markerfacealpha',0.1); hold on
                    end
                    xlabel('True \phi (deg.)'); ylabel('Bias \phi (deg.)');
                    xticks(-90:30:90); yticks(-90:30:90); xlim([-90 90]); ylim([-90 90]);
                    axis square; grid on; set(gca,'FontSize',fontsize)
                    
                    subplot(1,4,4)
                    histogram(phiErrorPhoton(:,id_photons,:)*180/pi,-45:0.5:45,'edgecolor','none','facecolor',[0.4 0.4 0.4])
                    xlabel('Bias \phi (deg.)'); ylabel('Occurence');
                    axis square
                end
                fprintf('done!\n')
            end
            
            % plot summary results bias and precision curves as a function of photons
            
            fprintf('Generating summary figure... ')
            figure('position',[100 100 700 220]);
            set(0,'DefaultAxesTitleFontWeight','normal');
            
            linewidth = 1;
            
            subplot(1,2,1)
            scatter(reshape(repmat(photonRange,numel(phiRange),1),1,[]),reshape(bias*180/pi,1,[]),...
                5,'k','filled','markerfacealpha',0.2); hold on
            plot(photonRange,meanBias,'-','color','r','linewidth',linewidth)
            % plot(photons,meanError-stdError,'-k'); hold on
            % plot(photons,meanError+stdError,'-k'); hold on
            xlim([0 inf]); ylim([-20 20])
            xlabel('Photons'); ylabel('Bias \phi (deg.)')
            title(['\theta = ',num2str(Theta),'\circ, \gamma = ',num2str(gamma)])
            set(gca,'FontSize',fontsize); grid on;
            % legend('mean','mean \pm std','min, max')
            legend('data','mean')
            
            subplot(1,2,2)
            scatter(reshape(repmat(photonRange,numel(phiRange),1),1,[]),reshape(precision*180/pi,1,[]),...
                5,'k','filled','markerfacealpha',0.2); hold on
            plot(photonRange,meanPrecision,'-','color','r','linewidth',linewidth)
            % plot([photonRange fliplr(photonRange)],[meanPrecision-stdPrecision fliplr(meanPrecision+stdPrecision)],'-k'); hold on
            % plot([photonRange fliplr(photonRange)],[minPrecision fliplr(maxPrecision)],'--k'); hold on
            xlim([0 inf]); ylim([0 20])
            xlabel('Photons'); ylabel('Precision \phi (deg.)')
            title(['\theta = ',num2str(Theta),'\circ, \gamma = ',num2str(gamma)])
            set(gca,'FontSize',fontsize); grid on;
            % legend('mean','mean \pm std','min, max')
            legend('data','mean')
            fprintf('done!\n')
        end
        
        function [results,simulationsNoise,simulations] = ...
                calculateBiasPrecision_theta_vs_photons(obj,thetaRange,photonRange,phi,repeats,...
                gamma,bkgnd,fov,conversionmethod,w,flagshow)
            %CALCULATEBIASPRECISION_THETA_VS_PHOTONS
            % Note: the molecule will be in focus at z = 0, Delta = 0
            
            % simulate images without noise or background (for current phi)
            Cam_nonoise = Camera('nonoise',obj.Cam.camPixSize,0,obj.Cam.bit,obj.Cam.polariserUnitMosaic);
            locsTable = obj.fixedPhiVaryTheta(thetaRange,phi,gamma,1);
            Sim = Simulation(locsTable,obj.Scope,Cam_nonoise,fov,obj.typeNormalisation);
            [simulations,~] = Sim.simulateSample(0,0);
            
            masks = Sim.masksFov;
            Proc = Process(obj.Scope,obj.Cam,masks,conversionmethod,w);
            
            % normalise and add noise in a loop and analyse
            simulationsNoise = nan(Sim.fov,Sim.fov,numel(thetaRange),numel(photonRange)); % save a few example images with noise
            thetaEst = nan(numel(thetaRange),numel(photonRange),repeats);
            count = 0;
            for id_theta=1:numel(thetaRange)
                for id_photons = 1:numel(photonRange)
                    for id_rep=1:repeats
                        count = count + 1;
                        if ~mod(count,100); fprintf('Analysing frame %d/%d\n',count,numel(thetaRange)*numel(photonRange)*repeats); end
                        
                        frame = abs(photonRange(id_photons)*simulations(:,:,id_theta) + bkgnd); % abs to avoid loosing whole image to 'nan' after adding poisson noise (poissrnd(negative number) gives 'nan')
                        frame = obj.Cam.cameraNoiseModel(frame);
                        if or(strcmp(obj.Cam.cameraType,'nonoise'),strcmp(obj.Cam.cameraType,'ideal'))
                            frame = frame - bkgnd;
                        elseif or(strcmp(obj.Cam.cameraType,'sCMOS'),strcmp(obj.Cam.cameraType,'EMCCD'))
                            frame = frame - bkgnd - obj.Cam.noiseParam.offset;
                        end
                        
                        if flagshow
                            if id_rep == 1
                                simulationsNoise(:,:,id_theta,id_photons) = frame;
                            end
                        end
                        
                        [~,theta,~,~] = Proc.analyseSimulation(frame);
                        thetaEst(id_theta,id_photons,id_rep) = theta;
                    end
                end
            end
            fprintf('done!\n')
            
            % show example images with and without noise
            if flagshow
                figure('Name','Example simulations without noise');
                for i=1:floor(sqrt(numel(thetaRange)))^2
                    subplot(floor(sqrt(numel(thetaRange))),floor(sqrt(numel(thetaRange))),i);
                    imshow(simulations(:,:,i),[])
                end
                for id_photons=1:numel(photonRange)
                    figure('Name',sprintf('Example simulations %d photons',photonRange(id_photons)));
                    for i=1:floor(sqrt(numel(thetaRange)))^2
                        subplot(floor(sqrt(numel(thetaRange))),floor(sqrt(numel(thetaRange))),i);
                        img = simulationsNoise(:,:,i,id_photons);
                        imshow(img,[]); colorbar; title(sprintf('sum(img(:)) = %.1f',sum(img(:))))
                    end
                end
            end
            
            fprintf('Calculating bias and precision... ')
            % calculate bias
            trueTheta = repmat(thetaRange'*pi/180,1,numel(photonRange));
            thetaErrorPhoton = thetaEst - trueTheta;
            bias = mean(thetaErrorPhoton,3);
            meanBias = mean(bias)*180/pi;
            stdBias = std(bias)*180/pi;
            minBias = min(bias)*180/pi;
            maxBias = max(bias)*180/pi;
            
            % calculate precision
            precision = std(thetaErrorPhoton,[],3);
            meanPrecision = mean(precision)*180/pi;
            stdPrecision = std(precision)*180/pi;
            minPrecision = min(precision)*180/pi;
            maxPrecision = max(precision)*180/pi;
            
            % group results
            results.bias          = bias;
            results.meanError     = meanBias;
            results.stdError      = stdBias;
            results.minError      = minBias;
            results.maxError      = maxBias;
            results.precision     = precision;
            results.meanPrecision = meanPrecision;
            results.stdPrecision  = stdPrecision;
            results.minPrecision  = minPrecision;
            results.maxPrecision  = maxPrecision;
            fprintf('done!\n')
            
            %+++++++++++++++++++++
            % Plotting
            fontsize = 10;
            if flagshow
                fprintf('Generating plots... ')
                for id_photons = 1:numel(photonRange)
                    % check by plotting true vs estimate in scatter plot (rawest results)
                    figure('position',[100 450 800 250],'Name',sprintf('%d photons',photonRange(id_photons)));
                    subplot(1,3,1)
                    for i=1:repeats
                        thetaEst_rep_i = thetaEst(:,id_photons,i);
                        scatter(thetaRange,thetaEst_rep_i(:)*180/pi,5,'k','filled','markerfacealpha',0.1); hold on
                    end
                    xlabel('True \theta (deg.)'); ylabel('Estimate \theta (deg.)');
                    xticks(0:15:90); yticks(0:15:90)
                    axis equal; axis tight; grid on; set(gca,'FontSize',fontsize)
                    
                    % check whether error makes sense
                    subplot(1,3,2)
                    for i=1:repeats
                        thetaErrorPhoton_rep_i = thetaErrorPhoton(:,id_photons,i);
                        scatter(thetaRange,thetaErrorPhoton_rep_i(:)*180/pi,5,'k','filled','markerfacealpha',0.1); hold on
                    end
                    xlabel('True \theta (deg.)'); ylabel('Bias \theta (deg.)');
                    xticks(0:15:90); yticks(0:15:90); xlim([0 90]); ylim([0 90]);
                    axis square; grid on; set(gca,'FontSize',fontsize)
                    
                    subplot(1,3,3)
                    histogram(thetaErrorPhoton(:,id_photons,:)*180/pi,-45:0.5:45,'edgecolor','none','facecolor',[0.4 0.4 0.4])
                    xlabel('Bias \theta (deg.)'); ylabel('Occurence');
                    axis square
                end
                fprintf('done!\n')
            end
            
            % plot summary results bias and precision curves as a function of photons
            
            fprintf('Generating summary figure... ')
            figure('position',[100 100 700 220]);
            set(0,'DefaultAxesTitleFontWeight','normal');
            
            linewidth = 1;
            
            subplot(1,2,1)
            scatter(reshape(repmat(photonRange,numel(thetaRange),1),1,[]),reshape(bias*180/pi,1,[]),...
                5,'k','filled','markerfacealpha',0.2); hold on
            plot(photonRange,meanBias,'-','color','r','linewidth',linewidth)
            % plot(photons,meanError-stdError,'-k'); hold on
            % plot(photons,meanError+stdError,'-k'); hold on
            xlim([0 inf]); ylim([-20 20])
            xlabel('Photons'); ylabel('Bias \theta (deg.)')
            title(['\phi = ',num2str(phi),'\circ, \gamma = ',num2str(gamma)])
            set(gca,'FontSize',fontsize); grid on;
            % legend('mean','mean \pm std','min, max')
            legend('data','mean')
            
            subplot(1,2,2)
            scatter(reshape(repmat(photonRange,numel(thetaRange),1),1,[]),reshape(precision*180/pi,1,[]),...
                5,'k','filled','markerfacealpha',0.2); hold on
            plot(photonRange,meanPrecision,'-','color','r','linewidth',linewidth)
            % plot([photonRange fliplr(photonRange)],[meanPrecision-stdPrecision fliplr(meanPrecision+stdPrecision)],'-k'); hold on
            % plot([photonRange fliplr(photonRange)],[minPrecision fliplr(maxPrecision)],'--k'); hold on
            xlim([0 inf]); ylim([0 20])
            xlabel('Photons'); ylabel('Precision \theta (deg.)')
            title(['\phi = ',num2str(phi),'\circ, \gamma = ',num2str(gamma)])
            set(gca,'FontSize',fontsize); grid on;
            % legend('mean','mean \pm std','min, max')
            legend('data','mean')
            fprintf('done!\n')
        end

        function [bias,precision] = plotBiasPrecision_Phi_vs_variedParam(obj,estimates,flag)
            %PLOTBIASPRECISION_PHI_VS_VARIEDPARAM
            fontsize = 10;
            variedParam = estimates.variedParam;
            paramRange  = estimates.variedParamRange;
            thetaRange  = estimates.thetaRange;
            phiRange    = estimates.phiRange;
            repeats     = estimates.repeats;
            
            % calculate bias and precision for phi
            truePhi = repmat(phiRange'*pi/180,1,numel(thetaRange));
            phiError       = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            phiEstRemapped = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            for id_param=1:numel(paramRange)
                for id_reps=1:repeats
                    [phiError_i, phiEstRemapped_i] = obj.remapPhiError(estimates.phi(:,:,id_param,id_reps),truePhi);
                    phiError(:,:,id_param,id_reps) = phiError_i;
                    phiEstRemapped(:,:,id_param,id_reps) = phiEstRemapped_i;
                end
            end
            bias = mean(phiError,4)*180/pi;
            precision = std(phiError,[],4)*180/pi;
            
            if flag.showLandscapes % plot the results
                for i=1:numel(paramRange)
                    fig = figure('name',sprintf('%s = %.2f',estimates.variedParam,estimates.variedParamRange(i)),'position',[100 100 600 230]);
                    
                    subplot(1,2,1)
                    imagesc('XData',thetaRange,'YData',phiRange,'CData',bias(:,:,i)); colorbar;
                    xlim([0 90]); ylim([-90 90]); xticks(0:30:90); yticks(-90:45:90); caxis([-45 45])
                    title(sprintf('Bias phi\n(%s = %.2f)',variedParam,paramRange(i)),'fontweight','normal')
                    xlabel('\theta (\circ)'); ylabel('\phi (\circ)');
                    hcb = colorbar; hcb.Title.String = '\Delta_\phi (\circ)';
                    set(gca,'fontsize',fontsize,'TickDir','out','TickLength',[0.03 0.03])
                    
                    subplot(1,2,2)
                    imagesc('XData',thetaRange,'YData',phiRange,'CData',precision(:,:,i)); colorbar;
                    xlim([0 90]); ylim([-90 90]); xticks(0:30:90); yticks(-90:45:90); caxis([0 45])
                    title(sprintf('Precision phi\n(%s = %.2f)',variedParam,paramRange(i)),'fontweight','normal')
                    xlabel('\theta (\circ)'); ylabel('\phi (\circ)');
                    hcb = colorbar; hcb.Title.String = '\sigma_\phi (\circ)';
                    set(gca,'fontsize',fontsize,'TickDir','out','TickLength',[0.03 0.03])
                    
                    if flag.saveLandscapes
                        figurename = sprintf('landscape_phi_vary_%s_%.1f',variedParam,paramRange(i));
                        saveas(fig,fullfile(flag.outputpath,strcat(figurename,'.fig')))
                        exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.png')),'Resolution',600,'BackgroundColor','none')
                        set(gcf,'renderer','Painters') % always go to painters mode before writing vector graphics
                        exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.eps')))
                        close(fig)
                    end
                end
            end
            
            if flag.showCurves % Get bias and precision phi vs varied param
                linewidth = 1;
                cols = parula(numel(thetaRange)+1);
                fig = figure('position',[100 100 750 230]);
                
                legendNames = {};
                for id_theta=1:numel(thetaRange)
                    legendNames{id_theta} = sprintf('\\theta = %d',thetaRange(id_theta));
                    
                    bias_i = squeeze(mean(bias(:,id_theta,:),1));
                    precision_i = squeeze(mean(precision(:,id_theta,:),1));
                    
                    subplot(1,2,1)
                    plot(paramRange,bias_i,'Color',cols(end - id_theta,:),'LineWidth',linewidth)
                    xlabel(variedParam); ylabel('Bias \phi (\circ)')
                    ylim([-20 20])
                    set(gca,'fontsize',fontsize)
                    hold on; grid on
                    
                    subplot(1,2,2)
                    plot(paramRange,precision_i,'Color',cols(end - id_theta,:),'LineWidth',linewidth)
                    xlabel(variedParam); ylabel('Precision \phi (\circ)')
                    ylim([0 20])
                    set(gca,'fontsize',fontsize)
                    hold on; grid on
                end
                legend(legendNames)
                
                if flag.saveCurves
                    figurename = sprintf('summary_phi_vary_%s_%.1f',variedParam,paramRange(i));
                    saveas(fig,fullfile(flag.outputpath,strcat(figurename,'.fig')))
                    exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.png')),'Resolution',600,'BackgroundColor','none')
                    set(gcf,'renderer','Painters') % always go to painters mode before writing vector graphics
                    exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.eps')))
                    close(fig)
                end
            end
        end
        
        function [bias,precision] = plotBiasPrecision_Theta_vs_variedParam(obj,estimates,flag)
            %PLOTBIASPRECISION_THETA_VS_VARIEDPARAM
            fontsize = 10;
            variedParam = estimates.variedParam;
            paramRange  = estimates.variedParamRange;
            thetaRange  = estimates.thetaRange;
            phiRange    = estimates.phiRange;
            repeats     = estimates.repeats;
            
            % calculate bias and precision for phi
            trueTheta = repmat(thetaRange'*pi/180,1,numel(phiRange))';
            thetaError = nan(numel(phiRange),numel(thetaRange),numel(paramRange),repeats);
            for id_param=1:numel(paramRange)
                for id_reps=1:repeats
                    thetaError(:,:,id_param,id_reps) = estimates.theta(:,:,id_param,id_reps) - trueTheta;
                end
            end
            bias = mean(thetaError,4)*180/pi;
            precision = std(thetaError,[],4)*180/pi;
            
            % plot the results
            if flag.showLandscapes % plot the results
                for i=1:numel(paramRange)
                    fig = figure('name',sprintf('%s = %.2f',estimates.variedParam,estimates.variedParamRange(i)),'position',[100 400 600 230]);
                    
                    subplot(1,2,1)
                    imagesc('XData',thetaRange,'YData',phiRange,'CData',bias(:,:,i)); colorbar;
                    xlim([0 90]); ylim([-90 90]); xticks(0:30:90); yticks(-90:45:90); caxis([-45 45])
                    title(sprintf('Bias theta\n(%s = %.2f)',variedParam,paramRange(i)),'fontweight','normal')
                    xlabel('\theta (\circ)'); ylabel('\phi (\circ)');
                    hcb = colorbar; hcb.Title.String = '\Delta_\theta (\circ)';
                    set(gca,'fontsize',fontsize,'TickDir','out','TickLength',[0.03 0.03])
                    
                    subplot(1,2,2)
                    imagesc('XData',thetaRange,'YData',phiRange,'CData',precision(:,:,i)); colorbar;
                    xlim([0 90]); ylim([-90 90]); xticks(0:30:90); yticks(-90:45:90); caxis([0 45])
                    title(sprintf('Precision theta\n(%s = %.2f)',variedParam,paramRange(i)),'fontweight','normal')
                    xlabel('\theta (\circ)'); ylabel('\phi (\circ)');
                    hcb = colorbar; hcb.Title.String = '\sigma_\theta (\circ)';
                    set(gca,'fontsize',fontsize,'TickDir','out','TickLength',[0.03 0.03])
                    
                    if flag.saveLandscapes
                        figurename = sprintf('landscape_theta_vary_%s_%.1f',variedParam,paramRange(i));
                        saveas(fig,fullfile(flag.outputpath,strcat(figurename,'.fig')))
                        exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.png')),'Resolution',600,'BackgroundColor','none')
                        set(gcf,'renderer','Painters') % always go to painters mode before writing vector graphics
                        exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.eps')))
                        close(fig)
                    end
                end
            end
            
            if flag.showCurves % Get bias and precision phi vs varied param
                
                % Get bias and precision theta vs varied param
                linewidth = 1;
                cols = parula(numel(thetaRange)+1);
                fig = figure('position',[100 100 750 230]);
                
                legendNames = {};
                for i=1:numel(thetaRange)
                    legendNames{i} = sprintf('\\theta = %d',thetaRange(i));
                    
                    bias_i = squeeze(mean(bias(:,i,:),1));
                    precision_i = squeeze(mean(precision(:,i,:),1));
                    
                    subplot(1,2,1)
                    plot(paramRange,bias_i,'Color',cols(end - i,:),'LineWidth',linewidth)
                    xlabel(variedParam); ylabel('Bias \theta (\circ)')
                    ylim([-20 20])
                    set(gca,'fontsize',fontsize)
                    hold on; grid on
                    
                    subplot(1,2,2)
                    plot(paramRange,precision_i,'Color',cols(end - i,:),'LineWidth',linewidth)
                    xlabel(variedParam); ylabel('Precision \theta (\circ)')
                    ylim([0 20])
                    set(gca,'fontsize',fontsize)
                    hold on; grid on
                end
                legend(legendNames)
                
                if flag.saveCurves
                    figurename = sprintf('summary_theta_vary_%s_%.1f',variedParam,paramRange(i));
                    saveas(fig,fullfile(flag.outputpath,strcat(figurename,'.fig')))
                    exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.png')),'Resolution',600,'BackgroundColor','none')
                    set(gcf,'renderer','Painters') % always go to painters mode before writing vector graphics
                    exportgraphics(fig,fullfile(flag.outputpath,strcat(figurename,'.eps')))
                    close(fig)
                end
            end
        end
    end
    
    methods (Static)
        
        function CRLB = getCRLB(stack,dz)
            % Assumes data is a stack and calculates crlb for 3rd dimension
            derivative = diff(stack,1,3);
            F = (1./stack(:,:,2:end)).*((derivative/dz).^2);
            F = sum(sum(F,1),2);
            CRLB = sqrt(1./F(:));
        end

        function locs = fixedPhiVaryTheta(phi,theta,gamma,photons)
            % FIXEDPHIVARYTHETA Generate a localisation list
            % (frame,x,y,z,phi,theta,gamma,photons) to simulate a molecule
            % of some fixed phi, but varying theta.
            theta = theta(:);
            phi   = phi*ones(size(theta));
            locs = table;
            locs.frame   = (1:numel(phi))';
            locs.x       = zeros(numel(phi),1);
            locs.y       = zeros(numel(phi),1);
            locs.z       = zeros(numel(phi),1);
            locs.phi     = phi;
            locs.theta   = theta;
            locs.gamma   = gamma*ones(size(locs.frame));
            locs.photons = photons*ones(size(locs.frame));
            locs.delta   = zeros(numel(phi),1);
        end

        function locs = fixedThetaVaryPhi(phi,theta,gamma,photons)
            % FIXEDTHETAVARYPHI Generate a localisation list
            % (frame,x,y,z,phi,theta,gamma,photons) to simulate a molecule
            % of some fixed theta, but varying phi.
            phi   = phi(:);
            theta = theta*ones(size(phi));
            locs = table;
            locs.frame   = (1:numel(phi))';
            locs.x       = zeros(numel(phi),1);
            locs.y       = zeros(numel(phi),1);
            locs.z       = zeros(numel(phi),1);
            locs.phi     = phi;
            locs.theta   = theta;
            locs.gamma   = gamma*ones(size(locs.frame));
            locs.photons = photons*ones(size(locs.frame));
            locs.delta   = zeros(numel(phi),1);
        end
        
        function [phiError, phiEstRemapped] = remapPhiError(phiEst,phiTrue)
            %REMAPPHIERROR
            % subtract or add an integer multiple of pi from phiEst to minimize error
            candidates = [-2*pi,-pi,0,pi,2*pi];
            distance = zeros(size(phiEst,1),size(phiEst,2),numel(candidates));
            for i=1:length(candidates)
                distance(:,:,i) = abs(phiEst - phiTrue + candidates(i));
            end
            [~,idx] = min(distance,[],3);
            phiEstRemapped = phiEst + candidates(idx);
            phiError = phiEstRemapped - phiTrue;
        end
        
    end
end
