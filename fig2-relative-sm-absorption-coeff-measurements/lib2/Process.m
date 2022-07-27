classdef Process
    %PROCESS object contains properties and methods to process
    %Resonator image data.
    
    properties
        Microscope % instance of the Microscope class
        bkgndEstimation % method for background estimation
        localisationParams % localisation parameters
        
        border
    end
    
    methods
        
        function obj = Process(Microscope,bkgndEstimation,localisationParams)
            %PROCESS Construct an instance of this class
            obj.Microscope = Microscope;
            obj.bkgndEstimation = bkgndEstimation;
            obj.localisationParams = localisationParams;
            
            obj.border = localisationParams.w + 1;
        end
        
        function locs = processFrame(obj,img,bkgnd)
            % PROCESSFRAME
            
            img = (img - bkgnd - obj.Microscope.camOffset)./(obj.Microscope.camGain*obj.Microscope.camQE); % convert to photons
            
            [x,y] = find2DlocalisationCandidates(obj,img);
            figure; imshow(img,[]); colorbar; hold on; scatter(y,x,'or')
            
            locs = fit2Dlocalisations(obj,img,x,y);
            hold on; scatter([locs.y]*1e-9,[locs.x]*1e-9,'xg')
            
%             locs = groupLocalisations()

            if isempty(locs); locs = nan; end
        end

        function bkgnd = estimateBackground(obj,filepath)
            % ESTIMATEBACKGROUND Different methods to estimate
            % fluorescence background in the four polarised channels.
            
            warning('off','imageio:tiffmexutils:libtiffWarning'); % turn off a warning caused by format thorcam tif files that blows up the command window
            switch obj.bkgndEstimation.method

                case 'none'; bkgnd = 0;
                
                case 'min of avg'
                    avg = obj.getAverageIntensityProjectionFromPath(filepath) - obj.Microscope.camOffset;
                    bkgnd = min(avg(:));

                case 'medfilt mip'
                    mip = obj.getMinimumIntensityProjectionFromPath(filepath) - obj.Microscope.camOffset;
                    params = obj.bkgndEstimation.params;
                    bkgnd = medfilt2(mip,[params.k params.k],'symmetric');
                    
                case 'medfilt avg'
                    avg = obj.getAverageIntensityProjectionFromPath(filepath) - obj.Microscope.camOffset;
                    params = obj.bkgndEstimation.params;
                    bkgnd = medfilt2(avg,[params.k params.k],'symmetric');

                case 'medfilt med'
                    tr = Tiff(filepath,'r'); % set up object for reading frame by frame
                    params = obj.bkgndEstimation.params;
                    for i=1:params.frames
                        tr.nextDirectory();
                        img = double(tr.read()) - obj.Microscope.camOffset;
                        if i==1
                            bkgnd = nan(size(img,1),size(Q00,2),params.frames);
                        end
                        bkgnd(:,:,i) = img;
                    end
                    % temporal median projection
                    bkgnd = median(bkgnd,3);
                    % spatial median filter
                    bkgnd = medfilt2(bkgnd,[params.k params.k],'symmetric');

                otherwise; disp('Unexpected value for variable "methodBkgndEstimation".'); return

            end
        end

        function [x,y] = find2DlocalisationCandidates(obj,img)
            % FIND2DLOCALISATIONCANDIDATES find local maxima in an image after a
            % difference-of-gaussian (DoG) filter.
            % 
            % DoG = (img guassian blurred with sigma1) - (img guassian blurred with sigma2)
            % 
            % input:
            %   img    - array, the image
            %   sigma1 - double, standard deviation of the first gaussian filter
            %   sigma1 - double, standard deviation of the second gaussian filter
            %   DoGmin - double, ignore all maxima in regions where DoG < DoGmin
            % output:
            %   x - nx1 array, the x-coordinates of identified local maxima (pixel resolution)
            %   y - nx1 array, the y-coordinates of identified local maxima (pixel resolution)

            % pre-processing filter
            imgMed = medfilt2(img,[3,3],'symmetric');
            G1 = imgaussfilt(imgMed,obj.localisationParams.DoGsigmaSmall);
            G2 = imgaussfilt(imgMed,obj.localisationParams.DoGsigmaLarge);
            DoG = G1 - G2;

            % find local maxima (localisation candidates)
            DoG(DoG < obj.localisationParams.DoGminThreshold) = 0;
            BW = imregionalmax(DoG);
            [x,y] = find(BW); % get coordinates of local maxima
            
            % remove candidates that are too close to image borders
            keep_x = (x > obj.border).*(x < size(BW,2) - obj.border);
            keep_y = (y > obj.border).*(y < size(BW,1) - obj.border);
            keep = logical(keep_x.*keep_y);
            x = x(keep);
            y = y(keep);
        end
        
        function locs = fit2Dlocalisations(obj,img,x,y)
            % FIT2DLOCALISATIONS
                                    
            % create meshgrid for centroid calculation in small 2*w+1 square roi
            x_mesh = -obj.localisationParams.w:obj.localisationParams.w;
            [x_mesh,y_mesh] = meshgrid(x_mesh,x_mesh);
            
            locs = {};
            switch obj.localisationParams.method
                
                case 'centroid'
                    for i=1:length(x)
                        roi = obj.cropToROI(img,x(i),y(i),obj.localisationParams.w);
                        roi = roi - min(roi(:));
                        intensity = sum(roi(:));
                        x_centroid = sum(x_mesh.*roi,'all')/intensity;
                        y_centroid = sum(y_mesh.*roi,'all')/intensity;
                        locs(i).intensity = intensity/2;
                        locs(i).x = (x(i) + y_centroid)*obj.Microscope.virtualPixelSize;
                        locs(i).y = (y(i) + x_centroid)*obj.Microscope.virtualPixelSize;
                    end
                    
                case 'symmetric gaussian'
                    for i=1:length(x)
                        roi = obj.cropToROI(img,x(i),y(i),obj.localisationParams.w);
                        intensity = sum(roi(:) - min(roi(:)));
                        [params,~,~] = obj.fit2DGaussian(roi,obj.localisationParams.w,'symmetric');
                        locs(i).intensity = intensity/2;
                        locs(i).amplitude = params.amplitude;
                        locs(i).x = (x(i) + params.y0)*obj.Microscope.virtualPixelSize;
                        locs(i).y = (y(i) + params.x0)*obj.Microscope.virtualPixelSize;
                        locs(i).sigma = params.sigma;
                        locs(i).bkgnd = params.bkgnd/2;
                        locs(i).photons = (2*pi*params.amplitude.*params.sigma.*params.sigma)/2;
                        
                        a = obj.Microscope.virtualPixelSize;
                        s = a*locs(i).sigma;
                        b = locs(i).bkgnd;
                        N = locs(i).photons;
                        [uncertainty_xy, uncertainty_photons] = obj.getThompsonUncertainty(s,a,b,N);
                        locs(i).uncertainty_xy = uncertainty_xy;
                        locs(i).uncertainty_photons = uncertainty_photons;
                    end
                    
                case 'asymmetric gaussian'
                    for i=1:length(x)
                        roi = obj.cropToROI(img,x(i),y(i),obj.localisationParams.w);
                        intensity = sum(roi(:) - min(roi(:)));
                        [params,~,~] = obj.fit2DGaussian(roi,obj.localisationParams.w,'asymmetric');
                        locs(i).intensity = intensity/2;
                        locs(i).amplitude = params.amplitude;
                        locs(i).x = (x(i) + params.y0)*obj.Microscope.virtualPixelSize;
                        locs(i).y = (y(i) + params.x0)*obj.Microscope.virtualPixelSize;
                        locs(i).sigmax = params.sigmax;
                        locs(i).sigmay = params.sigmay;
                        locs(i).bkgnd = params.bkgnd/2;
                        locs(i).photons = (2*pi*params.amplitude.*params.sigmax.*params.sigmay)/2;
                        
                        a = obj.Microscope.virtualPixelSize;
                        s = a*(locs(i).sigmax + locs(i).sigmay)/2;
                        b = locs(i).bkgnd;
                        N = locs(i).photons;
                        [uncertainty_xy, uncertainty_photons] = obj.getThompsonUncertainty(s,a,b,N);
                        locs(i).uncertainty_xy = uncertainty_xy;
                        locs(i).uncertainty_photons = uncertainty_photons;
                    end
                    
                case 'rotated asymmetric gaussian'
                    for i=1:length(x)
                        roi = obj.cropToROI(img,x(i),y(i),obj.localisationParams.w);
                        intensity = sum(roi(:) - min(roi(:)));
                        [params,~,~] = obj.fit2DGaussian(roi,obj.localisationParams.w,'rotated asymmetric');
                        locs(i).intensity = intensity/2;
                        locs(i).amplitude = params.amplitude;
                        locs(i).x = (x(i) + params.y0)*obj.Microscope.virtualPixelSize;
                        locs(i).y = (y(i) + params.x0)*obj.Microscope.virtualPixelSize;
                        locs(i).sigmax = params.sigmax;
                        locs(i).sigmay = params.sigmay;
                        locs(i).rot = params.theta;
                        locs(i).bkgnd = params.bkgnd/2;
                        locs(i).photons = (2*pi*params.amplitude.*params.sigmax.*params.sigmay)/2;
                        locs(i).sigmaRatio = params.sigmax./params.sigmay;
                        
                        a = obj.Microscope.virtualPixelSize;
                        s = a*(locs(i).sigmax + locs(i).sigmay)/2;
                        b = locs(i).bkgnd;
                        N = locs(i).photons;
                        [uncertainty_xy, uncertainty_photons] = obj.getThompsonUncertainty(s,a,b,N);
                        locs(i).uncertainty_xy = uncertainty_xy;
                        locs(i).uncertainty_photons = uncertainty_photons;
                    end
                    % get orientation estimate from orientation fit
                    angleGaus = params.theta;
                    angleGaus(params.sigmax./params.sigmay > 1) = angleGaus(params.sigmax./params.sigmay > 1) + pi/2;
                    angleGaus = (angleGaus - pi/2)*180/pi;
                    angleGaus(angleGaus < - 90) = angleGaus(angleGaus < - 90) + 180;
                    locs(i).angleGaus = angleGaus*pi/180;
                    
                otherwise; fprintf('Unexpected value for "method" in function "fit2Dlocalisations".')
            end
            
        end
                
    end
    
    methods (Static)
        
        function crop = cropToROICentered(img,fov)
            %CROPTOROICENTERED
            [nx,ny] = size(img);
            xc = ceil(nx/2);
            yc = ceil(ny/2);
            w = floor(fov/2);
            crop = img(xc-w:xc-w+fov-1,...
                yc-w:yc-w+fov-1,:);
        end
        
        function crop = cropToROI(img,xc,yc,w)
            %CROPTOROI
            crop = img(xc-w:xc+w,yc-w:yc+w,:);
        end

        function [params,resnorm,residual] = fit2DGaussian(img,w,type)
            % FITSYMMETRIC2DGAUSSIAN fits a symmetric 2d gaussian to a square image.
            % input:
            %   img  - square array of size [2*w+1, 2*w+1], the image
            %   w    - int, as defined above (could be calculated each time, but as it usually is
            %               a constant variable it can just be supplied to the function)
            %   type - string, 'symmetric', 'asymmetric' or 'rotated asymmetric' type
            %               of 2d gaussian
            % output:
            %   params - ...
            %   amplitude - double, the amplitude of the fitted gaussian
            %   x0 and y0 - double, the x and y coordinate of the center of the
            %               gaussian, relative to the center of the middle pixel which
            %               is at (0,0)
            %   sigma     - double, the standard deviation of the fitted gaussian
            %   resnorm   - double, sum of the square of the residuals, i.e. sum(residual(:).^2)
            %   residual  - array, the residual of the fit, same size as img

            % get (x,y) meshgrid
            x = -w:w;
            [x,y] = meshgrid(x,x);
            xdata = cat(3,x,y);

            switch type

                case 'symmetric'
                    % set lower and upper bounds and initial values of x = [amplitude,xo,yo,sigma]
                    lb = [0,-w,-w,0,-inf];
                    ub = [realmax('double'),w,w,w^2,inf];
                    x0 = [max(img(:)),0,0,1,min(img(:))]; % initial values
                    options = optimset('Display','off'); % don't flood command window with messages
                    
                    % symmetric 2D gaussian
                    f = @(x,xdata)x(1)*exp(-( (xdata(:,:,1)-x(2)).^2  + ...
                                              (xdata(:,:,2)-x(3)).^2  )/(2*x(4).^2)) + x(5);
                    [x,resnorm,residual,~] = lsqcurvefit(f,x0,xdata,img,lb,ub,options);
                    params.amplitude = x(1);
                    params.x0        = x(2);
                    params.y0        = x(3);
                    params.sigma     = x(4);
                    params.bkgnd     = x(5);

                case 'asymmetric'
                    % set lower and upper bounds and initial values of x = [amplitude,xo,yo,sigmax,sigmay]
                    lb = [0,-w,-w,0,0,-inf];
                    ub = [realmax('double'),w,w,w^2,w^2,inf];
                    x0 = [max(img(:)),0,0,1,1,min(img(:))]; % initial values
                    options = optimset('Display','off'); % don't flood command window with messages
                    
                    % asymmetric 2D gaussian
                    f = @(x,xdata)x(1)*exp(-( ((xdata(:,:,1)-x(2)).^2)/(2*x(4).^2) + ...
                                              ((xdata(:,:,2)-x(3)).^2)/(2*x(5).^2) )) + x(6);
                    [x,resnorm,residual,~] = lsqcurvefit(f,x0,xdata,img,lb,ub,options);
                    params.amplitude = x(1);
                    params.x0        = x(2);
                    params.y0        = x(3);
                    params.sigmax    = x(4);
                    params.sigmay    = x(5);
                    params.bkgnd     = x(6);

                case 'rotated asymmetric'
                    % set lower and upper bounds and initial values of x = [amplitude,xo,yo,sigmax,sigmay,theta]
                    lb = [0,-w,-w,0,0,-pi/4,-inf];
                    ub = [realmax('double'),w,w,w^2,w^2,pi/4,inf];
                    x0 = [max(img(:)),0,0,1,1,0,min(img(:))]; % initial values
                    options = optimset('Display','off'); % don't flood command window with messages
                    
                    % asymmetric 2d gaussian but using rotated coordinates
                    f = @(x,xdata)x(1)*exp(-(((xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6))) - (x(2)*cos(x(6)) - x(3)*sin(x(6)))).^2/(2*x(4)^2) + ...
                                             ((xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6))) - (x(2)*sin(x(6)) + x(3)*cos(x(6)))).^2/(2*x(5)^2))) + x(7);
                    
                    [x,resnorm,residual,~] = lsqcurvefit(f,x0,xdata,img,lb,ub,options);
                    params.amplitude = x(1);
                    params.x0        = x(2);
                    params.y0        = x(3);
                    params.sigmax    = x(4);
                    params.sigmay    = x(5);
                    params.theta     = x(6);
                    params.bkgnd     = x(7);

                otherwise
                    disp('Unexpected value for "type" in function "fit2DGaussian".')
            end

        end

        function F = symmetric2DGaussian(x,xdata)
            % SYMMETRIC2DGAUSSIAN function describing a symmetric 2d gaussian used as
            % input for lsqcurvefit during fitting.
            % input:
            %   x     - 1x4 array containing parameters [amplitude, x0, y0, sigma]
            %   xdata - NxNx2 matrix, containing spatial x and y values at which to
            %           calculate the function
            % output:
            %   F - NxN matrix, the calcualted 2d gaussian
            F = x(1)*exp(-(  (xdata(:,:,1)-x(2)).^2  + ...
                             (xdata(:,:,2)-x(3)).^2  )/(2*x(4).^2));
        end

        function F = asymmetric2DGaussian(x,xdata)
            % ASYMMETRIC2DGAUSSIAN function describing an asymmetric 2d gaussian used as
            % input for lsqcurvefit during fitting.
            % input:
            %   x     - 1x5 array containing parameters [amplitude, x0, y0, sigmax, sigmay]
            %   xdata - NxNx2 matrix, containing spatial x and y values at which to
            %           calculate the function
            % output:
            %   F - NxN matrix, the calcualted 2d gaussian
            F = x(1)*exp(-(  ((xdata(:,:,1)-x(2)).^2)/(2*x(4).^2) + ...
                             ((xdata(:,:,2)-x(3)).^2)/(2*x(5).^2)  ));
        end

        function F = rotatedAsymmetric2DGaussian(x,xdata)
            % ROTATEDASYMMETRIC2DGAUSSIAN function describing a rotated asymmetric 2d
            % gaussian used as input for lsqcurvefit during fitting.
            % input:
            %   x     - 1x6 array containing parameters [amplitude, x0, y0, sigmax, sigmay, theta]
            %   xdata - NxNx2 matrix, containing spatial x and y values at which to
            %           calculate the function
            % output:
            %   F - NxN matrix, the calcualted 2d gaussian

            % rotate x and y coordinate space by angle theta = x(6) around origin
            xdatarot(:,:,1) = xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
            xdatarot(:,:,2) = xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));

            % rotate centroid position by angle theta = x(6) around origin
            x0rot = x(2)*cos(x(6)) - x(3)*sin(x(6));
            y0rot = x(2)*sin(x(6)) + x(3)*cos(x(6));

            % asymmetric 2d gaussian but using rotated coordinates
            F = x(1)*exp(-((xdatarot(:,:,1)-x0rot).^2/(2*x(4)^2) + ...
                           (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2)));
        end
        
        function mip = getMinimumIntensityProjectionFromPath(filepath)
            % GETMINIMUMINTENSITYPROJECTIONFROMPATH Get a minimum intensity
            % projection of a tiff stack from file.
            tr = Tiff(filepath,'r'); % set up object for reading frame by frame
            mip = double(tr.read()); % read first frame
            if tr.lastDirectory()
                return
            else
                tr.nextDirectory();
                while true % read remaining frames
                    mip = min(cat(3,mip,double(tr.read())),[],3);
                    if tr.lastDirectory(); break;
                    else; tr.nextDirectory();
                    end
                end
            end
            close(tr)
        end
        
        function avg = getAverageIntensityProjectionFromPath(filepath)
            % GETAVERAGEINTENSITYPROJECTIONFROMPATH Get a minimum intensity
            % projection of a tiff stack from file.
            tr = Tiff(filepath,'r'); % set up object for reading frame by frame
            count = 1;
            avg = double(tr.read()); % read first frame
            if tr.lastDirectory()
                return
            else
                tr.nextDirectory();
                while true % read remaining frames
                    avg = avg + double(tr.read());
                    count = count + 1;
                    if tr.lastDirectory(); break;
                    else; tr.nextDirectory();
                    end
                end
            end
            close(tr)
            avg = avg/count;
        end
        
        function [uncertainty_xy, uncertainty_photons] = getThompsonUncertainty(s,a,b,N)
            % GETTHOMPSONUNCERTAINTY Calculate the theoretical localisation
            % uncertainty and photon number estimate uncertainty using the
            % equations 17 and 19 from Thompson, Larson and Webb, Biophys J
            % 82(5):2775-2783, 2002
            % s: sigma width of fitted gaussian (in nm)
            % a: virtual pixel size (in nm)
            % b: background/pixel in photons
            % N: total number of detected signal photons (without background)
            
            % uncertainty xy localisation
            term1 = (s.^2 + (a.^2/12))./N;
            term2 = (8*pi*(s.^4).*(b.^2))./(N*a).^2;
            uncertainty_xy = sqrt(term1 + term2);
            
            % uncertainty intensity estimate
            term1 = N;
            term2 = (4*pi*(s.^2).*(b.^2))./(a.^2);
            uncertainty_photons = sqrt(term1 + term2);
        end
        
    end
end
