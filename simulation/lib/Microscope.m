classdef Microscope
    %MICROSCOPE object contains microscope properties.
    
    properties
        detectionType % normal, 2Cpolarised, 4Cpolarised, polcam
        NA % numerical aperture of the objective
        fObj % focal length of the objective lens (fObj = M/fTubeLens)
        fTubeLens % focal length of the tube lens
        nImmersion % refractive index of the cover glass
        nMedium % refractive index of the medium in which the molecule has been embedded
        wavelength % wavelength of the emission
                
        numPixBFP % number of pixels in the BFP (for simulation purposes)
        
        k0 % wavenumber
        rhoMax % maximum collection angle of the objective lens
        bfpRadiusNative % radius of the native back focal plane
        bfpDiameterNative % diameter of the native back focal plane
        
        BFPphi
        BFPrho
        
        BFPapertureNative
        
        magnificationObjective % magnificationof the objective
    end
    
    methods
        
        function obj = Microscope(detectionType,NA,fObj,fTubeLens,...
                nImmersion,nMedium,wavelength,numPixBFP)
            %MICROSCOPE Construct an instance of this class
            obj.detectionType = detectionType;
            obj.NA = NA;
            obj.fObj = fObj;
            obj.fTubeLens = fTubeLens;
            obj.nImmersion = nImmersion;
            obj.nMedium = nMedium;
            obj.wavelength = wavelength;
                        
            obj.numPixBFP = numPixBFP;
            
            obj.k0 = 2*pi/wavelength;
            obj.rhoMax = getRhoMax(obj);
            obj.bfpRadiusNative = getNativeBFPRadius(obj);
            obj.bfpDiameterNative = getNativeBFPDiameter(obj);
            
            [BFPphi,BFPrho] = getPolarCoordinatesBFP(obj);
            obj.BFPphi = BFPphi;
            obj.BFPrho = BFPrho;
            
            obj.BFPapertureNative = getNativeBFPAperture(obj);
            
            obj.magnificationObjective = fTubeLens/fObj;
        end
        
        function rhoMax = getRhoMax(obj)
            %GETRHOMAX Calculate the maximum collection angle of the
            %objective lens (from NA definition and the Abbe sine condition)
            rhoMax = obj.NA/obj.nMedium;
        end
        
        function bfpRadiusNative = getNativeBFPRadius(obj)
            %GETNATIVEBFPRADIUS Calculate the radius of the native back focal plane
            bfpRadiusNative = obj.fObj*obj.NA;
        end
        
        function bfpDiameterNative = getNativeBFPDiameter(obj)
            %GETNATIVEBFPDIAMETER Calculate the radius of the native back focal plane
            bfpDiameterNative = 2*obj.fObj*obj.NA;
        end
        
        function [phi,rho] = getPolarCoordinatesBFP(obj)
            %GETPOLARCOORDINATESBFP Generate meshgrid in polar coordinates
            %(phi,rho) to sample the BFP
            x = linspace(-obj.rhoMax,obj.rhoMax,obj.numPixBFP);
            [x,y] = meshgrid(x,x);
            [phi,rho] = cart2pol(x,y);
        end
        
        function nativeBFPaperture = getNativeBFPAperture(obj)
            %GETNATIVEBFPAPERTURE Calculate the aperture (binary matrix) of
            %the native back focal plane
            nativeBFPaperture = obj.BFPrho < obj.rhoMax;
        end
        
    end    
    
    methods (Static)
        
        %...
        
    end
    
end
