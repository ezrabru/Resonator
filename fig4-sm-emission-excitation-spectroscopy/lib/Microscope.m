classdef Microscope
    %MICROSCOPE object contains microscope properties.
    
    properties
        NA % numerical aperture of the objective
        magnification % total system magnification
        nImmersion % refractive index immersion medium (e.g. 1.518 for oil)
        camOffset
        camGain
        camQE
        camPixSize

        virtualPixelSize % virtual pixel size in nm
    end
    
    methods
        
        function obj = Microscope(NA,magnification,nImmersion,camOffset,camGain,camQE,camPixSize)
            %MICROSCOPE Construct an instance of this class
            obj.NA = NA;
            obj.magnification = magnification;
            obj.nImmersion = nImmersion;
            obj.camOffset = camOffset;
            obj.camGain = camGain;
            obj.camQE = camQE;
            obj.camPixSize = camPixSize;
            
            obj.virtualPixelSize = (1e9*obj.camPixSize)/obj.magnification;
        end
        
    end
    
    methods (Static)
        
        %...
        
    end
    
end
