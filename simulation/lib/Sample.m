classdef Sample
    %SAMPLE object contains sample properties and methods to generate
    %ground-truth localisation lists for simulation.
    
    properties
        Microscope % instance of the Microscope class
        fov % field of view, (1x2) array
        %photobleaching, photophysics model
    end
    
    methods
        
        function obj = Sample(Microscope,fov)
            %SAMPLE Construct an instance of this class
            obj.Microscope = Microscope;
            obj.fov = fov;
        end
        
        function locs = immobilisedDyeOnGlassSample(obj,numEmitters,photons,numFrames)
            % IMMOBILISEDDYEONGLASSSAMPLE Generates (frame,x,y,z,phi,theta,gamma,photons)
            % for an immobilised dye sample. The dyes are randomly
            % positioned in XY (Z = 0), have random in-plane orientation,
            % no rotational mobility (gamma = 1) and no out-of-plane component.
            % numEmitters: int, number of dyes
            % fov: (1x2) array, field of view (fov x, fov y)
            % photons: float, number of photons/bead
            % numFrames: int, number of frames
            locs = table;
            frame = repmat((1:numFrames)',1,numEmitters)';
            locs.frame   = frame(:);
            locs.x       = repmat(obj.fov(1)*rand(1,numEmitters),1,numFrames)' - obj.fov(1)/2;
            locs.y       = repmat(obj.fov(2)*rand(1,numEmitters),1,numFrames)' - obj.fov(2)/2;
            locs.z       = zeros(numEmitters*numFrames,1);
            locs.phi     = repmat(rand(1,numEmitters)*pi - pi/2,1,numFrames)';
            locs.theta   = ones(numEmitters*numFrames,1)*pi/2;
            locs.gamma   = ones(numEmitters*numFrames,1); % immobilised, no rotational mobility
            locs.photons = photons*ones(numEmitters*numFrames,1);
        end
        
        function locs = beadSample(obj,numBeads,photons,numFrames)
            %BEADSAMPLE Generates (frame,x,y,z,photons) a bead sample
            % numBeads: int, number of beads
            % fov: (1x2) array, field of view (fov x, fov y)
            % photons: float, number of photons/bead
            % numFrames: int, number of frames
            locs = table;
            frame = repmat((1:numFrames)',1,numBeads)';
            locs.frame   = frame(:);
            locs.x       = repmat(obj.fov(1)*rand(1,numBeads),1,numFrames)' - obj.fov(1)/2;
            locs.y       = repmat(obj.fov(2)*rand(1,numBeads),1,numFrames)' - obj.fov(2)/2;
            locs.z       = zeros(numBeads*numFrames,1);
            locs.phi     = zeros(numBeads*numFrames,1);
            locs.theta   = zeros(numBeads*numFrames,1);
            locs.gamma   = zeros(numBeads*numFrames,1); % isotropic
            locs.photons = photons*ones(numBeads*numFrames,1);
        end
        
    end    
    
    methods (Static)
        
        function locs = gridPointsOnHalfHemisphere(dAngle,gamma,photons)
            % GRIDPOINTSONHALFHEMISPHERE Generate a localisation list
            % (frame,x,y,z,phi,theta,gamma,photons) to simulate a range of
            % orientations regularly spaced on a hemisphere.
            phi   = (0:dAngle:180-dAngle)*pi/180;
            theta = linspace(0,90,90/dAngle)*pi/180;
            [theta,phi] = meshgrid(theta,phi);
            phi = phi(:);
            theta = theta(:);

            locs = table;
            locs.frame   = (1:numel(phi))';
            locs.x       = zeros(size(locs.frame));
            locs.y       = zeros(size(locs.frame));
            locs.z       = zeros(size(locs.frame));
            locs.phi     = phi;
            locs.theta   = theta;
            locs.gamma   = gamma*ones(size(locs.frame));
            locs.photons = photons*ones(size(locs.frame));
            locs.delta   = zeros(size(locs.frame));
        end
        
        function locs = gridPointsOnFullHemisphere(dAngle,gamma,photons)
            % GRIDPOINTSONFULLHEMISPHERE Generate a localisation list
            % (frame,x,y,z,phi,theta,gamma,photons) to simulate a range of
            % orientations regularly spaced on a hemisphere.
            phi   = (0:dAngle:360-dAngle)*pi/180;
            theta = linspace(0,90,90/dAngle)*pi/180;
            [theta,phi] = meshgrid(theta,phi);
            phi = phi(:);
            theta = theta(:);

            locs = table;
            locs.frame   = 1:numel(phi);
            locs.x       = zeros(size(locs.frame));
            locs.y       = zeros(size(locs.frame));
            locs.z       = zeros(size(locs.frame));
            locs.phi     = phi;
            locs.theta   = theta;
            locs.gamma   = gamma*ones(size(locs.frame));
            locs.photons = photons*ones(size(locs.frame));
        end
        
    end
end
