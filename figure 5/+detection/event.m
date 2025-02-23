classdef event
    % EVENT Class to store information about detected events.
    % This class encapsulates the start times, stop times, and peak indices
    % of events detected in a signal (e.g., ripple events in neural data).
    %
    % Author: AmirMohammad Azarmehri

    properties
        starts  % Vector containing the start indices of events
        stops   % Vector containing the stop indices of events
        peaks   % Vector containing the peak indices of events
    end
    
    methods
        function obj = event(starts, stops, peaks)
            % EVENT Construct an instance of the event class.
            %   Inputs:
            %       starts - vector of start indices for each event
            %       stops  - vector of stop indices for each event
            %       peaks  - vector of peak indices for each event
            obj.starts = starts;
            obj.stops = stops;
            obj.peaks = peaks;
        end
    end
end
