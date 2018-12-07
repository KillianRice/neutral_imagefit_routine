function [onInd] = scopeTimeInd(expTime,timeVec,timeInc)
% Function to serve as a lookup table for given exposure times. 
% This function will return the points where the offset value and 
% intensity values are calculated.
%
% INTPUTS:
%   expTime - Set exposure time for scan
%   timeVec - Time vector for the scan currently being analyzed
%   timeInc - Time increment of the data in the time vector
%
% OUTPUTS:
%   onInd  - Indices of points used to calculate the mean photodiode reading

%% Initialize variables
percDrop        = .01; % Percentage of the total window for dropping area around t = 0;
[offInd, onInd] = deal(zeros(1,length(timeVec)));

%% Calculate size of buffer regions
dropBuffer = floor(length(timeVec)*percDrop);

%% Find where timeVec crosses zero
dTime   = floor(abs(timeVec/(10^(floor(log(abs(timeInc))./log(10))))));

% %% Get the off time indices
% offStart = 1;
% offEnd   = find(dTime == min(dTime),1,'first') - dropBuffer;
% 
% offInd(offStart:offEnd) = 1;
% offInd = logical(offInd);

%% Get the on time indices
onStart = find(dTime == min(dTime),1,'last') + dropBuffer;
onEnd   = length(timeVec);
% If exposure time ends within the window then redefine onEnd
onNumPnts = expTime/timeInc - dropBuffer;
if onNumPnts <= length(timeVec) - dropBuffer;
    onEnd = onStart + onNumPnts - dropBuffer;
end

onInd(onStart:onEnd) = 1;
onInd = logical(onInd);

%% In case special cases are needed
switch expTime 
end
end