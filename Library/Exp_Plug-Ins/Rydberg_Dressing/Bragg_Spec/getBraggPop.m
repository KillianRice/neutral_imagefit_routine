function [indVar,braggFrac,totPop,expTime,label,numPeak] = getBraggPop(analyVar,indivDataset,avgDataset)
% Function to find the population of the Bragg peaks. If there was enough room to specify separate windows and
% fit everything at once this function will simply return the number already calculated. However, in the case
% of shallow angle Bragg Spectroscopy the kicked population cannot be clearly spatially resolved. To overcome
% this challenge we will subtract the main condensate fit from the 2D image and numerically integrate the
% signal in boxes defined around the kicked populations.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit routines
%   indivDataset - Cell of structures containing all scan/batch specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%   indVar   - Cell of independent variables for each scan
%   braggPop - Cell of bragg populations in each scan
%   totPop   - Cell of the total population in each scan
%   expTime  - Exposure time of each scan
%   label    - Label used to enumerate scans in plots

%% Decide to use averaged data or individual scans
avgAutoFlag = length(cell2mat(analyVar.posOccurUniqVar)) > length(analyVar.posOccurUniqVar) & analyVar.uniqScanList ~= 0;

%% Extract independent variables and number from correct structure
indField   = {'imagevcoAtom' 'simScanIndVar'};  % fields of independent variables in indivDataset and avgDataset
labelField = {'timevectorAtom' 'uniqScanList'}; % label used for each plot

% Selected datatype picks the data to be analyzed
if avgAutoFlag
    curDataset = avgDataset;   fieldIndx = 2;
else
    curDataset = indivDataset; fieldIndx = 1;
end

% Save variables from cell of structures
indVar  = cellfun(@(x) analyVar.funcDataScale(x.(indField{fieldIndx})),curDataset,'UniformOutput',0);
label   = analyVar.(labelField{fieldIndx});


%% Get peak populations
if analyVar.braggAngle == 180
    [braggFrac,totPop]  = largeAngleBragg(analyVar,indivDataset,avgDataset);
elseif analyVar.braggAngle == 21
    [braggFrac,totPop]  = shallowAngleBragg(analyVar,indivDataset,avgDataset);
end

%% Scale Bragg Population by recorded power if required
if analyVar.scaleIntFlag && ~sum(avgAutoFlag) % averaged data is scaled before averaging
    for basenameNum = 1:analyVar.numBasenamesAtom
        % Scale by Bragg intensities\
        braggFrac{basenameNum} = intenNormNum(analyVar,indivDataset,basenameNum,braggFrac{basenameNum});
    end 
end

%% Save Exposure times for each scan (or group)
if avgAutoFlag
    % Initialize loop variables
    expTime = nan(length(analyVar.uniqScanList),1);
    for uniqScanIter = 1:length(analyVar.uniqScanList);
        % Get exposure times of all members of the group
        tmpTime = analyVar.expTime(analyVar.posOccurUniqVar{uniqScanIter});
        
        % Check that the window number is the same for all scans in a group
        if range(tmpTime) ~= 0
            error('expTime is not the same for each member of the group labeled %g',analyVar.uniqScanList(uniqScanIter))
        end
        
        % Save exposure time 
        expTime(uniqScanIter) = tmpTime(1)*1e-3;
    end
else
    expTime = analyVar.expTime*1e-3;
end      

%% Save Exposure times for each scan (or group)
if avgAutoFlag
    % Initialize loop variables
    numPeak = nan(length(analyVar.uniqScanList),1);
    for uniqScanIter = 1:length(analyVar.uniqScanList);
        % Get window ID of all members of the group
        tmpPeak = analyVar.winID(analyVar.posOccurUniqVar{uniqScanIter});
        
        % Save number of peaks
        numPeak(uniqScanIter) = abs(tmpPeak(1));
    end
else
    numPeak = abs(analyVar.winID);
end   
end