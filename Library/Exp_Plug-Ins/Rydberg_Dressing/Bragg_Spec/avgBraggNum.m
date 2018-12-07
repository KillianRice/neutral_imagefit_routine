function  [avgBraggFrac, avgTotNum] = avgBraggNum(analyVar,indivDataset,avgDataset,uniqScanIter,braggNum,totNum)
% This function will normalize and average the atom number when applying
% averaging to Bragg spectroscopy
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit routines
%   indivDataset - Cell of structures containing all scan/batch specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%   uniqScanIter - 
%   braggNum     - 
%   totPop       - 

%% Preallocate loop variables
[avgBraggFrac, avgTotNum]...
    = deal(NaN(length(avgDataset{uniqScanIter}.simScanIndVar),...
    length(analyVar.posOccurUniqVar{uniqScanIter})));

%% Loop through all scans that share the current value of the averaging variable
for simScanIter = 1:length(analyVar.posOccurUniqVar{uniqScanIter})
    % Assign the current scan to be opened
    basenameNum = analyVar.posOccurUniqVar{uniqScanIter}(simScanIter);
    
    % Find the intersection of the scanned variable with the list of all possible values
    [~,idxSharedWithSim,idxSharedInBatch] = intersect(avgDataset{uniqScanIter}.simScanIndVar,...
        indivDataset{basenameNum}.imagevcoAtom);
    
    % Find number in the peaks of interest
    avgBraggFrac(idxSharedWithSim,simScanIter) = braggNum{simScanIter}(idxSharedInBatch)./totNum{simScanIter}(idxSharedInBatch);
    avgTotNum(idxSharedWithSim,simScanIter)    = totNum{simScanIter}(idxSharedInBatch);
    
    % Scale Bragg intensities
    if analyVar.scaleIntFlag
        avgBraggFrac(sort(idxSharedWithSim),simScanIter) = intenNormNum(analyVar,indivDataset,basenameNum,avgBraggFrac(sort(idxSharedWithSim),simScanIter)');
    end
end

% Average together all non nan values in matrix
avgBraggFrac = nanmean(avgBraggFrac,2)';
avgTotNum    = nanmean(avgTotNum,2)';
end