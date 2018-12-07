function [braggFrac,totPop]  = largeAngleBragg(analyVar,indivDataset,avgDataset)
% Extracts the populations from spatially resolved fitting windows
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit routines
%   indivDataset - Cell of structures containing all scan/batch specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%   braggPop - Cell of bragg populations in each scan
%   totPop   - Cell of the total population in each scan

%% Decide to use averaged data or individual scans
avgAutoFlag = length(cell2mat(analyVar.posOccurUniqVar)) > length(analyVar.posOccurUniqVar) & analyVar.uniqScanList ~= 0;

%% Extract number and normalize to total number
if avgAutoFlag
    % Initialize loop variables
    [braggFrac, totPop] = deal(cell(length(analyVar.uniqScanList)));
    for uniqScanIter = 1:length(analyVar.uniqScanList);
        % Build temporary cell with Bragg number of a certain peak for each member of the group
        tmpBragg = cellfun(@(x) x.winTotNum(2,:), indivDataset(analyVar.posOccurUniqVar{uniqScanIter}),'UniformOutput',0);
        tmpTot   = cellfun(@(x) sum(x.winTotNum,1), indivDataset(analyVar.posOccurUniqVar{uniqScanIter}),'UniformOutput',0);
        % Average Bragg number and total number
        [braggFrac{uniqScanIter} totPop{uniqScanIter}] = avgBraggNum(analyVar,indivDataset,avgDataset,uniqScanIter,tmpBragg,tmpTot);
    end
else
    braggFrac = cellfun(@(x) x.winTotNum(2,:)./sum(x.winTotNum,1), indivDataset,'UniformOutput',0);
    totPop    = cellfun(@(x) sum(x.winTotNum,1), indivDataset,'UniformOutput',0);
end
end