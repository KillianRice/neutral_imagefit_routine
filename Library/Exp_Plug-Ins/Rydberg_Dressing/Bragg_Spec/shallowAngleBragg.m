function [braggFrac,totPop]  = shallowAngleBragg(analyVar,indivDataset,avgDataset)
% Function to find the population of Bragg peaks that are not clearly spatially resolved.
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit routines
%   indivDataset - Cell of structures containing all scan/batch specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%   braggPop - Cell of bragg populations in each scan
%   totPop   - Cell of the total population in each scan

%% Error if more than central window is defined
if length(analyVar.winToFit) > 1
    error('Disable other windows except central window')
end

%% Decide if data should be averaged
avgAutoFlag = length(cell2mat(analyVar.posOccurUniqVar)) > length(analyVar.posOccurUniqVar) & analyVar.uniqScanList ~= 0;

%% Initialize Variables
[braggNum] = deal(cell(analyVar.numBasenamesAtom,1));

%% Loop through each scan and get number
for basenameNum = 1:analyVar.numBasenamesAtom
    % Determine which peaks need to be fit
    if analyVar.winID(basenameNum) == 2
        numPeak = [-1 1];
    else
        numPeak = 1;
    end
    
    % Numerical integration
    [braggNum{basenameNum},iPkSub,jPkSub] = numIntBraggPop(analyVar,indivDataset,basenameNum,numPeak);
    
    % Show numerical integration region if fitEval enabled
    if analyVar.plotFitEval
        drawIntReg(analyVar,indivDataset,basenameNum,iPkSub,jPkSub,numPeak)
    end
end

%% Calculate the total number of atoms
totNum = cellfun(@(x,y) sum(cell2mat(x),1) + y.winTotNum,braggNum,indivDataset,'UniformOutput',0);

%% Average data if needed
if avgAutoFlag
    % Initialize loop variables
    [braggFrac, totPop] = deal(cellfun(@(x) cell(length(numPeak),1), cell(length(analyVar.uniqScanList),1),'UniformOutput',0));
    for uniqScanIter = 1:length(analyVar.uniqScanList);
        % Check that the window number is the same for all scans in a group
        if range(analyVar.winID(analyVar.posOccurUniqVar{uniqScanIter})) ~= 0
            error('winID is not the same for each member of the group labeled %g',analyVar.uniqScanList(uniqScanIter))
        end
        
        % Loop through Bragg peaks
        for i = 1:length(numPeak)
            % Build temporary cell with Bragg number of a certain peak for each member of the group
            tmpBragg = cellfun(@(x) x{i}, braggNum(analyVar.posOccurUniqVar{uniqScanIter}),'UniformOutput',0);
            tmpTot   = totNum(analyVar.posOccurUniqVar{uniqScanIter});
            % Average Bragg number and total number
            [braggFrac{uniqScanIter}{i} totPop{uniqScanIter}{i}] = avgBraggNum(analyVar,indivDataset,avgDataset,uniqScanIter,tmpBragg,tmpTot);
        end
    end
else
    % Make both braggFrac and totPop the same size for iterating through and plotting
    braggFrac = braggNum;
    totPop    = cellfun(@(x,y) repmat({x},size(y)), totNum, braggNum,'UniformOutput',0);
    % Replace the Bragg peak number with fraction
    for i = 1:length(braggNum)
        for j = 1:length(braggNum{i})
            braggFrac{i}{j} = braggNum{i}{j}./totNum{i};
        end
    end
end
end