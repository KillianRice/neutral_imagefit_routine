function KDOut = KapitzaDirac(analyVar,indivDataset,avgDataset)
% Function to analyze Kapitza-Dirac scattering of the lattice. This process
% is used when calibrating the lattice depth.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:

%% Preallocate and initialize variables
skipStart  = 1; % Skip this many points from beginning of data to fit freq.
skipEnd    = 0; % Skip this many points from end of data to fit freq.

%% Find axis being analyzed
latAx = analyVar.LatticeAxesFit{1}(2:2:end); % Find axes being analyzed

%% Check if Kapitza-Dirac can be applied
% This analysis requires a lattice sampleType and windows outside of the central window to be defined
if not(strcmpi(analyVar.sampleType,'Lattice') && sum(analyVar.LatticeAxesFit{1}) > 1)
    if ~strcmpi(analyVar.sampleType,'Lattice')
        errType = 'sampleType = Lattice';
    elseif ~sum(analyVar.LatticeAxesFit) > 0
        errType = 'specifying Lattice axis window';
    end
    error('KD_Analysis:NotApplicable',['Kapitza-Dirac scattering analysis enabled without %s. Please check ',...
            'AnalysisVariables and run again.'],errType)
end
% Should also only apply Kapitza-Dirac in 1D lattice
if sum(latAx) > 1
    error('KD_Analysis:NotApplicable',['Kapitza-Dirac scattering analysis only applicable with 1-D lattice. ',...
        'Please check AnalysisVariables and run again.'])
end

%% Plot spectrum 
% Call the function that plots spectrum of +/- 2hk peaks
plot_num_2hk(analyVar,indivDataset)

%% Find average number
% Average number is used for fitting spectrum
% Loop through all unique values used for averaging
for uniqScanIter = 1:length(analyVar.uniqScanList);
    %% Preallocate nested loop variables
    [avg2hkPeakNum, avg0hkPeakNum]...
        = deal(NaN(length(avgDataset{uniqScanIter}.simScanIndVar),...
                length(analyVar.posOccurUniqVar{uniqScanIter})));
    clear smoothAvg
  
    %% Loop through all scans that share the current value of the averaging variable
    for simScanIter = 1:length(analyVar.posOccurUniqVar{uniqScanIter})
        % Assign the current scan to be opened
        basenameNum = analyVar.posOccurUniqVar{uniqScanIter}(simScanIter);
        
        % Use local variables for window atom numbers for convenience
        winTotNum  = indivDataset{basenameNum}.winTotNum;
        
        % Convert numbers to integer to remove any floating point residual
        % differences prior to comparison.
        tempInt_simScanIndVar = int32(avgDataset{uniqScanIter}.simScanIndVar.*analyVar.compPrec);
        tempInt_imagevcoAtom = int32(indivDataset{basenameNum}.imagevcoAtom.*analyVar.compPrec);
        
        
        % Find the intersection of the scanned variable with the list of all possible values
        [~,idxSharedWithSim,idxSharedInBatch] = intersect(tempInt_simScanIndVar,tempInt_imagevcoAtom);
        
%         Convert numbers back to doubles so they can be used with regular, e.g. exp(), functions.
%         avgDataset{uniqScanIter}.simScanIndVar = double(tempInt_simScanIndVar./analyVar.compPrec);
%         indivDataset{basenameNum}.imagevcoAtom = double(tempInt_imagevcoAtom./analyVar.compPrec);
        
               
        % Find number in the peaks of interest
        avg2hkPeakNum(idxSharedWithSim,simScanIter) = sum(winTotNum(2:3,idxSharedInBatch))./sum(winTotNum(:,idxSharedInBatch));
        avg0hkPeakNum(idxSharedWithSim,simScanIter) = winTotNum(1,idxSharedInBatch)./sum(winTotNum(:,idxSharedInBatch));
    end
    
    % Use local variables for similar independent variable for convenience
    indVar = analyVar.funcDataScale(avgDataset{uniqScanIter}.simScanIndVar);
    
    % Average together all non nan values in matrix
    avg2hkPeakNum = nanmean(avg2hkPeakNum,2);
    avg0hkPeakNum = nanmean(avg0hkPeakNum,2);
    
    % Smooth data for fitting
    smoothAvg(1,:) = smooth(avg0hkPeakNum,'sgolay');
    smoothAvg(2,:) = smooth(avg2hkPeakNum,'sgolay');
    
    
    figure; hold on; grid on
    [Harmfit, latticedepth] = HarmFitAndLatticedepth(indVar(skipStart:end-skipEnd)'*1e-3,...
        smoothAvg(2,skipStart:end-skipEnd),analyVar);
    
    % Parameters for population oscillation plot
    arm = analyVar.winToFit{2};
    exp = 0;
    
    plot_PopOscillation(indVar*1e3,smoothAvg,{arm,analyVar.uniqScanList,exp,latticedepth})
    
% Pack workspace into a structure for output
% If you don't want a variable output prefix it with lcl_
KDOut = who();
KDOut = v2struct(cat(1,'fieldNames',KDOut(cellfun('isempty',regexp(KDOut,'\<lcl_')))));
end
    