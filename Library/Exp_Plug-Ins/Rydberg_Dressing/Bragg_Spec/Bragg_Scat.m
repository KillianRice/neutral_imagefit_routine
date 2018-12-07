function braggOut = Bragg_Scat(lcl_analyVar,lcl_indivDataset,lcl_avgDataset)
%Function to calculate Bragg spectra values
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%   braggOut     - structure containing the workspace of Bragg_Scat compiled at the end
%                  of the function (same behavior as AnalysisVariables)

%% Initialize variables
figSpec    = figure; % Composite figure of all peaks
specFitFig = figure; % Fubplot figure showing fits

%% Calculate Bragg Peak Populations
if lcl_analyVar.winID == 0
    error('Bragg peak populations cannot be calculated since no window types are specified.')
else
    [indVar,braggFrac,totPop,expTime,legData,numPeak] = getBraggPop(lcl_analyVar,lcl_indivDataset,lcl_avgDataset);
end

% Initialize loop variables
coeffOut = cellfun(@(x) cell(size(x)), braggFrac,'UniformOutput',0);
for iterVar = 1:length(indVar)
    for i = 1:numPeak(iterVar)
        %% Plotting
        % Define plot label and title
        axLabel = {{'Number in' 'Bragg Peak'}, {'Sum of all Peaks' '(Total Num)'}};
        numTitle = {'Bragg Populations'};
        % Plot call is not general and will only work with both peaks analyzed at the same time
        if i == numPeak(iterVar); maxIter = iterVar; else maxIter = 0; end
        % Loop through peaks
        [plotHan, axHan] = default_plot(lcl_analyVar,[iterVar maxIter],...
            figSpec,axLabel,numTitle,legData,...
            repmat(indVar{iterVar},1,2)',[braggFrac{iterVar}{i}; totPop{iterVar}{i}],...
            'Hide Legend',2);
        set(figSpec,'CurrentAxes',axHan(1)); axis tight
        
        %% Fit number spectrum
        if expTime % only fit if enabled for all scans
            coeffOut{iterVar}{i} = fit_bragg_spectra(lcl_analyVar,indVar{iterVar},braggFrac{iterVar}{i},...
                legData(iterVar),expTime(iterVar),figSpec,specFitFig,iterVar,length(indVar));
        else
            if ishandle(specFitFig)
                close(specFitFig)
            end
        end
    end
end

%% Pack workspace into a structure for output
% If you don't want a variable output prefix it with lcl_
braggOut = who();
braggOut = v2struct(cat(1,'fieldNames',braggOut(cellfun('isempty',regexp(braggOut,'\<lcl_')))));
end