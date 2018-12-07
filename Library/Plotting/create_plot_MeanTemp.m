function create_plot_MeanTemp(analyVar,avgDataset)
% Function to plot the spectrum of atomic number.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%   avgDataset   - Cell of structures similar to indivDataset but
%                  containing only averaged values across similar scans
%                  (scans which share the same experimental parameters)
%
% OUTPUTS:
%   Creates plots showing the mean cloud radius

%% Loop through each batch file and image
%%%%%%%-----------------------------------%%%%%%%%%%
for uniqScanIter = 1:length(analyVar.uniqScanList);
    % Reference variables in structure by shorter names for convenience
    % (will not create copy in memory as long as the vectors are not modified)
    indVar       = analyVar.funcDataScale(avgDataset{uniqScanIter}.simScanIndVar);
    avgAtomTempX = avgDataset{uniqScanIter}.avgAtomTempX; stdAtomTempX = avgDataset{uniqScanIter}.stdAtomTempX;
    avgAtomTempY = avgDataset{uniqScanIter}.avgAtomTempY; stdAtomTempY = avgDataset{uniqScanIter}.stdAtomTempY;
    
%% Plot of radius in X & Y
%%%%%%%-----------------------------------%%%%%%%%%%
    figNum = analyVar.figNum.meanTemp; 
    tempLabel = {'Avg. Temperature [nK]', 'Avg. Temperature [nK]'}; 
    tempTitle = {'','Y-axis','X-axis'};
    default_plot(analyVar,[uniqScanIter length(analyVar.uniqScanList)],...
        figNum,tempLabel,tempTitle,analyVar.uniqScanList,...
         repmat(indVar,1,2)',[avgAtomTempY'; avgAtomTempX'],'Errorbar',[stdAtomTempY'; stdAtomTempX']);
    hFig = gcf;
    hFig.CurrentAxes = hFig.Children(end - 1); axis tight
    hFig.CurrentAxes = hFig.Children(end); axis tight
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
end
end