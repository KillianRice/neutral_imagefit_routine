function create_plot_MeanNum(analyVar,avgDataset)
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
%   Creates plots showing the mean atomic number. Valid plots are total number, 
%   BEC number, or lattice peak number
%

%% Loop through each batch file and image
%%%%%%%-----------------------------------%%%%%%%%%%
for uniqScanIter = 1:length(analyVar.uniqScanList);
    % Reference variables in structure by shorter names for convenience
    % (will not create copy in memory as long as the vectors are not modified)
    indVar    = avgDataset{uniqScanIter}.simScanIndVar;
    avgTotNum = avgDataset{uniqScanIter}.avgTotNum; stdTotNum = avgDataset{uniqScanIter}.stdTotNum;
    avgBECNum = avgDataset{uniqScanIter}.avgBECNum; stdBECNum = avgDataset{uniqScanIter}.stdBECNum;
    avgCondFrac = avgDataset{uniqScanIter}.avgCondFrac; stdCondFrac = avgDataset{uniqScanIter}.stdCondFrac;
    
%% Spectrum of Total number 
%%%%%%%-----------------------------------%%%%%%%%%%
% Decide whether to normalize the spectrum or not
if analyVar.normMeanNum
    numNormFactor = mean([avgTotNum(1:5)' avgTotNum(end-5:end)']);
else
    numNormFactor = 1;
end

    figNum = analyVar.figNum.meanNum; 
    numLabel = {'Average Total Number'}; 
    numTitle = {};
    default_plot(analyVar,[uniqScanIter length(analyVar.uniqScanList)],...
        figNum,numLabel,numTitle,analyVar.uniqScanList,...
        analyVar.funcDataScale(indVar)',avgTotNum'./numNormFactor,'Errorbar',stdTotNum'./numNormFactor);
    axis tight
    curYLim = get(gca,'YLim');
    ylim([0 curYLim(2)*1.2]);
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

%% BEC statistics for bimodal distributions    
    if strcmpi(analyVar.InitCase,'Bimodal')
    %% Spectrum of BEC Number
    %%%%%%%-----------------------------------%%%%%%%%%%
        figNum = analyVar.figNum.meanBEC; 
        numLabel = {'Average Condensate Number'}; 
        numTitle = {};
        default_plot(analyVar,[uniqScanIter length(analyVar.uniqScanList)],...
            figNum,numLabel,numTitle,analyVar.uniqScanList,...
            analyVar.funcDataScale(indVar)',avgBECNum','Errorbar',stdBECNum');
        curYLim = get(gca,'YLim');
        ylim([0 curYLim(2)]);

    %% Plot of condensate fraction
    %%%%%%%-----------------------------------%%%%%%%%%%
        figNum = analyVar.figNum.meanFrac; 
        numLabel = {'Average Condensate Fraction [%]'}; 
        numTitle = {};
        default_plot(analyVar,[uniqScanIter length(analyVar.uniqScanList)],...
            figNum,numLabel,numTitle,analyVar.uniqScanList,...
            analyVar.funcDataScale(indVar)',avgCondFrac','Errorbar',stdCondFrac');
        ylim([0 100]);
    end
end
end