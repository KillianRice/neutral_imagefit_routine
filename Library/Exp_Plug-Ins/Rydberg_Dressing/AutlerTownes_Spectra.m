function autTownOut = AutlerTownes_Spectra(lcl_analyVar,lcl_indivDataset,lcl_avgDataset)
%AUTLERTOWNES_SPECTRA - Fits double peaked Autler-Townes data
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%   autTownOut      - Structure containing all variables in the workspace of Spectrum_Fit
%                     This is the same behavior as AnalysisVariables and can be used to
%                     facilitate calling as a generalized routine.

%% Initialize Variables
autTownFitFig = figure;
cntrPosFig    = figure;
ampFig        = figure;
fwhmFig       = figure;

%% Define Physical Functions Used in Calculations
% coeffs has elements coeffs = [offset, amp_1, cntr_1, sigma_1, amp_2, cntr_2, sigma_2]
singPeakFit  = @(coeffs,x) coeffs(1)*(1 - abs(coeffs(2))*exp(-(x-coeffs(3)).^2/(2*coeffs(4)^2)));
autTownFit   = @(coeffs,x) coeffs(1)*(1 - abs(coeffs(2))*exp(-(x-coeffs(3)).^2/(2*coeffs(4)^2)) - coeffs(5)*exp(-(x-coeffs(6)).^2/(2*coeffs(7)^2)));
hyperFitFunc = @(coeffs,x) sqrt((x-coeffs(1)).^2 + coeffs(2)^2); % Fit of the energy difference between dressed states

% Full Width Half Max
fwhmFunc = @(sigma) abs(sigma*2*sqrt(2*log(2)));

%% Decide to use averaged data or individual scans
avgAutoFlag = length(cell2mat(lcl_analyVar.posOccurUniqVar)) > length(lcl_analyVar.posOccurUniqVar) & lcl_analyVar.uniqScanList ~= 0;

%% Extract independent variables and number from correct structure
indField   = {'imagevcoAtom' 'simScanIndVar'};  % fields of independent variables in indivDataset and avgDataset
numField   = {'winTotNum' 'avgTotNum'};         % number data in indivDataset and avgDataset
labelField = {'timevectorAtom' 'uniqScanList'}; % label used for each plot

% Selected datatype picks the data to be analyzed
if avgAutoFlag
    curDataset = lcl_avgDataset;   fieldIndx = 2;
else
    curDataset = lcl_indivDataset; fieldIndx = 1;
end

% Save variables from cell of structures
indVarCell  = cellfun(@(x) x.(indField{fieldIndx}),curDataset,'UniformOutput',0);
partNumCell = cellfun(@(x) x.(numField{fieldIndx}),curDataset,'UniformOutput',0);
labelVec    = lcl_analyVar.(labelField{fieldIndx});

%% Initialize loop variables
[autTownOut.amp_left, autTownOut.amp_right, autTownOut.cntr_left,...
    autTownOut.cntr_right, autTownOut.fwhm_left, autTownOut.fwhm_right] = deal(zeros(length(indVarCell),1));

autTownFitCell = cell(1,length(indVarCell));

%% Fitting Loop
for iterVar = 1:length(indVarCell)
    % Reference variables in structure by shorter names for convenience
    indVar = nan(1,length(indVarCell{iterVar}));  indVar(:) = lcl_analyVar.funcDataScale(indVarCell{iterVar});
    totNum = nan(1,length(partNumCell{iterVar})); totNum(:) = partNumCell{iterVar};

    % Initial Guesses
    smoothNum            = smooth(totNum,'sgolay',1);  % smooth for helping make guesses
    [minVal_1, minLoc_1] = min(smoothNum(indVar < 0)); % find location of minimum 1
    [minVal_2, minLoc_2] = min(smoothNum(indVar > 0)); % find location of minimum 2
    
    initOffset   = mean(totNum([1:5 end-5 end]));   % Mean of first 5 and last 5 points
    initAmp_1    = (1 - minVal_1/initOffset);       % Value below guessed offset
    initAmp_2    = (1 - minVal_2/initOffset);       % Value below guessed offset
    initCenter_1 = indVar(length(indVar(indVar > 0)) + minLoc_1);                % Location of minimum in smoothed data
    initCenter_2 = indVar(minLoc_2);
    initWidth_1  = sqrt(abs(sum((indVar(indVar < 0) - initCenter_1).^2.*totNum(indVar < 0))./sum(totNum(indVar < 0))))/3;
    initWidth_2  = sqrt(abs(sum((indVar(indVar > 0) - initCenter_2).^2.*totNum(indVar > 0))./sum(totNum(indVar > 0))))/3;
    
    % Fit peaks indivually
    leftPeakFit  =  NonLinearModel.fit(indVar',totNum,singPeakFit,[initOffset initAmp_1 initCenter_1 initWidth_1],...
        'CoefficientNames',{'Offset', 'Amp (Left)','Center (Left)','Sigma (Left)'});
    rightPeakFit  =  NonLinearModel.fit(indVar',totNum,singPeakFit,[initOffset initAmp_2 initCenter_2 initWidth_2],...
        'CoefficientNames',{'Offset', 'Amp (Right)','Center (Right)','Sigma (Right)'});
    
    % Use the individual fits as guesses for complete fit
    meanOffset = mean([leftPeakFit.Coefficients.Estimate(1) rightPeakFit.Coefficients.Estimate(1)]);
    beta0      = [meanOffset leftPeakFit.Coefficients.Estimate(2:end)' rightPeakFit.Coefficients.Estimate(2:end)'];

    % Composite function fiitting routine
    autTownFitCell{iterVar} = NonLinearModel.fit(indVar',totNum,autTownFit,beta0,...
        'CoefficientNames',{'Offset', 'Amp (Left)','Center (Left)','Sigma (Left)'...
        ,'Amp (Right)','Center (Right)','Sigma (Right)'});

    % Calculate output quantities
    autTownOut.amp_left(iterVar)     = autTownFitCell{iterVar}.Coefficients.Estimate(2);
    autTownOut.cntr_left(iterVar)    = autTownFitCell{iterVar}.Coefficients.Estimate(3);
    autTownOut.fwhm_left(iterVar)    = fwhmFunc(autTownFitCell{iterVar}.Coefficients.Estimate(4));
    autTownOut.amp_right(iterVar)    = autTownFitCell{iterVar}.Coefficients.Estimate(5);
    autTownOut.cntr_right(iterVar)   = autTownFitCell{iterVar}.Coefficients.Estimate(6);
    autTownOut.fwhm_right(iterVar)   = fwhmFunc(autTownFitCell{iterVar}.Coefficients.Estimate(7));
end

% Try to fit the change in energy difference between the dressed states
    if avgAutoFlag & length(lcl_analyVar.uniqScanList) > 2
        % Data for fitting
        enDiff     = autTownOut.cntr_right - autTownOut.cntr_left;
        uvDetuning = (lcl_analyVar.uniqScanList - 76.19074)*-8;
        
        % Initial guesses
        [initRabi, minLoc] = min(enDiff);
        initCntr           = uvDetuning(minLoc);
        
        % Fitting routine
        hyperFit = NonLinearModel.fit(uvDetuning,enDiff,hyperFitFunc,[initCntr initRabi],'CoefficientNames',{'Center' 'Rabi Freq.'});
        
        % Save output quantities
        autTownOut.hyperCntr = hyperFit.Coefficients.Estimate(1);
        autTownOut.hyperRabi = hyperFit.Coefficients.Estimate(2);
    end

%% Plotting Routine
[subPlotRows,subPlotCols] = optiSubPlotNum(length(indVarCell));
for iterVar = 1:length(indVarCell)
    % Reference variables in structure by shorter names for convenience
    indVar    = nan(1,length(indVarCell{iterVar}));  indVar(:) = lcl_analyVar.funcDataScale(indVarCell{iterVar});
    totNum    = nan(1,length(partNumCell{iterVar})); totNum(:) = partNumCell{iterVar};
    rangeInd  = minmax(lcl_analyVar.funcDataScale(cell2mat(indVarCell(:)))');
    fitIndVar = linspace(rangeInd(1),rangeInd(2),1e3)';
    
    % Plot number vs. fit for inspection
    figure(autTownFitFig); 
    subplot(subPlotRows,subPlotCols,iterVar); 
    dataHan   = plot(indVar,totNum,fitIndVar,autTownFitCell{iterVar}.predict(fitIndVar));
    title(num2str(labelVec(iterVar)),...
        'FontWeight'    ,   'bold'  ,...
        'FontSize'      ,   11      );
    xlabel('Detuning [MHz]',...
        'FontWeight'    ,   'bold'  ,...
        'FontSize'      ,   15      ); 
    grid on; 
    axis tight; 
    curYLim = get(gca,'YLim'); 
    set(gca,'YLim',[0 max(curYLim)])
    set(dataHan(1)        ,...
        'LineStyle'       ,   'none'  ,...
        'Marker'          ,   'o'     ,...
        'MarkerSize'      ,   7      ,...
        'MarkerFaceColor' ,   'b'     );
    set(dataHan(2),...
        'LineWidth' ,   2.5   )
    set(gca,...
        'FontSize'  ,   15  ,...
        'LineWidth' ,   2   )
    if iterVar == length(indVarCell);
        set(gcf,...
            'Name'  ,   'Autler-Townes Fits'  );
    end
end

if avgAutoFlag & length(lcl_analyVar.uniqScanList) ~= 1
% Plot the line centers vs. detuning
    figure(cntrPosFig);
    dataHan = plot(repmat(lcl_analyVar.uniqScanList*1e3,[1,2]), [autTownOut.cntr_left autTownOut.cntr_right]);
    grid on
    xlabel('UV Detuning [kHz]',...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    ylabel('^3P_1 Line Center [MHz]',...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    set(dataHan           ,...
        'LineStyle'       ,   'none'  ,...
        'Color'           ,   'k'     ,...
        'Marker'          ,   'o'     ,...
        'MarkerSize'      ,   15      );
    set(dataHan(1)        ,...
        'MarkerFaceColor' ,   'r'     );
    set(dataHan(2)        ,...
        'MarkerFaceColor' ,   'b'     );
    set(gca,...
        'FontSize'  ,   30      ,...
        'FontWeight',   'bold'  ,...
        'LineWidth' ,   2       )
    legend({'Left Peak' 'Right Peak'},...
        'Location'  ,   'NorthWest' )
    set(gcf,...
        'Name'  ,   'Line Center vs. UV Detuning'  );
    
% Plot the amplltudes vs. detuning
    figure(ampFig);
    dataHan = plot(repmat(lcl_analyVar.uniqScanList*1e3,[1,2]), [autTownOut.amp_left autTownOut.amp_right]);
    grid on
    xlabel('UV Detuning [kHz]',...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    ylabel('Peak Amplitude',...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    set(dataHan           ,...
        'LineStyle'       ,   'none'  ,...
        'Color'           ,   'k'     ,...
        'Marker'          ,   'o'     ,...
        'MarkerSize'      ,   15      );
    set(dataHan(1)        ,...
        'MarkerFaceColor' ,   'r'     );
    set(dataHan(2)        ,...
        'MarkerFaceColor' ,   'b'     );
    set(gca,...
        'FontSize'  ,   30      ,...
        'FontWeight',   'bold'  ,...
        'LineWidth' ,   2       )
    legend({'Left Peak' 'Right Peak'})
    set(gcf,...
        'Name'  ,   'Amplitude vs. UV Detuning'  );
    
% Plot the FWHM vs. detuning
    figure(fwhmFig);
    dataHan = plot(repmat(lcl_analyVar.uniqScanList*1e3,[1,2]), [autTownOut.fwhm_left autTownOut.fwhm_right]);
    grid on
    xlabel('UV Detuning [kHz]',...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    ylabel('^3P_1 FWHM [MHz]',...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    set(dataHan           ,...
        'LineStyle'       ,   'none'  ,...
        'Color'           ,   'k'     ,...
        'Marker'          ,   'o'     ,...
        'MarkerSize'      ,   15      );
    set(dataHan(1)        ,...
        'MarkerFaceColor' ,   'r'     );
    set(dataHan(2)        ,...
        'MarkerFaceColor' ,   'b'     );
    set(gca,...
        'FontSize'  ,   30      ,...
        'FontWeight',   'bold'  ,...
        'LineWidth' ,   2       )
    legend({'Left Peak' 'Right Peak'})
    set(gcf,...
        'Name'  ,   'FWHM vs. UV Detuning'  );
else
    close(cntrPosFig)
    close(ampFig)
    close(fwhmFig)
end

if avgAutoFlag & length(lcl_analyVar.uniqScanList) > 2
% Plot the difference in the energy of the dressed states
    figure;
    theoInd = linspace(min(uvDetuning),max(uvDetuning),1e3);
    dataHan = plot(uvDetuning, enDiff,theoInd,hyperFit.predict(theoInd'));
    grid on
    xlabel('UV Detuning [MHz]',...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    ylabel({'Energy difference between' 'dressed states [MHz]'},...
        'FontSize'  ,   35      ,...
        'FontWeight',   'bold'  );
    set(dataHan(1)                               ,...
        'LineStyle'       ,   'none'             ,...
        'Color'           ,   'k'                ,...
        'Marker'          ,   'o'                ,...
        'MarkerSize'      ,   15                 ,...
        'MarkerFaceColor' ,   [216 41 0]/255     );
    set(dataHan(2)                                  ,...
        'LineStyle'       ,   '--'                  ,...
        'LineWidth'       ,   5                     ,...
        'Color'           ,   [128 128 128]/255     );
    set(gca,...
        'FontSize'  ,   30      ,...
        'FontWeight',   'bold'  ,...
        'LineWidth' ,   2       )
    legHan = legend({'Data' 'Fit: $\Delta E = \hbar \sqrt{(\delta_2)^2+(\Omega_2)^2}$'});
    set(legHan, 'interpreter', 'latex')
    set(gcf,...
        'Name'  ,   'Energy Difference Fit'  );
end
end