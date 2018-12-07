function intenOut = picoIntenAnalysis(lcl_analyVar,lcl_indivDataset)
% Function to load intensity fluctuation data and create plots for anaylzing 
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit routines
%   indivDataset - Cell of structures containing all scan/batch specific data
%
% OUTPUTS:
%   intenOut    - Structure containing all variables in the workspace of Spectrum_Fit
%                 This is the same behavior as AnalysisVariables and can be used to
%                 facilitate calling Spectrum_Fit as a generalized routine.
%

%% Initialize variables
binWidthFunc = @(vec) round(std(vec)/8,3);
figVarNum    = figure;
figHistNum   = figure;
    
%% Load intensity fluctuation data
intenOut = loadIntenData(lcl_analyVar,lcl_indivDataset);

%% Loop through each batch file and image
%%%%%%%-----------------------------------%%%%%%%%%%

for basenameNum = 1:lcl_analyVar.numBasenamesAtom
    % Reference variables in structure by shorter names for convenience
    % (will not create copy in memory as long as the vectors are not modified)
    indVar    = lcl_indivDataset{basenameNum}.imagevcoAtom;
    intenVar  = intenOut{basenameNum}.p2p;
    
%% Plot Intensity variation data
    figYLabel = {'Photodiode Voltage [V]'};
    figTitle = {};
    default_plot(lcl_analyVar,[basenameNum lcl_analyVar.numBasenamesAtom],...
        figVarNum,figYLabel,figTitle,lcl_analyVar.timevectorAtom,...
        lcl_analyVar.funcDataScale(indVar)',intenVar');
    
%% Plot Histogram of fluctuations
    figure(figHistNum);
    histogram(intenVar,'BinWidth',binWidthFunc(intenVar))
    ylabel(sprintf('Occurences (Bin Width = %g V',binWidthFunc(intenVar)),...
        'FontSize'  ,   lcl_analyVar.labelFontSize  ,...
        'FontWeight',   'bold'                  );
    xlabel(lcl_analyVar.xDataLabel,...
        'FontSize'  ,   lcl_analyVar.labelFontSize  ,...
        'FontWeight',   'bold'                  );
    set(gca,...
        'LineWidth' ,   1                       ,...
        'FontSize'  ,   lcl_analyVar.axisFontSize   ,...
        'FontWeight',   'bold'                  ,...
        'Box'       ,   'on'                    );
    hold on

end
end