function rydDressOut = Rydberg_Dress(lcl_analyVar,lcl_indivDataset,lcl_avgDataset)
% Rydberg_Dress is a specialized function for analyzing data associated with
% the rydberg dressing and rydberg molecule experiments studied in 2014.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%

%% Decide to use averaged data or individual scans
avgAutoFlag = length(cell2mat(lcl_analyVar.posOccurUniqVar)) > length(lcl_analyVar.posOccurUniqVar) & lcl_analyVar.uniqScanList ~= 0;

%% Fit images and return parameters
% v2Struct will transfer the workspace of Spectrum_Fit into the current workspace
if lcl_analyVar.flagSpecFit
    v2struct(Spectrum_Fit(lcl_analyVar,lcl_indivDataset,lcl_avgDataset));
end

%% Fit and Plot DC Stark Shift data
% Used for Rydberg Molecule study - added 07/18/2014
if avgAutoFlag & lcl_analyVar.flagDC_Stark & lcl_analyVar.flagSpecFit
    DC_Stark_shift(lcl_analyVar,lineCenter)
elseif avgAutoFlag & lcl_analyVar.flagDC_Stark & ~lcl_analyVar.flagSpecFit
    error('Please enable spectrum fitting before continuing')
end

%% Fit and Plot AC Stark Shift data
% Used for EIT study - added 10/13/2014
if avgAutoFlag & lcl_analyVar.flagAC_Stark & lcl_analyVar.flagSpecFit
    AC_Stark_shift(v2struct)
elseif avgAutoFlag & lcl_analyVar.flagAC_Stark & ~lcl_analyVar.flagSpecFit
    error('Please enable spectrum fitting before continuing')
end

%% Plot Rydberg Dressing specific data
% Used for Rydberg Dressing study - added 07/18/2014
if lcl_analyVar.flagRydPlot
    Rydberg_Dress_Spectra(lcl_analyVar,fieldIndx,rabiFreq,lineCenter,fullWidth)
elseif avgAutoFlag & ~lcl_analyVar.flagDC_Stark & lcl_analyVar.flagRydPlot
    error('Please enable DC Stark shift fitting before continuing')
end

%% Bragg Spectroscopy Analysis
if lcl_analyVar.flagBraggSc
    braggOut = Bragg_Scat(lcl_analyVar,lcl_indivDataset,lcl_avgDataset);
end

%% Autler-Townes Analysis
if lcl_analyVar.flagAutTown
    autTownOut = AutlerTownes_Spectra(lcl_analyVar,lcl_indivDataset,lcl_avgDataset);
end

%% Plot data with Offset
if lcl_analyVar.flagPlotOff & avgAutoFlag
    plotWithOffset(lcl_analyVar,lcl_avgDataset);
end

%% Pack workspace into a structure for output
% If you don't want a variable output prefix it with lcl_
rydDressOut = who();
rydDressOut = v2struct(cat(1,'fieldNames',rydDressOut(cellfun('isempty',regexp(rydDressOut,'\<lcl_')))));
end