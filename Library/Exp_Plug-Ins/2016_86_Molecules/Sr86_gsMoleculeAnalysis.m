function sr86MoleOut = Sr86_gsMoleculeAnalysis(lcl_analyVar,lcl_indivDataset,lcl_avgDataset)
% Plugin for analyzing ground state molecule data using strontium 86.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%   moleOut      - Structure containing all variables in the workspace of Sr86_gsMoleculeAnalysis
%                  This is the same behavior as AnalysisVariables and can be used to
%                  facilitate calling this function as a generalized routine.
%

%% Variable conversion
% Replace powers and conversions for data input into fitting routine
% All calculations are done in SI then modified for plotting

%% Decide to use averaged data or individual scans
avgAutoFlag = length(cell2mat(lcl_analyVar.posOccurUniqVar)) > length(lcl_analyVar.posOccurUniqVar) & lcl_analyVar.uniqScanList ~= 0;

%% Extract independent variables and number from correct structure
indField   = {'imagevcoAtom'   'simScanIndVar'};  % fields of independent variables in indivDataset and avgDataset
numField   = {'winTotNum'      'avgTotNum'};      % number data in indivDataset and avgDataset
labelField = {'timevectorAtom' 'uniqScanList'};   % label used for each plot

% Selected datatype picks the data to be analyzed
if avgAutoFlag
    curDataset = lcl_avgDataset;   fieldIndx = 2;
else
    curDataset = lcl_indivDataset; fieldIndx = 1;
end

% Save variables from cell of structures
% Also apply any scaling functions
lcl_indVar  = cellfun(@(x) x.(indField{fieldIndx}),curDataset,'UniformOutput',0);
lcl_partNum = cellfun(@(x) x.(numField{fieldIndx}),curDataset,'UniformOutput',0);
lcl_label   = lcl_analyVar.(labelField{fieldIndx});

%% Load intensity data from picoscope
if lcl_analyVar.flagGSMoleInten
    intenOut = picoIntenAnalysis(lcl_analyVar,lcl_indivDataset);
end

%% Pack workspace into a structure for output
% If you don't want a variable output prefix it with lcl_
sr86MoleOut = who();
sr86MoleOut = v2struct(cat(1,'fieldNames',sr86MoleOut(cellfun('isempty',regexp(sr86MoleOut,'\<lcl_')))));
end