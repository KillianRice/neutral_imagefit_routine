function specOut = Spectrum_Fit(lcl_analyVar,lcl_indivDataset,lcl_avgDataset)
% Spectrum fit is a parent function to fit spectra and extract the rabi
% frequencies from single scan datasets as well as average datasets.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%   avgDataset   - Cell of structures containing grouped data for averaging
%
% OUTPUTS:
%   specOut      - Structure containing all variables in the workspace of Spectrum_Fit
%                  This is the same behavior as AnalysisVariables and can be used to
%                  facilitate calling Spectrum_Fit as a generalized routine.
%
% NOTE: Remember to update your x-axis conversion and time beneath variable
%       conversion

%% Variable conversion
% Replace powers and conversions for data input into fitting routine
% All calculations are done in SI then modified for plotting
scaleFuncs.ExpTime = @(expTime) expTime*1e-3; % data input in ms

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
% Also apply any scaling functions
indVar  = cellfun(@(x) x.(indField{fieldIndx}),curDataset,'UniformOutput',0);
partNum = cellfun(@(x) x.(numField{fieldIndx}),curDataset,'UniformOutput',0);
label   = lcl_analyVar.(labelField{fieldIndx});

%% Fit spectra to loss equation
% Use single body loss equation
[amplitude, lineCenter, fullWidth] = one_body_spectra(indVar,partNum,label)

%% Calculate Rabi Frequencies 
% Rabi Frequency - derivation in folder
if isempty(find(avgAutoFlag == 0, 1))
    simExpTimes = cellfun(@(x) x(1),lcl_analyVar.posOccurUniqVar);
    expTime     = scaleFuncs.ExpTime(lcl_analyVar.expTime(simExpTimes));
else
    expTime     = scaleFuncs.ExpTime(lcl_analyVar.expTime);
end
% rabi freq. in s^(-1) (not Hz)
rabiFreq      = (2*pi)^(-1/4).*sqrt(amplitude(:,1)./expTime);
rabiFreq(:,2) = (1/2).*abs(rabiFreq).*amplitude(:,2)./abs(amplitude(:,1))

%% Pack workspace into a structure for output
% If you don't want a variable output prefix it with lcl_
specOut = who();
specOut = v2struct(cat(1,'fieldNames',specOut(cellfun('isempty',regexp(specOut,'\<lcl_')))));
end