function TripletFitOut = TripletP1Fit(lcl_analyVar,lcl_indivDataset,lcl_avgDataset)
% Function to fit the 3P1 m_j = {-1,0,1} Zeeman levels for use in zeroing the magnetic field.
%
% Master Batch file column 4 acting as a flag for the munber of peaks the fit function should expect.
% 
% INPUTS: 
% 
% OUTPUTS: 
%
% Josh Hill Dec. 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTANTS
%%-----------------------------------------------------------------------%%
% e     = 1.6022e-19; % In Coulombs
% hbar  = 1.0546e-34; % In J*s
% h     = hbar*2*pi;
% me    = 9.1094e-31; % In kg
% muB   = e*hbar/(2*me); %Bohr magneton in J/T


% PARAMETERS
%%-----------------------------------------------------------------------%%
% S      = 1;
% L      = 1;
% J      = 1;
% mj     = 1;

% % Currents
% Ix     =
% Iy     =
% Iz     =
% 
% 
% B      = sqrt((Bx-Bx0)^2+(By-By0)^2+(Bz-Bz0)^2);
% 
% 
% gj     = 3/2+(S*(S+1)-L*(L+1))/(2*J*(J+1));
% muz    = -mj*gj*muB; % Z component of dipole
% deltaf = mj*gj*muB*B/h % Frequency Splitting in Hz



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
% Input x, applies the function x.field to every cell in the array curDataseet. I.e. it constructs curDataset(x).indField{fieldIndx} for all x elements.
indVar  = cellfun(@(x) x.(indField{fieldIndx}),curDataset,'UniformOutput',0);  
partNum = cellfun(@(x) x.(numField{fieldIndx}),curDataset,'UniformOutput',0);
label   = lcl_analyVar.(labelField{fieldIndx});

%% Fit spectra to loss equation
numPeaks                = lcl_analyVar.synthGrossFreq(1); % Master Batch file column 4 acting as a flag for the munber of peaks the fit function should expect. Listed in ~line 241 in AnalysisVariables
% Use single body loss equation
[lineCenter, fullWidth] = one_body_spectra3P1Summed(indVar,partNum,label,numPeaks)


%% Pack workspace into a structure for output
% If you don't want a variable output, prefix it with lcl_
TripletFitOut = who();
TripletFitOut = v2struct(cat(1,'fieldNames',TripletFitOut(cellfun('isempty',regexp(TripletFitOut,'\<lcl_')))));
end