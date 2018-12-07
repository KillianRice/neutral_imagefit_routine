function scopeTrace = read_TDS2024(file)
% Function to read scope traces from TDS2024 scope
%
% INPUTS:
%   file     - full file path to open the trace file
%
% OUTPUTS:
%   scopeTrace - structure containing the following fields:
%       time: sub-structure with fields zero, inc, and axis. All units in seconds.
%       ch# : sub-structure with fields offset, scale, and value for each
%               channel enabled on the scope. Offset unit is arbitrary, scale converts to Volts
%
% NOTE:
%   1. This file expects output as written from Labview by tds2024_4channel_mod.vi
%   2. Since each channel can have a different scale, the scope returns an 8-bit signed integer
%       that represents the measured voltage. This integer varies from -128 -> 127 and is converted 
%       to real voltage by:  Volts = (scopeVal + Offset)*scale. The file tds2024_4channel_mod does
%       this before saving.

%% Initialize variables
chanSpec = 1:4; %Channel specifier
numParam = 4;   % Number of parameters return for each channel

%% Read file
% Expects 1 header line
try
    scopeRead = importdata(file,'\t',1);
catch
    error(['Error importing file ' file]) 
end

% Find which channels contain data
chanOn = chanSpec(~isnan(scopeRead.data(:,1)));

% Scan header and save channel parameters
headCell = textscan(scopeRead.textdata{:}, ['%s' repmat('%f',1,numParam*length(chanSpec))],...
    'CollectOutput',1);
headData = reshape(headCell{2},[length(chanSpec),numParam])';

%% Construct Time axis
scopeTrace.time.zero = unique(nonzeros(headData(~isnan(headData(:,1)),1))); 
scopeTrace.time.inc  = unique(headData(~isnan(headData(:,2)),2));
scopeTrace.time.axis = (0:length(scopeRead.data(chanOn(1),:))-1)*scopeTrace.time.inc + scopeTrace.time.zero;

%% Construct channel data
for i = chanOn
    scopeTrace.( ['ch' num2str(i)]).offset = headData(i,3);
    scopeTrace.( ['ch' num2str(i)]).scale  = headData(i,4);
    scopeTrace.( ['ch' num2str(i)]).value  = scopeRead.data(i,:);
end
end