function numBragg = intenNormNum(analyVar,indivDataset,basenameNum,numBragg)
% Function to read, manipulate, and average the intensity
% traces from the scope and use them to normalize
% number fluctuations of the Bragg spectroscopy data

% INPUTS:
%   analyVar     - Structure of all pertinent variables for the imagefit routines
%   indivDataset - Cell of structures containing all scan/batch specific data
%   avgDataset   - Cell of structures containing grouped and averaged data 
%   numBragg     - Number in the excited Bragg peak to be scaled
%
% OUTPUTS:
%   numBragg - Scaled number in the excited Bragg peak
%
% NOTE:
%   1. ch2 is dedicated to 
%   2. ch3 is dedicated to

% To Do
% 3. save scaling vector - possibly use try, catch block

%% Initialize Variables
PD_1_calib = 1; % TBD calibration factor [Intensity/volt]
PD_2_calib = 1; % TBD calibration factor [Intensity/volt]
fileSuffix = 'Traces.4scope'; % Filename suffix appended to scope traces

expTime = analyVar.expTime(basenameNum)*1e-3; % Get exposure time 

%% Define formulas
omega = @(inten) analyVar.gamma_3P1_Hz.*sqrt(inten./(2.*analyVar.satInten_3P1));

%% Loop through each point in scan
% Filenames for scaope data
traceFile = cellfun(@(x) [analyVar.dataDir char(x) fileSuffix],indivDataset{basenameNum}.fileAtom,'UniformOutput',0);
% Initialize loop variables
[inten1, inten2] = deal(zeros(1,indivDataset{basenameNum}.CounterAtom));
parfor k = 1:indivDataset{basenameNum}.CounterAtom;
    % Retrieve data from scope
    scopeData = read_TDS2024(traceFile{k});
    % Get on and off logical indexes
    [onInd] = scopeTimeInd(expTime,scopeData.time.axis,scopeData.time.inc);
    % Calculate measured intensities
    inten1(k) = PD_1_calib*mean(scopeData.ch1.value(onInd));
    inten2(k) = PD_2_calib*mean(scopeData.ch2.value(onInd));
end

%% Calculate measured Rabi frequency
omgBragg1 = omega(inten1); omgBragg2 = omega(inten2);
omgEff    = omgBragg1.*omgBragg2./(2*analyVar.braggDetune);
%% Calculate scale factor
scale = sin(mean(omgEff)*expTime/2)^2./(sin(omgEff.*expTime/2).^2);

%% Scale upper state population
if iscell(numBragg)
    numBragg = cellfun(@(x) x.*scale,numBragg,'UniformOutput',0);
else
    numBragg = numBragg.*scale;
end

%% Plot intensities 
if analyVar.plotBraggInt
    figure(analyVar.figNum.braggInt)
    intHan = plot([inten1; inten2]',analyVar.MARKERS{basenameNum});
    hold on
    grid on
    set(intHan,...
        'Color'     ,   analyVar.COLORS(basenameNum,:)   ,...
        'MarkerSize',   analyVar.markerSize              );
    xlabel('Number of Observations',...
        'FontSize'  ,   analyVar.labelFontSize  ,...
        'FontWeight',   'bold'                  );
    ylabel('Photodiode Voltage [V]',...
        'FontSize'  ,   analyVar.labelFontSize  ,...
        'FontWeight',   'bold'                  );
    set(gca,...
        'LineWidth' ,   2                       ,...
        'FontSize'  ,   analyVar.axisFontSize   ,...
        'FontWeight',   'bold'                  );
    set(gcf,...
        'Name'  ,   'Bragg Photodiode Voltages' )
    set(get(get(intHan(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','Off')
    if basenameNum == analyVar.numBasenamesAtom
        legend(num2str(analyVar.timevectorAtom),...
            'Location'  ,   'Best'  );
    end
end
end