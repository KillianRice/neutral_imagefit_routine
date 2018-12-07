function StructOut = Wavemeter(lcl_analyVar,lcl_indivDataset,lcl_avgDataset)
% Performs analysis on wavemeter readings taken from each individual scan.
% Works in conjunction with param_extract_ind_var_to_wavemeter
%
%STATUS: 2018.08.22: Does not work with averaging or normalization of atom number, so turn off the
%"normMeanNum" variable in AnalysisVariables when using.
% INPUTS:
%
% OUTPUTS:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variale extraction from batch files
%%%%%%%%%%%%%%

% Decide to use averaged data or individual scans
avgAutoFlag = length(cell2mat(lcl_analyVar.posOccurUniqVar)) > length(lcl_analyVar.posOccurUniqVar) & lcl_analyVar.uniqScanList ~= 0;

% Extract independent variables and number from correct structure
indField   = {'imagevcoAtom' 'simScanIndVar'};  % fields of independent variables in indivDataset and avgDataset
numField   = {'winTotNum' 'avgTotNum'};         % number data in indivDataset and avgDataset
labelField = {'timevectorAtom' 'uniqScanList'}; % label used for each plot


curDataset = lcl_indivDataset; fieldIndx = 1;

% Save variables from cell of structures, also apply any scaling functions.
% Input x, applies the function x.field to every cell in the array curDataseet.
%I.e. it constructs curDataset(x).indField{fieldIndx} for all x elements.
indVar  = cellfun(@(x) x.(indField{fieldIndx}),curDataset,'UniformOutput',0);
partNum = cellfun(@(x) x.(numField{fieldIndx})',curDataset,'UniformOutput',0);
label   = lcl_analyVar.(labelField{fieldIndx});

for i=1:length(indVar)
    indVar{i} = indVar{i}-21698.46;
end

%% Lorentzian fitting and packaging of results
%%%%%%%%%%%%%%
fitLorz = 0; % Flag for choosing to try and fit the Lorentzians to data, or not.
if fitLorz
    numLorz = 4;         % Number of expected loss peaks to be fit
    for i = 1:length(indVar)
        [indivCoeffMtx,indivFitStruct,indivWaveNum] = multiLorzFit(indVar{i},partNum{i},numLorz);
        coeffStruct(i)        = {indivCoeffMtx};
        fitStruct(i)          = {indivFitStruct};
        indivWaveNumStruct(i) = {indivWaveNum};
    end
    
    % Fit Coefficient Analysis
    coeffAvgMtx = zeros(numLorz,5);
    coeffSTDMtx = zeros(numLorz,5);
    scanVars    = zeros(1,5);
    for k = 1:numLorz               % Loop over number of Lorentzians
        for j = 1:5                 % Loop over number of fit coefficients
            for i = 1:length(indVar)% Loop over scans
                scanVars(i)  = coeffStruct{i}(k,j);
            end
            coeffAvgMtx(k,j) = mean(scanVars);
            coeffSTDMtx(k,j) = std(scanVars);
        end
    end
end

%% PLOTTING and Data Visualization
%%%%%%%%%%%%%%
figure % Plot each scan separately, each in a seperate subplot.
for i = 1:length(indVar)
    subplot(ceil(length(indVar)/2),2,i);
    plot(indVar{i},partNum{i},'-','LineWidth',2);
    %plot(smooth(indVar{i}),smooth(partNum{i}),'-','LineWidth',2);
    hold on
    plot(indVar{i},partNum{i},'.','LineWidth',2);
    ylabel('Atom Number','FontSize',13,'FontWeight','bold');
    xlabel('Wavemeter Wavelength (cm^-1)','FontSize',10,'FontWeight','bold');
    grid on;  grid minor;
    set(gca,'FontSize',18);
    axis tight
    %     coeffText = sprintf('Fit Coeff \n%5.4e\n%5.4e\n%5.4e',p(1),p(2),p(3));
    %     fig=subplot(ceil(numScans/3),3,scanNum);
    %     text(0.5,0.2,coeffText      ,...
    %         'Units'            , 'normalized'   ,...
    %         'Parent'           , fig            ,...
    %         'LineStyle'        , '-'            ,...
    %         'BackgroundColor'  , [1,1,0.84]     ,...
    %         'EdgeColor'        , [0,0,0]        ,...
    %         'FontSize'         , 10             );
    str   = sprintf('Each Scan Separately, Scan # %d',i);
    title(str);
    if fitLorz  % If fitting Lorentzians, plot them on the data
        for k=1:numLorz
            hold on
            plot(indivWaveNumStruct{i}{k},fitStruct{i}{k},'LineWidth',4);
        end
    end
end

figure % Plot all scans on same figure (but individually) for comparison.
for i = 1:length(indVar)
    plot(indVar{i},partNum{i},'-','LineWidth',2);
    %plot(smooth(indVar{i}),smooth(partNum{i}),'-','LineWidth',2);
    hold on
end
ylabel('Atom Number','FontSize',23,'FontWeight','bold');
xlabel('Wavemeter Wavelength (cm^-1)','FontSize',20,'FontWeight','bold');
grid on;  grid minor;
set(gca,'FontSize',20);
axis tight
str   = sprintf('All Individual Scans, Distinct Colors');
title(str);

% figure % Concatenate, sort, and plot all data on top of each other.
% bigindVar = []; bigpartNum = [];
% for i = 1:length(indVar)
%     bigindVar  = vertcat(bigindVar,indVar{i});
%     bigpartNum = vertcat(bigpartNum,partNum{i});
%     %bigpartNum = vertcat(bigpartNum,partNum{i}-(mean(partNum{i})));
%     % bigindVar = smooth(vertcat(bigindVar,indVar{i}),'sgolay');
%     % bigpartNum = smooth(vertcat(bigpartNum,partNum{i}-(mean(partNum{i}))),'sgolay');
%     % bigpartNum = smooth(vertcat(bigpartNum,partNum{i}),'sgolay');
% end
% M = sortrows([bigindVar bigpartNum]);
% plot(M(:,1),M(:,2),'b-','LineWidth',2);
% hold on
% plot(M(:,1),M(:,2),'.','LineWidth',2);
% ylabel('Atom Number','FontSize',18,'FontWeight','bold');
% xlabel('Wavemeter Wavelength (cm^-1)','FontSize',15,'FontWeight','bold');
% grid on;  grid minor;
% set(gca,'FontSize',20);
% str   = sprintf('Data From All Scans Aggregated');
% title(str);


figure % Concatenate all scans into large vectors {X,Y}.
bigIndVar = []; bigPartNum = [];
for i = 1:length(indVar)
    indVar{i}  = indVar{i}+21698.46; % Convert from detuning to bare wavenumbers
    bigIndVar  = vertcat(bigIndVar,indVar{i});
    bigPartNum = vertcat(bigPartNum,partNum{i});
end
% Bin and average the combined data 
binSize    = 30; % In MHz
nBins      = ceil((bigIndVar(end)-bigIndVar(1))/(binSize*(0.001/29.979))); % Number of bins such that they are binSize wide.
[N,edges,bins] = histcounts(bigIndVar,nBins); % The binning
avgBinX = zeros(nBins,1); % Preallocate average vectors 
avgBinY = zeros(nBins,1);
stdBinY = zeros(nBins,1);
for i=1:nBins
    avgBinX(i) = mean(bigIndVar(bins==i));
    avgBinY(i) = mean(bigPartNum(bins==i));
    stdBinY(i) = std(bigPartNum(bins==i));
end
avgBinX(isnan(avgBinX)) = 0; % If no points in that bin, set the value to zero.
avgBinY(isnan(avgBinY)) = 0;
stdBinY(isnan(avgBinY)) = 0;
hold on
errorbar(avgBinX,avgBinY,stdBinY,'b*-','LineWidth',2);
%plot(avgBinX,avgBinY,'b*-','LineWidth',2);
ylabel('Binned Atom Number','FontSize',18,'FontWeight','bold');
xlabel('Binned Wavemeter Wavelength (cm^-1)','FontSize',15,'FontWeight','bold');
grid on;  grid minor;
set(gca,'FontSize',20);
str   = sprintf('Binned and Averaged Total Dataset');
title(str)

% Add vertical lines at positions of 84Sr resonances 
line1=21698.014; 
line2=21698.049; 
line3=21698.078; 
yL = get(gca,'YLim');
line([line1 line1],yL,'Color','r'); 
line([line2 line2],yL,'Color','r'); 
line([line3 line3],yL,'Color','r'); 


% Plot the wavemeter data and interpolations
figure
for basenameNum = 1:lcl_analyVar.numBasenamesAtom
    ind_var = lcl_indivDataset{basenameNum}.rawimagevcoAtom;
    wavemeter = lcl_indivDataset{basenameNum}.imagevcoAtom;
    p = lcl_indivDataset{basenameNum}.indCalib;
   
    subplot(ceil(lcl_analyVar.numBasenamesAtom/3),3,basenameNum);
    hold on
    plot(ind_var,wavemeter,'k*','LineWidth',5)
    plot(ind_var,polyval(p,ind_var),'r-','LineWidth',3)
    title(['Wavemeter Interpolation: Scan #',num2str(basenameNum)]);
    ylabel('Wavemeter Wavelength','FontSize',13,'FontWeight','bold');
    xlabel('Laser Control Voltage (V)','FontSize',10,'FontWeight','bold');
    grid on;  grid minor;
    set(gca,'FontSize',18);
    coeffText = sprintf('Interpolation Coefficients\n%5.4e\n%4.3e\n%5.4e',p(1),p(2),p(3));
    fig=subplot(ceil(length(indVar)/3),3,basenameNum);
    text(0.5,0.2,coeffText      ,...
        'Units'            , 'normalized'   ,...
        'Parent'           , fig            ,...
        'LineStyle'        , '-'            ,...
        'BackgroundColor'  , [1,1,0.84]     ,...
        'EdgeColor'        , [0,0,0]        ,...
        'FontSize'         , 10             );
end


%% Pack workspace into a structure for output
%%%%%%%%%%%%%%
% If you don't want a variable output, prefix it with lcl_
StructOut = who();
StructOut = v2struct(cat(1,'fieldNames',StructOut(cellfun('isempty',regexp(StructOut,'\<lcl_')))));
end


