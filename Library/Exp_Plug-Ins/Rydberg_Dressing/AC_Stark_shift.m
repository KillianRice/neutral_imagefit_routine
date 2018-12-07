function UV_detuning_plot_6(data)
%% Initialize Variables
a            = length(data.avgDataset);
[~,orderVec] = sort(data.analyVar.uniqScanList);
nlmC         = cell(1,a);
offset       = -0.3;
[dataHanNum,dataHanFit] = deal(cell(1,a));
% legStr     = {
%     '-5 MHz'
%     '-4 MHz'
%     '-3 MHz'
%     '-2 MHz'
%     '-1 MHz'
%     '0 MHz'
%     '1 MHz'
%     '2 MHz'
%     '3 MHz'
%     '4 MHz'
%     '5 MHz'
%     '6 MHz'
%     '7 MHz'
%     '8 MHz'};

% Fitted data from spectrum fit
detune = data.analyVar.uniqScanList;
AoC    = data.amplitude(:,1);
cntr   = data.lineCenter(:,1);
width  = data.fullWidth(:,1);

% Plotting variables
COLORS     = {'k' 'r' 'b' 'm' 'g' [0 0.75 0.75] [0.6 0.2 1] [1 0.588 0] [0.5 0.5 0.5]};
MARKER     = {'o' 's' 'p' '>' 'd' '*' '<' 'h' '^'};
numFig     = figure;
fitFig     = figure;
AoCFig     = figure;
widthFig   = figure;
cntrFig    = figure(104);

% Fitting function
gauss      = @(coeffs,x) coeffs(4)*exp(-coeffs(1)/coeffs(3)*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2)));
dGauss     = @(coeffs,x) coeffs(4) - coeffs(1)*(x-coeffs(3))/(coeffs(2))^2.*exp(-(x - coeffs(3)).^2./(2*coeffs(2).^2));
dLorentz   = @(coeffs,x) coeffs(4) - 2*coeffs(1)*(x-coeffs(3))*(coeffs(2))^2./((x - coeffs(3)).^2 + coeffs(2).^2).^2;

%% Guessing and Fitting
for i = 1:a
    ind = data.avgDataset{i}.simScanIndVar*1e6;
    dep = data.avgDataset{i}.avgTotNum;
    
    % Initial Guesses
    smoothNum        = smooth(dep,'sgolay',1); % smooth for helping make guesses
    [minVal, minLoc] = min(smoothNum);         % find location of minimum
    
    initOffset = mean(dep([1:5 end-5 end])); % Mean of first 5 and last 5 points
    initCenter = ind(minLoc);                % Location of minimum in smoothed data
    initWidth  = sqrt(abs(sum((ind-initCenter).^2.*dep)./sum(dep)))/10;
    initAmp    = (initOffset - minVal);      % Value below guessed offset
    beta0      = [initAmp initCenter initWidth initOffset];
    
    % Fitting
    nlmC{i} = NonLinearModel.fit(ind,dep,gauss,beta0,...
        'CoefficientNames',{'Amplitude','Line Center','Sigma','Number Offset'});
end

%% Plotting
for j = 1:a
    i = orderVec(j);
    % Rename data into more concise variables
    ind      = data.avgDataset{i}.simScanIndVar;
    dep      = data.avgDataset{i}.avgTotNum;
    setCntr  = 90.2381;
    N0       = nlmC{i}.Coefficients.Estimate(4);
    rangeInd = minmax(cell2mat(cellfun(@(x) x.simScanIndVar,data.avgDataset,'UniformOutput',0)')');
    fitInd   = linspace(rangeInd(1),rangeInd(2),1e3);
    fitDep   = nlmC{i}.predict(fitInd'*1e6)./N0;
    depNorm  = (dep./N0);
    
    legStr{j} = [num2str(data.analyVar.uniqScanList(i)) ' MHz'];
    
    % Pick value of color and marker to use
    indIter = mod(j,length(COLORS)); if indIter == 0; indIter = length(COLORS); end
    
    % Data Curves with fit
    figure(numFig);
    dataHanNum{j} = plot(ind-setCntr,depNorm + (j-1)*offset,fitInd-setCntr,fitDep + (j-1)*offset); hold on;
    set(dataHanNum{j}(1),...
        'Color'             ,   COLORS{indIter} ,...
        'LineStyle'         ,   'none'          ,...
        'Marker'            ,   MARKER{indIter} ,...
        'MarkerSize'        ,   13              ,...
        'MarkerFaceColor'   ,   COLORS{indIter} )
    set(dataHanNum{j}(2),...
        'Color'     ,   COLORS{indIter} ,...
        'LineStyle' ,   ':'             ,...
        'LineWidth' ,   2.5             )
    set(get(get(dataHanNum{j}(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','Off')
    
    % Only show Fit curves
    figure(fitFig)
    dataHanFit{j} = plot(fitInd-setCntr,fitDep); hold on
    set(dataHanFit{j},'Color',COLORS{indIter},'LineWidth',2.5)
end

figure(numFig)
axis tight; 
grid on
set(gca,...
    'FontSize'      ,   35              ,...
    'FontWeight'    ,   'Bold'          ,...
    'YTickLabel'    ,   []              ,...
    'LineWidth'     ,   2               )
xlabel('689 Detuning [MHz]')
ylabel('Relative Atom Loss')
legHan = legend(cellfun(@(x) x(1),dataHanNum),legStr,...
    'location'  ,   'NorthWest'  );
set(legHan,...
    'FontSize'  ,   16  )

figure(fitFig)
axis tight; grid on
curYLim = get(gca,'YLim');
set(gca,...
    'YLim'          ,   [0 curYLim(2)]  ,...
    'FontSize'      ,   35              ,...
    'FontWeight'    ,   'Bold'          ,...
    'LineWidth'     ,   2               )
xlabel('689 Detuning [MHz]')
ylabel({'Normalized Loss'})
% legHan = legend(cellfun(@(x) x(1),dataHanFit),legStr,...
%     'location'  ,   'Best'  );
% set(legHan,...
%     'FontSize'  ,   16  )


figure(widthFig)
widthHan = plot(detune,width*1e-3);
hold on
set(widthHan(1),...
    'Marker'            ,   'o'     ,...
    'MarkerSize'        ,   15      ,...
    'MarkerFaceColor'   ,   [1 0.588 0]     ,...
    'LineStyle'         ,   'none'  ,...
    'Color'             ,   'k'     )
set(gca,...
    'FontSize'      ,   35              ,...
    'FontWeight'    ,   'Bold'          ,...
    'LineWidth'     ,   2               )
xlabel('UV Detuning [MHz]')
ylabel('FWHM [kHz]')
legend(widthHan,{'11mW UV'},...
    'location'  ,   'NorthEast'  );
grid on

figure(AoCFig)
AoCHan = plot(detune,AoC);
set(AoCHan(1),...
    'Marker'            ,   'o'     ,...
    'MarkerSize'        ,   15      ,...
    'MarkerFaceColor'   ,   'b'     ,...
    'LineStyle'         ,   'none'  ,...
    'Color'             ,   'k'     )
set(gca,...
    'FontSize'      ,   35              ,...
    'FontWeight'    ,   'Bold'          ,...
    'LineWidth'     ,   2               )
xlabel('UV Detuning [MHz]')
ylabel('Area Under the Curve')
legend(AoCHan,{'11mW UV'},...
    'location'  ,   'SouthWest'  );
grid on


figure(cntrFig)
cntrHan  = plot(detune,cntr*1e-6 - setCntr);
set(cntrHan(1),...
    'Marker'            ,   'o'     ,...
    'MarkerSize'        ,   15      ,...
    'MarkerFaceColor'   ,   [1 0.588 0]     ,...
    'LineStyle'         ,   'none'  ,...
    'Color'             ,   'k'     )
set(gca,...
    'FontSize'      ,   35              ,...
    'FontWeight'    ,   'Bold'          ,...
    'LineWidth'     ,   2               )
xlabel('UV Detuning [MHz]')
ylabel('Line Center [MHz]')
legend(cntrHan,{'11mW UV'},...
    'location'  ,   'NorthWest'  );
grid on

% Fit to derivative of lorentzian
if 0
    nlmL = NonLinearModel.fit(detune,cntr*1e-6 - setCntr,dLorentz,[1 4 2.5 0.05],...
        'CoefficientNames',{'Amplitude','HWHM','Line Center','Offset'})
    %nlmG = NonLinearModel.fit(detune,cntr*1e-6 - setCntr,dGauss,[1 4 0 -0.1]);
    hold on; axis tight
    curXLim = get(gca,'XLim');
    fitInd  = linspace(curXLim(1)-0.25,curXLim(2)+0.25,1e3);
    derivHan = plot(fitInd,nlmL.predict(fitInd'));%,fitInd,nlmG.predict(fitInd'),'r');
    set(derivHan(:),...
        'Color'             ,   [1 0.588 0]   ,...
        'LineStyle'         ,   '--'    ,...
        'LineWidth'         ,   3       )
    set(get(get(derivHan,'Annotation'),'LegendInformation'),'IconDisplayStyle','Off')
    axis tight; set(gca,'XLim',[curXLim(1)-0.25 curXLim(2)+0.25])
end
end