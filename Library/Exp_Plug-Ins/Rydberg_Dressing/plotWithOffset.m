function plotWithOffset(analyVar,avgDataset)
% Plot data with offset for clarity

%% Initialize Variables
a            = length(avgDataset);
[~,orderVec] = sort(analyVar.uniqScanList); orderVec = flipud(orderVec);
numFig       = figure;
offset       = 1;
[dataHanNum,legStr] = deal(cell(1,a));
xRange       = [-2 2];

% Plotting variables
COLORS     = {'k' 'r' 'b' 'm' 'g' [0 0.75 0.75] [0.6 0.2 1] [1 0.588 0] [0.5 0.5 0.5]};
MARKER     = {'o' 's' 'p' '>' 'd' '*' '<' 'h' '^'};

%% Plotting
for j = 1:a
    i = orderVec(j);
    % Rename data into more concise variables
    ind      = analyVar.funcDataScale(avgDataset{i}.simScanIndVar);
    dep      = avgDataset{i}.avgTotNum;
    N0       = max(dep);
    depNorm  = (dep./N0);
    
    legStr{j} = [num2str(analyVar.uniqScanList(i)) ' \mus'];
    
    % Pick value of color and marker to use
    indIter = mod(j,length(COLORS)); if indIter == 0; indIter = length(COLORS); end
    
    % Data Curves with fit
    figure(numFig); subplot(1,2,1);
    vertPos = (a - j - 1)*offset;
    
    dataHanNum{j} = plot(ind,depNorm + vertPos); hold on;
    
    vertLine = depNorm + vertPos; vertLine = mean(vertLine(1:5));
    line(xRange,[1 1].*vertLine,...
        'Color'     ,   COLORS{indIter} ,...
        'LineStyle' ,   '--'            ,...
        'LineWidth',    2)
    set(dataHanNum{j}(1),...
        'Color'             ,   COLORS{indIter} ,...
        'LineStyle'         ,   'none'          ,...
        'Marker'            ,   MARKER{indIter} ,...
        'MarkerSize'        ,   10              ,...
        'MarkerFaceColor'   ,   COLORS{indIter} )
end

figure(numFig)
grid on
set(gca,...
    'FontSize'      ,   35              ,...
    'FontWeight'    ,   'Bold'          ,...
    'YTickLabel'    ,   []              ,...
    'LineWidth'     ,   2               )
xlabel('^3P_1 Detuning [MHz]')
ylabel('Relative Ground State Population')
legHan = legend(cellfun(@(x) x(1),dataHanNum),legStr,...
    'location'  ,   'NorthWest'  );
set(legHan,...
    'FontSize'  ,   16  )
xlim(xRange)
ylim([0 a+0.15]-1)
end