function [plotHan axHan] = default_plot(analyVar,iter,figNum,label,titleCell,legData,xData,yData,varargin)
% Function to create default plots of instantaneous variables. Mostly used
% to cleanup code and eliminate redundancy of code.
%
% INPUTS:
%   analyVar - structure of all pertinent variables for the imagefit
%              routines
%   iter     - 2 element vector [i j] where 
%                   i - value of current iteration (used for finding data in vectors)
%                   j - limit of iterations (used for plotting legend and title)
%   figNum   - figure where the data will be plotted
%   xData    - independent variable (along x)   
%   yData    - dependent variable (along y)
%   label    - label along y-axis
%   legData  - Vector to extract legend information from
%   title    - Cell containing 
%
% OUTPUTS:
%   plotHan  - handle to the axes plotted in case modification is required
%   Will create plot on figure figNum

%% Set booleans on Properties called
flagHideLeg = 0; plotErr = 0;
for propIter = 1:length(varargin)
    if ischar(varargin{propIter})
        switch varargin{propIter}
            case 'Hide Legend'
                flagHideLeg = 1;
                dataHideLeg = propIter + 1;
            case 'Errorbar'
                plotErr = 1;
                errData = varargin{propIter + 1};
        end
    end
end

%% Call figure
figure(figNum);
% Make background white
set(gcf,'color','white');

% Set the update function so the data cursor has more digits
dcm_obj = datacursormode(figNum);
set(dcm_obj,'UpdateFcn',@default_cursorUpdateFcn);
% Unfortunately, I am not sure how to make the datacursor still show the
% values for error bars. Perhaps modifying the updatefcn is how but I am
% not sure :(

%% Check size of title
% If title is 1 element make this the main title
% if title is equal to the number of data vectors then use all
% else ignore title
if length(titleCell) == 1 
    mainTitle = titleCell{1};
    titleCell = cell(1,size(yData,1));
elseif length(titleCell) == size(yData,1)+1
    mainTitle = titleCell{1};
    titleCell = titleCell(2:end);
else
    mainTitle = '';
    titleCell = cell(1,size(yData,1));
end

%% Initialize variables
[plotHan axHan] = deal(zeros(1,size(yData,1)));

%% Set subplot based on number of rows of yData
for i = 1:size(yData,1)
    % Developed to compare two parameters in single plot (i.e. 2hk lattice number)
    axHan(i) = subplot(size(yData,1),1,i); hold on; grid on;

    %% Plot data
    indMarker  = mod(iter(1) - 1,length(analyVar.MARKERS)) + 1;
    indColor   = mod(iter(1) - 1,length(analyVar.COLORS)) + 1;
    plotHan(i) = plot(xData(i,:),yData(i,:),analyVar.MARKERS{indMarker},...
        'Color'     ,   analyVar.COLORS(indColor,:)  ,...
        'MarkerSize',   analyVar.markerSize         );
    if plotErr == 1
        try
        hErr = errorbar(xData(i,:),yData(i,:),errData(i,:)  ,...
            'Color'     ,   analyVar.COLORS(indColor,:)     ,...
            'MarkerSize',   analyVar.markerSize             );
        catch
            i;
        end
        set(get(get(hErr,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end

    %% Axes labels and Properties
    xlabel(analyVar.xDataLabel,...
        'FontSize'  ,   analyVar.labelFontSize  ,...
        'FontWeight',   'bold'                  );
    ylabel(label{i},...
        'FontSize'  ,   analyVar.labelFontSize  ,...
        'FontWeight',   'bold'                  );
    set(gca,...
        'LineWidth' ,   1                       ,...
        'FontSize'  ,   analyVar.axisFontSize   ,...
        'FontWeight',   'bold'                  ,...
        'Box'       ,   'on'                    );

    %% Axes limits
    % Trim mean trims away the TrimPercent and returns the average of the
    % remaining data. This helps to focus on the relevant data instead of
    % throwing off the window due to a poor fit.
    try
        ylim(trimmean(yData(i,:),analyVar.ylimMeanTrimPercent*1e2)...
            .*[(1 - analyVar.yPlotLimBounds),(1 + analyVar.yPlotLimBounds)]);
    catch
        fprintf('There was a problem changing the y-limits for %s. \t',label{i})
        disp('Skipping changing yLim\n')
    end
    
    %% Subtitle
    title(titleCell{i},...
        'FontSize'  ,   analyVar.titleFontSize   ,...
        'FontWeight',   'bold'                  );

    %% On last scan file insert legend enumerating all scans shown
    if iter(1) == iter(2) && i == size(yData,1)
        % Create legend
        if flagHideLeg
            axChild = sort(get(gca,'Children'));
            set(get(get(axChild(dataHideLeg),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        end
        legend(num2str(legData(1:iter(1))),'Location','Best');
        
        % Set Figure window title
        if iscell(label{i})
            name = label{i}{1};
        else
            name = label{i};
        end
        set(gcf,'Name',[name ': Time = ' num2str(analyVar.timevectorAtom(iter(1)))])
        
        % Set main title
        if ~isempty(mainTitle)
            mtit([mainTitle ' using ' strrep(analyVar.fitModel,analyVar.InitCase,[analyVar.InitCase ' '])],...
                'FontSize',analyVar.titleFontSize,'zoff',0,'xoff',-.01)
        end
    end
end
end