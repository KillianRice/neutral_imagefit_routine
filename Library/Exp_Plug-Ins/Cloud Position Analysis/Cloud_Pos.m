function posOut = Cloud_Pos(analyVar,indivDataset,avgDataset)
% Allow plotting of the peak position found through fitting the 2D distribution

%% Decide to average data or use individual scans
avgAutoFlag = length(cell2mat(analyVar.posOccurUniqVar)) > length(analyVar.posOccurUniqVar) & ~any(analyVar.uniqScanList == 0);

%% Initialize loop variables
[cntrXCell, cntrYCell] = deal(cell(1,analyVar.numBasenamesAtom));
indivFigNum = figure; if avgAutoFlag; avgFigNum = figure; end
histoFigNum = figure;

%% Loop through each batch file listed in basenamevectorAtom
for basenameNum = 1:analyVar.numBasenamesAtom
    % Preallocate nested loop variables
    [cntrX, cntrY]...
        = deal(NaN(1,indivDataset{basenameNum}.CounterAtom));
    
    for k = 1:indivDataset{basenameNum}.CounterAtom;
        %% Extract cloud position
        cntrX(:,k) = (analyVar.roiWinRadAtom(basenameNum) - (analyVar.funcFitWin(basenameNum) - 1)/2)...
            + indivDataset{basenameNum}.All_fitParams{k}{1}.xCntr;
        cntrY(:,k) = (analyVar.roiWinRadAtom(basenameNum) - (analyVar.funcFitWin(basenameNum) - 1)/2)...
            + indivDataset{basenameNum}.All_fitParams{k}{1}.yCntr;
    end
    
    %% Show statistics if enabled
    if analyVar.fitStatFlag
        fprintf('\nStatistics for %g\n####################',analyVar.timevectorAtom(basenameNum))
        fprintf('\nMean in X: %g um',mean(cntrX)*analyVar.pixelsize*1e6)
        fprintf('\nStd. Dev. in X: %g um',std(cntrX)*analyVar.pixelsize*1e6)
        fprintf('\nInitial Velocity in X: %g mm/s \n',std(cntrX)*analyVar.pixelsize/analyVar.droptimeAtom(basenameNum)*1e6)
        
        fprintf('\nMean in Y: %g um',mean(cntrY)*analyVar.pixelsize*1e6)
        fprintf('\nStd. Dev. in Y: %g um',std(cntrY)*analyVar.pixelsize*1e6)
        fprintf('\nInitial Velocity in Y: %g mm/s \n',std(cntrY)*analyVar.pixelsize/analyVar.droptimeAtom(basenameNum)*1e6)
    end
    
    %% Plotting
    % Plot the center Positions
    %%%%%%%-----------------------------------%%%%%%%%%%
    indVar = indivDataset{basenameNum}.imagevcoAtom; % saved for convenience
    cntrLabel = {'Center Position', 'Center Position'};
    cntrTitle = {'', 'X - axis', 'Y - axis'};
    [plotHan axHan] = default_plot(analyVar,[basenameNum analyVar.numBasenamesAtom],...
        indivFigNum,cntrLabel,cntrTitle,analyVar.timevectorAtom,...
        repmat(indVar,1,2)',[cntrX; cntrY]);
    
    
    % Histogram of all the positions
    %%%%%%%-----------------------------------%%%%%%%%%%
    figure(histoFigNum); subplot(2,1,1); hold on;
    hanHistX = histogram((cntrX - mean(cntrX)) * analyVar.pixelsize*1e6          ,...
                  'FaceColor'   ,   analyVar.COLORS(basenameNum,:)               ,...
                  'BinWidth'    ,   analyVar.xbinWidth * analyVar.pixelsize*1e6  );
    xlabel('X axis - Deviation from mean position [\mum]'   ,...
        'FontSize'      ,   analyVar.labelFontSize          ,...
        'FontWeight'    ,   'bold'                          );
    ylabel({'Occurences' [sprintf('(bin size = %g',analyVar.xbinWidth * analyVar.pixelsize*1e6) '\mum)']} ,...
        'FontSize'      ,   analyVar.labelFontSize  ,...
        'FontWeight'    ,   'bold'                  );
    set(gca             ,...
        'LineWidth'     ,   1                       ,...
        'FontSize'      ,   analyVar.axisFontSize   ,...
        'FontWeight'    ,   'bold'                  ,...
        'Box'           ,   'on'                    );
    
    subplot(2,1,2); hold on;
    hanHistY = histogram((cntrY - mean(cntrY)) * analyVar.pixelsize*1e6         ,...
                  'FaceColor'   ,   analyVar.COLORS(basenameNum,:)              ,...
                  'BinWidth'    ,   analyVar.ybinWidth * analyVar.pixelsize*1e6 );
    xlabel('Y axis - Deviation from mean position [\mum]'   ,...
        'FontSize'      ,   analyVar.labelFontSize          ,...
        'FontWeight'    ,   'bold'                          );
    ylabel({'Occurences' [sprintf('(bin size = %g',analyVar.ybinWidth * analyVar.pixelsize*1e6) '\mum)']} ,...
        'FontSize'      ,   analyVar.labelFontSize  ,...
        'FontWeight'    ,   'bold'                  );
    set(gca             ,...
        'LineWidth'     ,   1                       ,...
        'FontSize'      ,   analyVar.axisFontSize   ,...
        'FontWeight'    ,   'bold'                  ,...
        'Box'           ,   'on'                    );
    
    
    %% Oscillation Fitting (if flagged)
    if analyVar.oscFitFlag & ~avgAutoFlag
        figure; 
        if strcmpi('x',analyVar.oscAxis)
            title('X - Axis'); fitCntr = cntrX;
        else
            title('Y - Axis'); fitCntr = cntrY;
        end  
        Spatail_Osc_Fit(indVar*1e-3,fitCntr');
        set(gca,'FontSize',16)
    end
    
    %% Pixel Size Fitting (if flagged)
    if analyVar.pixelSizeFit & ~avgAutoFlag
        % assumes that cloud will fall along the Y-axis
        accelFit = Pixel_Size_Fit(indVar*1e-3,cntrY',analyVar.accelGrav,analyVar.pixelsize);
        theoTime = linspace(min(indVar),max(indVar),1e3);
        figure(indivFigNum); hold on
        set(get(axHan(2),'Children'),...
            'LineStyle' ,   'none'  )
        accelLineHan = plot(theoTime,accelFit.predict(theoTime'*1e-3));
        set(accelLineHan,...
            'LineWidth' ,   2                              ,...
            'Color'     ,   analyVar.COLORS(basenameNum,:) )
    end
    
    %% Save for averaging
    if avgAutoFlag
        % If averaging, then save into cell
        cntrXCell{basenameNum} = cntrX;
        cntrYCell{basenameNum} = cntrY;
    end
end

if avgAutoFlag
    for uniqScanIter = 1:length(analyVar.uniqScanList);
        %% Preallocate nested loop variables
        [matCntrX, matCntrY]...
            = deal(NaN(length(avgDataset{uniqScanIter}.simScanIndVar),...
            length(analyVar.posOccurUniqVar{uniqScanIter})));
        
        %% Loop through all scans that share the current value of the averaging variable
        for simScanIter = 1:length(analyVar.posOccurUniqVar{uniqScanIter})
            % Assign the current scan to be opened
            basenameNum = analyVar.posOccurUniqVar{uniqScanIter}(simScanIter);
            
            cntrX = cntrXCell{basenameNum};
            cntrY = cntrYCell{basenameNum};
            
            % Find the intersection of the scanned variable with the list of all possible values
            % idxShrdInBatch - index of the batch file ind. variables that intersect with
            %                  the set of all ind. variables of similar scans
            % idxShrdWithSim - index of all ind. variables of similar scans that intersect
            %                  with the set of current batch file ind. variables
            % Look at help of intersect if this is unclear
            [~,idxSharedWithSim,idxSharedInBatch] = intersect(avgDataset{uniqScanIter}.simScanIndVar,...
                double(int32(indivDataset{basenameNum}.imagevcoAtom*analyVar.compPrec))*1/analyVar.compPrec);
            
            %% Compute number of atoms
            % Matrix containing the various measured pnts for each scan image
            matCntrX(idxSharedWithSim,simScanIter) = cntrX(:,idxSharedInBatch);
            matCntrY(idxSharedWithSim,simScanIter) = cntrY(:,idxSharedInBatch);
        end
        
        %% Average all like points together
        % Function below takes the mean of the matrix of mixed NaN's and values
        avgCntrX = nanmean(matCntrX,2);
        avgCntrY = nanmean(matCntrY,2);
        
        %% Take the standard deviation of the data
        % Function below takes the standard deviation of the matrix of mixed NaN's and values
        stdCntrX = nanstd(matCntrX,0,2);
        stdCntrY = nanstd(matCntrY,0,2);
        
        %% Plotting
        %%%%%%%-----------------------------------%%%%%%%%%%
        indVar = avgDataset{uniqScanIter}.simScanIndVar;
        
        cntrLabel = {'Avg. Center Position', 'Avg. Center Position'};
        cntrTitle = {'','X-axis','Y-axis'};
        [plotHan axHan] = default_plot(analyVar,[uniqScanIter length(analyVar.uniqScanList)],...
            avgFigNum,cntrLabel,cntrTitle,analyVar.uniqScanList,...
            repmat(indVar,1,2)',[avgCntrX'; avgCntrY'],'Errorbar',[stdCntrX'; stdCntrY']);
        hFig = gcf;
        hFig.CurrentAxes = hFig.Children(2); axis tight
        hFig.CurrentAxes = hFig.Children(3); axis tight
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        
        %% Fitting (if flagged)
        if analyVar.oscFitFlag & avgAutoFlag
            figure;
            if strcmpi('x',analyVar.oscAxis)
                title('X - Axis'); fitCntr = avgCntrX;
            else
                title('Y - Axis'); fitCntr = avgCntrY;
            end
            Spatail_Osc_Fit(indVar*1e-3,fitCntr);
            set(gca,'FontSize',16)
        end
        
        %% Pixel Size Fitting (if flagged)
        if analyVar.pixelSizeFit & avgAutoFlag
            % assumes that cloud will fall along the Y-axis
            accelFit = Pixel_Size_Fit(indVar*1e-3,avgCntrY',analyVar.accelGrav,analyVar.pixelsize);
            theoTime = linspace(min(indVar),max(indVar),1e3);
            figure(avgFigNum); hold on
            set(get(axHan(2),'Children'),...
                'LineStyle' ,   'none'  )
            accelLineHan = plot(theoTime,accelFit.predict(theoTime'*1e-3));
            set(accelLineHan,...
                'LineWidth' ,   2                              ,...
                'Color'     ,   analyVar.COLORS(uniqScanIter,:) )
        end
    end
end
%% Pack workspace into a structure for output
% If you don't want a variable output prefix it with lcl_
posOut = who();
posOut = v2struct(cat(1,'fieldNames',posOut(cellfun('isempty',regexp(posOut,'\<lcl_')))));
end