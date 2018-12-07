function [braggPop,iPkSub,jPkSub] = numIntBraggPop(analyVar,indivDataset,basenameNum,numPeak)
% Function to load and subtract the central peak fit from OD
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit routines
%   indivDataset - Cell of structures containing all scan/batch specific data
%
% OUTPUTS:
%   braggPop - Vector of the numerical integrated population in the Bragg peak

%% Initialize variables
[braggPop, iPkSub, jPkSub] = deal(cellfun(@(x) nan(1,indivDataset{basenameNum}.CounterAtom), cell(length(numPeak),1),'UniformOutput',0));

%% Loop through each batch file listed in basenamevectorAtom
% Processes all the image files in the current batch
for k = 1:indivDataset{basenameNum}.CounterAtom;
    
    %% Initialize roi image to size
    %%%%%%%-----------------------------------%%%%%%%%%%
    roiDistImage = zeros(analyVar.funcPlotWin(basenameNum));
    
    %% Retrieve windowed image cell
    %%%%%%%-----------------------------------%%%%%%%%%%
    OD_Fit_ImageCell = cellfun(@(x) reshape(indivDataset{basenameNum}.All_OD_Image{k}(x),[1 1]*analyVar.funcFitWin(basenameNum)),...
        cellfun(analyVar.fitWinLogicInd,indivDataset{basenameNum}.roiWin_Index,'UniformOutput',0),'UniformOutput',0);
    
    %% Retrieve fit coefficients
    %%%%%%%-----------------------------------%%%%%%%%%%
    PCell = indivDataset{basenameNum}.All_PCell{k};
    
    %% Get fit center in full image and find one sigma bounds
    %%%%%%%-----------------------------------%%%%%%%%%%
    CloudCntr  = round((analyVar.roiWinRadAtom(basenameNum) - (analyVar.funcFitWin(basenameNum) - 1)/2)...
            + [indivDataset{basenameNum}.All_fitParams{k}{1}.xCntr indivDataset{basenameNum}.All_fitParams{k}{1}.yCntr]);
    oneSigBnds = [round(CloudCntr - cellfun(@(x) indivDataset{basenameNum}.All_fitParams{k}{1}.(x),{'sigX' 'sigY'})); ...
                  round(CloudCntr + cellfun(@(x) indivDataset{basenameNum}.All_fitParams{k}{1}.(x),{'sigX' 'sigY'}))];
    
    %% Weighting (to use fitModel function correctly)
    %%%%%%%-----------------------------------%%%%%%%%%%
    % Return a weighting cell if it was used in fitting
    errCell = cellfun(@(x,y) get_OD_weight(x(end),y),PCell,OD_Fit_ImageCell,'UniformOutput',0);
    
    %% Calculate the estimated distribution from the specified model
    %%%%%%%-----------------------------------%%%%%%%%%%
    % Generate spatial coordinates (independent variables)
    [Xgrid,Ygrid] = meshgrid(1:analyVar.funcFitWin(basenameNum));
    
    % Call fit model and return distribution
    fitDistCell = cellfun(@(x,y)(reshape(feval(str2func(analyVar.fitModel),x(1:length(analyVar.InitCondList)),...
        [Xgrid(:) Ygrid(:) y(:)]), [1 1].*analyVar.funcFitWin(basenameNum))),PCell,errCell,'UniformOutput',0);
    
    %% Collect individual fit windows into single image
    %%%%%%%-----------------------------------%%%%%%%%%%
    % Temporary variable assignment
    roiDistTmp = indivDataset{basenameNum}.roiWin_Index{1};
    
    % Index the pnts in the ROI window that define the fit window and assign the fit values
    roiDistTmp(analyVar.fitWinLogicInd(roiDistTmp)) = fitDistCell{1};
    
    % Sum all cloud images (fit windows) together
    roiDistImage = roiDistImage + roiDistTmp;
    
    %% Subtract the estimated distribution from the experimental data
    %%%%%%%-----------------------------------%%%%%%%%%%
    remCloudOD = indivDataset{basenameNum}.All_OD_Image{k} - roiDistImage;
    
    %% Ignore condensate region and smooth data for fitting
    %%%%%%%-----------------------------------%%%%%%%%%%
    % Get axes of the condensate 
    a = indivDataset{basenameNum}.All_fitParams{k}{1}.sigX;
    b = indivDataset{basenameNum}.All_fitParams{k}{1}.sigY;
    [X,Y] = meshgrid(1:size(remCloudOD,2),1:size(remCloudOD,1));
    % Create logical mask of area around condensate
    maskCloud = ((X-CloudCntr(1))/a).^2 + ((Y-CloudCntr(2))/b).^2 > 1;
    % Smooth the area around condensate to identify the bragg peak position
    smoothRemCloudOD = conv2(remCloudOD.*maskCloud,fspecial('gaussian',analyVar.gaussFiltAmp,analyVar.gaussFiltSig),'same');
    
    % Loop through peaks
    for braggIter = 1:length(numPeak)
        %% Restrict data to bragg sign specified
        if analyVar.winID(basenameNum) == 1 || (analyVar.winID(basenameNum) == 2 && numPeak(braggIter) == 1)
            regFiltOD    = remCloudOD(1:CloudCntr(2),oneSigBnds(1,1):oneSigBnds(2,1));
            smoothFiltOD = smoothRemCloudOD(1:CloudCntr(2),oneSigBnds(1,1):oneSigBnds(2,1));
        end
        
        if analyVar.winID(basenameNum) == -1 || (analyVar.winID(basenameNum) == 2 && numPeak(braggIter) == -1)
            regFiltOD    = remCloudOD(CloudCntr(2):end,oneSigBnds(1,1):oneSigBnds(2,1));
            smoothFiltOD = smoothRemCloudOD(CloudCntr(2):end,oneSigBnds(1,1):oneSigBnds(2,1));
        end

        %% Estimate position of Bragg peak
        %%%%%%%-----------------------------------%%%%%%%%%%
        sizeFilt = size(smoothFiltOD);
        % Guess peak position
        [~,peakInd] = max(smoothFiltOD(:));
        [iPkSub{braggIter}(k),jPkSub{braggIter}(k)] = ind2sub(sizeFilt,peakInd);
        % Check where max was found as it should be approximately fixed away from the central peak
        if analyVar.winID(basenameNum) == 1 || (analyVar.winID(basenameNum) == 2 && numPeak(braggIter) == 1)
            if abs(analyVar.braggCalib*analyVar.droptimeAtom(basenameNum) - iPkSub{braggIter}(k)) > 10
                iPkSub{braggIter}(k) = CloudCntr(2) - analyVar.braggCalib*analyVar.droptimeAtom(basenameNum);
                jPkSub{braggIter}(k) = round(sizeFilt(2))/2; % Force the x position to be near where the CloudCntr was if no good max
            end
        end
        
        if analyVar.winID(basenameNum) == -1 || (analyVar.winID(basenameNum) == 2 && numPeak(braggIter) == -1)
            if abs(CloudCntr(2) - analyVar.braggCalib*analyVar.droptimeAtom(basenameNum) - iPkSub{braggIter}(k)) > 10
                iPkSub{braggIter}(k) = analyVar.braggCalib*analyVar.droptimeAtom(basenameNum);
                jPkSub{braggIter}(k) = round(sizeFilt(2))/2; % Force the x position to be near where the CloudCntr was if no good max
            end
        end

        %% Calculate Bragg mask and numerically integrate
        %%%%%%%-----------------------------------%%%%%%%%%%
        % Compute indices within a radius of the center positions
        [X,Y] = meshgrid(1:sizeFilt(2),1:sizeFilt(1));
        r     = analyVar.braggWinRad;
        mask  = (X-jPkSub{braggIter}(k)).^2 + (Y-iPkSub{braggIter}(k)).^2 < r^2;
        % Computes number based of summed OD
        braggPop{braggIter}(k) = analyVar.sizefactor^2*sum(regFiltOD(mask))/analyVar.AbsCross;
    end
end