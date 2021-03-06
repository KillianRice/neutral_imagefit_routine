% Create new figures for plotting in
rawFig    = figure;
pcaFig    = figure;
compFig   = figure;

% Create 2D matricies for image plotting
roiShape = [1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1);
[bordBack, bordBackPCA, bordAtom, bordBackApprox] = deal(zeros(roiShape));

% Logical indexing to specify bord region
bordInd  = indivDataset{1}.roiWin_Index{1} == 0;

% Plot the naive atom image
figure(compFig);

bordAtom(bordInd) = atomNotCloudNaive;
subaxis(1,3,1,...
            'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
            'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
            'Margin'        ,   analyVar.subAxis.margin,...
            'MarginTop'     ,   analyVar.subAxis.marginTop,...
            'Padding'       ,   analyVar.subAxis.padding);
pcolor(bordAtom); shading flat; axis off; title('Original');
cAx = get(gca,'CLim');

bordBackApprox(bordInd) = sum(bsxfun(@times,pcBasis(:,1:varLim),A),2);
subaxis(1,3,2,...
            'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
            'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
            'Margin'        ,   analyVar.subAxis.margin,...
            'MarginTop'     ,   analyVar.subAxis.marginTop,...
            'Padding'       ,   analyVar.subAxis.padding);
pcolor(bordBackApprox); shading flat; axis off; title('Reconstruction');
set(gca,'CLim',cAx)

subaxis(1,3,3,...
            'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
            'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
            'Margin'        ,   analyVar.subAxis.margin,...
            'MarginTop'     ,   analyVar.subAxis.marginTop,...
            'Padding'       ,   analyVar.subAxis.padding);
pcolor(bordAtom - bordBackApprox); shading flat; axis off; colorbar
title(sprintf('Sum Mean Squared Difference: %g',sum((atomNotCloudNaive - sum(bsxfun(@times,pcBasis(:,1:varLim),A),2)).^2)/length(atomNotCloudNaive)));

% Loop through naive background and PCA basis to plot
for k = 1:varLim
    % Get image regions
    bordBack(bordInd)    = backNotCloudNaiveSet(:,k);
    bordBackPCA(bordInd) = pcBasis(:,k);
    
    % Plot backgrounds
    figure(rawFig)
    subaxis(indivDataset{basenameNum}.SubPlotRows,indivDataset{basenameNum}.SubPlotCols,k,...
        'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
        'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
        'Margin'        ,   analyVar.subAxis.margin,...
        'MarginTop'     ,   analyVar.subAxis.marginTop,...
        'Padding'       ,   analyVar.subAxis.padding);
    pcolor(bordBack); shading flat;
    hold on; grid off; axis off
    % set to the same CLim
%     set(gca,'CLim',cAx)
    if k == 1; cAx = get(gca,'CLim');
    else set(gca,'CLim',cAx); end
    
    % Plot PCA basis
    figure(pcaFig)
    subaxis(indivDataset{basenameNum}.SubPlotRows,indivDataset{basenameNum}.SubPlotCols,k,...
        'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
        'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
        'Margin'        ,   analyVar.subAxis.margin,...
        'MarginTop'     ,   analyVar.subAxis.marginTop,...
        'Padding'       ,   analyVar.subAxis.padding);
    pcolor(bordBackPCA); shading flat;
    hold on; grid off; axis off
    % set to the same CLim
%     set(gca,'CLim',cAx)
%     if k == 1; cAx = get(gca,'CLim');
%     else set(gca,'CLim',cAx); end
end

% Plotting the cloud images for diagnosis
        roiShape = [1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1);
        figure
        subaxis(1,2,1,...
            'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
            'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
            'Margin'        ,   analyVar.subAxis.margin,...
            'MarginTop'     ,   analyVar.subAxis.marginTop,...
            'Padding'       ,   analyVar.subAxis.padding);
        pcolor(reshape(atomCloudNaive ,roiShape)); 
        shading flat; axis off; title('Original');
        cAx = get(gca,'CLim');
        
        subaxis(1,2,2,...
            'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
            'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
            'Margin'        ,   analyVar.subAxis.margin,...
            'MarginTop'     ,   analyVar.subAxis.marginTop,...
            'Padding'       ,   analyVar.subAxis.padding);
        pcolor(reshape(BackCloudApproxState, roiShape)); 
        shading flat; axis off; title('Reconstruction');
        set(gca,'CLim',cAx)