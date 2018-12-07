function PCA_constructBasis(flagSubDCComp)
% Want to be able to see how PCA is doing constructing the basis of the
% background images. Going to reload all of the images and get the
% background part of the ROI (avoiding splitting out the atom and not atom
% parts)
%
% INPUTS
%   flagSubDCComp - will toggle between using the traditional way of 
%                   background analysis or the new trial version

%% Initialize variables
analyVar = AnalysisVariables;
indivDataset = get_indiv_batch_data(analyVar);

basenameNum = 1;

roiWinSize = [1 1].*(2*analyVar.roiWinRadAtom(basenameNum) + 1);
optimParmas = optimoptions('fminunc','Display','notify-detailed','PlotFcns',{@optimplotx,@optimplotfval,@optimplotstepsize},'Algorithm','quasi-newton');

%% Load backgrounds
BackgroundROI    = zeros(prod(roiWinSize),indivDataset{basenameNum}.CounterBack);
for k = 1:indivDataset{basenameNum}.CounterBack
    t = [analyVar.dataDir char(indivDataset{basenameNum}.fileBack(k)) analyVar.dataBack]; tFID = fopen(t,'r','ieee-be');
    fullRawImageBack = double(fread(tFID,analyVar.matrixSize,'*int16')); fclose(tFID);
    
    % Separate Background from entire image
    % This matrix will include all the background images from the scan
    BackgroundROI(:,k)    = fullRawImageBack(indivDataset{basenameNum}.image_Index ~= 0);
end

BackgroundROI = bsxfun(@minus,BackgroundROI, mean(BackgroundROI,2));

% Pick out the first background as our test background
% This is the background image we want to try and match (typically the atom background)
testImg = BackgroundROI(:,1);

%% Choose whether you want to try and subtract off the DC component or not
if flagSubDCComp
    % partner is the nearest background image in time for subtracting DC component
    partnerBack   = BackgroundROI(:,2);
    
    % Backgrounds to use for matching are everything except the partner
    BackgroundROI = bsxfun(@minus,BackgroundROI(:,3:end),partnerBack);
    
    % The image to match is the test image with the DC components removed
    imgToMatch    = testImg - partnerBack; 
    imgToMatch    = imgToMatch - mean(BackgroundROI,2);
else
    % Backgrounds to use for matching are all available images
    BackgroundROI = BackgroundROI(:,2:end);
    
    % The image to match is the raw test image
    imgToMatch    = testImg;
    %imgToMatch    = imgToMatch - mean(imgToMatch);
end

%% Construct principal components
% Compute principal components from background images
[pcCoeffs, pcBasis, pcEigenVals] = pca(BackgroundROI);

%% Curious to see if I can reconstruct the first background image
% take the inital conditions of the nth image as the projection onto the pcBasis
InitialCondition = pcCoeffs(1,1:15);
%InitialCondition = sum(bsxfun(@times,imgToMatch,pcBasis))/(sum(sum(bsxfun(@times,imgToMatch,pcBasis))))
%InitialCondition = imgToMatch'*pcBasis;

% Minimize the nth state in the original basis to the new pcBasis
A  = fminunc(@(A) WeightedBackgroundFunction( A,...
                  imgToMatch,...
                  pcBasis(:,1:15)),InitialCondition,optimParmas)

% Coefficients define how to transform original basis into pcBasis
% Need BackCloud in PCA basis to construct the nth image background in
% terms of the PCA basis of the cloud background
%pcBGCloud = (pcCoeffs'*BackgroundROI')';
% Construct linear approximation of cloud background using PCA basis of not cloud backgrounds
backApproxState = sum(bsxfun(@times,pcBasis(:,1:15),A),2);

%% Plot images
rawBackFig = figure;
pcaBackFig = figure;
compFig    = figure;

for k = 1:size(BackgroundROI,2);
% Plot the raw backgrounds
    figure(rawBackFig);
    subaxis(indivDataset{basenameNum}.SubPlotRows,indivDataset{basenameNum}.SubPlotCols,k,...
                'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
                'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
                'Margin'        ,   analyVar.subAxis.margin,...
                'MarginTop'     ,   analyVar.subAxis.marginTop,...
                'Padding'       ,   analyVar.subAxis.padding);
    pcolor(reshape(BackgroundROI(:,k),roiWinSize)); shading flat; axis off
    title(sprintf('%g',k));

% Plot the principal components
    figure(pcaBackFig);
    subaxis(indivDataset{basenameNum}.SubPlotRows,indivDataset{basenameNum}.SubPlotCols,k,...
                'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
                'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
                'Margin'        ,   analyVar.subAxis.margin,...
                'MarginTop'     ,   analyVar.subAxis.marginTop,...
                'Padding'       ,   analyVar.subAxis.padding);
    pcolor(reshape(pcBasis(:,k),roiWinSize)); shading flat; axis off
    title(sprintf('%g',pcEigenVals(k)/sum(pcEigenVals)));

end
% Compare the first background, its reconstruction, and their difference
    figure(compFig);
    subaxis(1,3,1,...
                'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
                'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
                'Margin'        ,   analyVar.subAxis.margin,...
                'MarginTop'     ,   analyVar.subAxis.marginTop,...
                'Padding'       ,   analyVar.subAxis.padding);
    pcolor(reshape(imgToMatch,roiWinSize)); shading flat; axis off
    title('Original');
    
    subaxis(1,3,2,...
                'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
                'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
                'Margin'        ,   analyVar.subAxis.margin,...
                'MarginTop'     ,   analyVar.subAxis.marginTop,...
                'Padding'       ,   analyVar.subAxis.padding);
    pcolor(reshape(backApproxState,roiWinSize)); shading flat; axis off
    title('Reconstruction');
    
    subaxis(1,3,3,...
                'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
                'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
                'Margin'        ,   analyVar.subAxis.margin,...
                'MarginTop'     ,   analyVar.subAxis.marginTop,...
                'Padding'       ,   analyVar.subAxis.padding);
    pcolor(reshape(imgToMatch - backApproxState,roiWinSize)); shading flat; axis off
    title(sprintf('Sum Mean Squared Difference: %g',sum((imgToMatch - backApproxState).^2/length(imgToMatch))));
end