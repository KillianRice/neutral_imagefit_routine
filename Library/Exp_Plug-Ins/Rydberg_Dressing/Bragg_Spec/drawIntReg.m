function drawIntReg(analyVar,indivDataset,basenameNum,iPkSub,jPkSub,numPeak)
% Draws an ellipse showing the integration region used for numerical integration

%% Choose where to draw the ellipses
%%%%%%%-----------------------------------%%%%%%%%%%
figure(analyVar.figNum.atomEvol + basenameNum)
figChild = sort(get(gcf,'Children'));

for k = 1:indivDataset{basenameNum}.CounterAtom;
    %% Get fit center in full image and find one sigma bounds
    %%%%%%%-----------------------------------%%%%%%%%%%
    CloudCntr  = round((analyVar.roiWinRadAtom(basenameNum) - (analyVar.funcFitWin(basenameNum) - 1)/2)...
            + [indivDataset{basenameNum}.All_fitParams{k}{1}.xCntr indivDataset{basenameNum}.All_fitParams{k}{1}.yCntr]);
    oneSigBnds = [round(CloudCntr - cellfun(@(x) indivDataset{basenameNum}.All_fitParams{k}{1}.(x),{'sigX' 'sigY'})); ...
                  round(CloudCntr + cellfun(@(x) indivDataset{basenameNum}.All_fitParams{k}{1}.(x),{'sigX' 'sigY'}))];
              
    %% Set radii for the main cloud and Bragg peak
    %%%%%%%-----------------------------------%%%%%%%%%%
    rBragg = analyVar.braggWinRad;
    aBEC = indivDataset{basenameNum}.All_fitParams{k}{1}.sigX;
    bBEC = indivDataset{basenameNum}.All_fitParams{k}{1}.sigY;
    
    %% Draw ellipses on the false color image
    %%%%%%%-----------------------------------%%%%%%%%%%
    set(gcf,'CurrentAxes',figChild(k))
    % Loop through peaks
    for braggIter = 1:length(numPeak)
        x = oneSigBnds(1,1) + jPkSub{braggIter}(k) - 1;
        if analyVar.winID(basenameNum) == 1 || (analyVar.winID(basenameNum) == 2 && numPeak(braggIter) == 1)
            y = iPkSub{braggIter}(k);
        end
        
        if analyVar.winID(basenameNum) == -1 || (analyVar.winID(basenameNum) == 2 && numPeak(braggIter) == -1)
            y = iPkSub{braggIter}(k) + CloudCntr(2);
        end
        
        % Draw Bragg area
        rectangle('Position',[x-rBragg y-rBragg [2 2]*rBragg],'Curvature',[1 1],'EdgeColor','r','LineWidth',2)
    end
    % Draw condensate area
    rectangle('Position',[CloudCntr(1)-aBEC CloudCntr(2)-bBEC 2*aBEC 2*bBEC],'Curvature',[1 1],'EdgeColor','r','LineWidth',2)
end
end