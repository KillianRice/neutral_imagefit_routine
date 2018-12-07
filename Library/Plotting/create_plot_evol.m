function evolAxH = create_plot_evol(analyVar,indivDataset)
% Function to plot the OD images of a scan on a subplot.
% This function does no parameter extraction from the fit but serves only
% to show the evolution of the cloud through a scan.
%
% INPUTS:
%   analyVar     - structure of all pertinent variables for the imagefit
%                  routines
%   indivDataset - Cell of structures containing all scan/batch
%                  specific data
%
% OUTPUTS:
%   evolAxH - Cell of vectors containing handles to the subplot axes. 
%             Used for setting standard color limits on pcolor plots.
%

%% Preallocate loop variables
evolAxH = cell(1,analyVar.numBasenamesAtom);

%% Loop through each batch file listed in basenamevectorAtom
for basenameNum = 1:analyVar.numBasenamesAtom
    %% Preallocate loop variables
    evolAxH{basenameNum} = zeros(1,indivDataset{basenameNum}.CounterAtom);
    
        % Processes all the image files in this batch/scan
       figure(analyVar.figNum.atomEvol + basenameNum); 
        for k = 1:indivDataset{basenameNum}.CounterAtom;
            %% Plot Evolution
            evolAxH{basenameNum}(k) = subaxis(indivDataset{basenameNum}.SubPlotRows,indivDataset{basenameNum}.SubPlotCols,k,...
                'SpacingVert'   ,   analyVar.subAxis.spacingVert,...
                'SpacingHoriz'  ,   analyVar.subAxis.spacingHoriz,...
                'Margin'        ,   analyVar.subAxis.margin,...
                'MarginTop'     ,   analyVar.subAxis.marginTop,...
                'Padding'       ,   analyVar.subAxis.padding);
            pcolor(indivDataset{basenameNum}.All_OD_Image{k}); shading flat;
            %%% Plot axis details
            title(strcat(num2str(indivDataset{basenameNum}.imagevcoAtom(k)),[' ' analyVar.xDataUnit]));
			 hold on; grid off; axis off
            if k == indivDataset{basenameNum}.CounterAtom;
                set(gcf,'Name',['Cloud Evolution: Time = ' num2str(analyVar.timevectorAtom(basenameNum))]);
                mtit(['Cloud Evolution: Time = ' num2str(analyVar.timevectorAtom(basenameNum))],'FontSize',16,'zoff',.05,'xoff',-.01)
            end
        end
        set(analyVar.figNum.atomEvol + basenameNum,'color','white')
end