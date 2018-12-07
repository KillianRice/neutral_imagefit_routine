function coeffOut = fit_bragg_spectra(analyVar,indVar,totNum,label,expTime,figSpec,specFitFig,iterVar,maxIter)
% Fits the number spectrum to a one body loss equation
% Derivation included in folder
% Lineshape function - Lorentzian with broadening term
%
% INPUTS:
%   indVarCell  - Cell of vectors containing the independent variables
%   partNumCell - Cell of vectors containing the particle number that is
%                 the specctrum to be analyzed
%   labelVec    - Vector of identifying labels for each scan that is analyzed
%   expTime     - Vector of exposure times for each dataset
%
% OUTPUTS:
%   coeffOut    - Structure of fit coefficients 

%% Got subplot area
[subPlotRows,subPlotCols] = optiSubPlotNum(maxIter);

%% Initialize loop variables
% Define the coefficient fields that will be saved to the output structure
fitFields = {{'Rabi_Freq' 'Line_Center'} {'Amplitude' 'Line_Center' 'FWHM'}};

% Reference variables in structure by shorter names for convenience
origIVar = linspace(min(indVar),max(indVar),1e4)'; %ind. variable before modification

switch analyVar.fitFuncBragg
    case 'Sinc'
        %% Define Physical Functions Used in Calculations
        % coeffs has elements coeffs = [rabi_freq, line center, offset]
        % Eq. 7.27 from Foot
        specFit = @(coeffs,x) coeffs(1)^2./(coeffs(1)^2+((2*pi)*(x-coeffs(2))).^2).*(sin(sqrt(coeffs(1)^2+((2*pi)*(x-coeffs(2))).^2)*expTime/2)).^2;
        
        %% Initial Guesses
        smoothNum = smooth(totNum,'sgolay',1); % smooth for helping make guesses
        [maxVal, minLoc] = max(smoothNum); % find location of minimum
        
        %initOffset = mean(totNum([1:5 end-5 end])); % Mean of first 5 and last 5 points
        initCenter = indVar(minLoc); % Location of minimum in smoothed data
        peakLoss   = maxVal;             % Value below guessed offset
        
        % Function to estimate rabi frequency based on peak loss
        initRabi = 2*asin(sqrt(peakLoss))/expTime;
        
        %% Fitting routine
        specFitModel = NonLinearModel.fit(indVar',totNum,specFit,[initRabi initCenter],...
            'CoefficientNames',{'Rabi Frequency','Line Center'});
        
        % Calculate output quantities
        for i = 1:length(fitFields{1})
            coeffOut.(fitFields{1}{i}) = specFitModel.Coefficients.Estimate(i);
        end
        
        %% Plot number vs. fit for inspection
        figure(specFitFig); subplot(subPlotRows,subPlotCols,iterVar);
        fitIndVar = linspace(min(indVar),max(indVar),1e4)';
        dataHan   = plot(indVar,totNum,fitIndVar,specFitModel.predict(fitIndVar));
        %%% Plot axis details
        title(num2str(label(iterVar)));
        xlabel('kHz'); grid on; axis tight; hold on
        set(dataHan(1),'LineStyle','none','Marker','*'); set(dataHan(2),'LineWidth',2)
        if iterVar == length(indVar);
            set(gcf,'Name','Bragg Spectra Fits');
        end
        
    case 'Gaussian'
        %% Gaussian Case
        specFit = @(coeffs,x) coeffs(4) + coeffs(1)*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2));
        
        % Initial Guesses
        smoothNum = smooth(totNum,'sgolay',1); % smooth for helping make guesses
        [maxVal, minLoc] = max(smoothNum); % find location of minimum
        
        initOffset = mean(totNum([1:5 end-5 end])); % Mean of first 5 and last 5 points
        initCenter = indVar(minLoc);                % Location of minimum in smoothed data
        initWidth = sqrt(abs(sum((indVar'-initCenter).^2.*totNum)./sum(totNum)));
        initAmp    = (initOffset + maxVal);         % Value below guessed offset
        
        % Fitting routine
        specFitModel = NonLinearModel.fit(indVar',totNum,specFit,[initAmp initCenter initWidth initOffset],...
            'CoefficientNames',{'Amplitude','Line Center','Sigma','Number Offset'});
        
        % Calculate output quantities
        for i = 1:length(fitFields{2})
            coeffOut.(fitFields{2}{i}) = specFitModel.Coefficients.Estimate(i);
            
            if strcmpi(fitFields{2}{i},'FWHM')
                coeffOut.(fitFields{2}{i}) = 2*sqrt(2*log(2))*specFitModel.Coefficients.Estimate(i);
            end
        end
        
        % Calculate output errors
        for i = 1:length(fitFields{2})
            coeffOut.([fitFields{2}{i} '_StanErr']) = specFitModel.Coefficients.SE(i);
            
            if strcmpi(fitFields{2}{i},'FWHM')
                coeffOut.([fitFields{2}{i} '_StanErr']) = 2*sqrt(2*log(2))*specFitModel.Coefficients.SE(i);
            end
        end
        
        %% Plot number vs. fit for inspection
        figure(specFitFig); subplot(subPlotRows,subPlotCols,iterVar);
        fitIndVar = linspace(min(indVar),max(indVar),1e4)';
        dataHan   = plot(indVar,totNum,fitIndVar,specFitModel.predict(fitIndVar));
        %%% Plot axis details
        title(num2str(label));
        xlabel('Detuning [kHz]'); grid on; axis tight; hold on
        set(dataHan(1),'LineStyle','none','Marker','*'); set(dataHan(2),'LineWidth',2)
        if iterVar == maxIter;
            set(gcf,'Name','Spectra Fits');
        end
end

%% Plot fit values on the original plot
figure(figSpec);
figChld = sort(get(figSpec,'Children')); set(figSpec,'CurrentAxes',figChld(1)); hold on
fitHan = plot(origIVar,specFit(specFitModel.Coefficients.Estimate,fitIndVar));
set(fitHan,'LineWidth',3,'LineStyle',':','Color',analyVar.COLORS(iterVar ,:)); axis tight
end