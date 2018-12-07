function [lineCenters, fullWidths] = one_body_spectra3P1Summed(indVarCell,partNumCell,labelVec,numPeaks)
% Fits the THREE PEAKS in number spectrum (m_j = {-1,0,1}) to a one body loss equation.
% Note: Width ("FWHM") is currently hard-coded.
% Lineshape function - Lorentzian with broadening term
%
% INPUTS:
%   indVarCell  - Cell of vectors containing the independent variables
%   partNumCell - Cell of vectors containing the particle number that is
%                 the specctrum to be analyzed
%   labelVec    - Vector of identifying labels for each scan that is analyzed
%   numPeaks    - Number of m_j peaks that the funciton should expect.
%                 Flagged in the master batch file column 4.
%
% OUTPUTS:
%   amplitude   - Vector the length of numBasenamesAtom containing the fit
%                 amplitudes and standard error
%   lineCenter  - Vector the length of numBasenamesAtom containing the fit
%                 line center positions (Hz) and standard error
%   fullWidth   - Vector the length of numBasenamesAtom containing the fit
%                 full widths of the spectra (Hz) and standard error



%% Initialize loop variables
[lineCenters, fullWidths] = deal(zeros(length(indVarCell),2));

%% Define Plot parameters
specFitFig = figure;
[subPlotRows,subPlotCols] = optiSubPlotNum(length(indVarCell));

%%
% Cases are the number of m_j peaks (numPeaks = {1,2,3})
%------------------------------------------------------------------------
switch numPeaks
    
    case 3
        % Define Physical Functions Used in Calculations
        % coeffs has elements coeffs = [amplitude, line center, halfwidth, num. offset]
        %specFit = @(coeffs,x) coeffs(4)*(1 - coeffs(1)/coeffs(3)*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2))); % Small time approx.
        %specFit  = @(coeffs,x)
        %initOffset*exp(-coeffs(1)/(2*pi*coeffs(3))*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2)))
        specFit   = @(coeffs,x) coeffs(10)*(exp(-coeffs(1)/(2*pi*coeffs(7))*exp(-(x-coeffs(4)).^2/(2*coeffs(7)^2))...
            + (-coeffs(2)/(2*pi*coeffs(8))*exp(-(x-coeffs(5)).^2/(2*coeffs(8)^2)))...
            + (-coeffs(3)/(2*pi*coeffs(9))*exp(-(x-coeffs(6)).^2/(2*coeffs(9)^2)))));
        
        % Loop through each batch file or average scan value
        for iterVar = 1:length(indVarCell) % Indexed by iterVar
            % Reference variables in structure by shorter names for convenience
            indVar = nan(1,length(indVarCell{iterVar}));  indVar(:) = indVarCell{iterVar}*1e6; %The nan function initializes the array with nan values.
            totNum = nan(1,length(partNumCell{iterVar})); totNum(:) = partNumCell{iterVar};
            
            % Initial Guesses
            initOffset = mean(totNum([1:5 end-5:end])); % Mean of first 5 and last 5 points, common to all peaks.
            
            [NumMtx(:,1),NumMtx(:,2)] = findpeaks(-totNum,indVar);%Center Position Guesses are next few rows
            NumMtx(:,1)   = -NumMtx(:,1);
            incNum        = sortrows(NumMtx); % Increasing in number
            ordered       = sortrows(incNum(1:3,:),2); % Picks off expected number of peaks in order in increasing mj
            initCenter    = ordered(:,1);
            initCenterPos = ordered(:,2);
            %initCenterPos(3) = 82.25e6 ;%%%%%% TEMPORARY, hand guessing the
            %initial center position of some lines when the fit clearly fails
            
            FWHM          = 0.03e6; % Intuitive width, guessed by eye. The initial guess is applied to all peaks.
            initWidth     = FWHM/(2*sqrt(2*log(2))); % Fit parameter sigma
            
            intuitiveAmp  = fliplr([(initOffset - initCenter(1)), (initOffset - initCenter(2)), (initOffset - initCenter(3))]); % Amplitude as defined as the size of the peaks
            initAmp       = -2*pi*initWidth/2.*log(intuitiveAmp./initOffset); % Amplitude as defined in the fit function, differernt from the above
            
            % Fitting routine
            specFitModel = NonLinearModel.fit(indVar',totNum,specFit,[initAmp(1) initAmp(2) initAmp(3) initCenterPos(1) initCenterPos(2) initCenterPos(3) [1 1 1]*initWidth initOffset],...
                'CoefficientNames',{'AmplitudeN1','Amplitude0','AmplitudeP1','Line CenterN1','Line Center0','Line CenterP1','SigmaN1','Sigma0','SigmaP1' 'Offset'});
            
            % Calculate output quantities
            lineCenterN1           = double(specFitModel.Coefficients{'Line CenterN1',{'Estimate', 'SE'}}); % m_j=-1
            lineCenter0            = double(specFitModel.Coefficients{'Line Center0',{'Estimate', 'SE'}});  % m_j=0
            lineCenterP1           = double(specFitModel.Coefficients{'Line CenterP1',{'Estimate', 'SE'}}); % m_j=+1
            lineCenters            = horzcat(lineCenterN1,lineCenter0,lineCenterP1);
            
            fullwidthFunc         = @(sigma) sigma.*(2*sqrt(2*log(2))); %Converts fit parameter sigma to meaningfull widths (FWHM)
            fullWidthN1           = [fullwidthFunc(double(specFitModel.Coefficients{'SigmaN1','Estimate'})), fullwidthFunc(double(specFitModel.Coefficients{'SigmaN1','SE'}))];
            fullWidth0            = [fullwidthFunc(double(specFitModel.Coefficients{'Sigma0','Estimate'})), fullwidthFunc(double(specFitModel.Coefficients{'Sigma0','SE'}))];
            fullWidthP1           = [fullwidthFunc(double(specFitModel.Coefficients{'SigmaP1','Estimate'})), fullwidthFunc(double(specFitModel.Coefficients{'SigmaP1','SE'}))];
            fullWidths            = horzcat(fullWidthN1,fullWidth0,fullWidthP1);
            
            fitIndVar = linspace(min(indVar),max(indVar),1e4)'; % Independent variable to plot fit guesses agianst.
            plotcheck = specFit([initAmp(1) initAmp(2) initAmp(3) initCenterPos(1) initCenterPos(2) initCenterPos(3) [1 1 1]*initWidth initOffset],fitIndVar);
            
            % Plot number vs. fit for inspection
            figure(specFitFig); subplot(subPlotRows,subPlotCols,iterVar);
            dataHan   = plot(indVar*1e-6,totNum,fitIndVar*1e-6,specFitModel.predict(fitIndVar));
            
            %%% Plot axis details
            title(num2str(labelVec(iterVar)));
            stringTitle = sprintf('%u (A)',labelVec(iterVar));
            ylabel('Atom Number','FontSize',15,'FontWeight','bold');
            xlabel('Spec Slave Frequency (MHz)','FontSize',15,'FontWeight','bold');
            set(gca,'FontSize',23);
            set(gca,'LineWidth',1);
            grid on; axis tight;
            curYLim = get(gca,'YLim'); set(gca,'YLim',[0 max(curYLim)])
            set(dataHan(1)        ,...
                'LineStyle'       ,   'none'  ,...
                'Marker'          ,   'o'     ,...
                'MarkerSize'      ,   7      ,...
                'MarkerFaceColor' ,   'b'     );
            set(dataHan(2),...
                'LineWidth' ,   2.5   )
            if iterVar == length(indVarCell);
                set(gcf,'Name','Spectra Fits');
            end
            
            
            hold on
            scatter(initCenterPos*10^-6,initCenter,'filled','d','MarkerFaceColor','red') % Plot the datapoints used as peak initial guesses.
            hold on
            plot(fitIndVar*1e-6,plotcheck) % Plot the number function with the initial guesses
            sprintf('\nData Columns Arranged as m_j=-1|SE|m_j=0|SE|m_j=+1|SE')
            legend('Data','Fit', 'Center Init. Guess', 'Fitfunc plot w/init. guesses')
        end
        
        
    case 2
        % Define Physical Functions Used in Calculations
        % coeffs has elements coeffs = [amplitude, line center, halfwidth, num. offset]
        %specFit = @(coeffs,x) coeffs(4)*(1 - coeffs(1)/coeffs(3)*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2))); % Small time approx.
        %specFit  = @(coeffs,x)
        %initOffset*exp(-coeffs(1)/(2*pi*coeffs(3))*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2)))
        specFit   = @(coeffs,x) coeffs(7)*(exp(-coeffs(1)/(2*pi*coeffs(5))*exp(-(x-coeffs(3)).^2/(2*coeffs(5)^2))...
            + (-coeffs(2)/(2*pi*coeffs(6))*exp(-(x-coeffs(4)).^2/(2*coeffs(6)^2)))));
        
        % Loop through each batch file or average scan value
        for iterVar = 1:length(indVarCell) % Indexed by iterVar
            % Reference variables in structure by shorter names for convenience
            indVar = nan(1,length(indVarCell{iterVar}));  indVar(:) = indVarCell{iterVar}*1e6; %The nan function initializes the array with nan values.
            totNum = nan(1,length(partNumCell{iterVar})); totNum(:) = partNumCell{iterVar};
            
            % Initial Guesses
            initOffset = mean(totNum([1:5 end-5:end])); % Mean of first 5 and last 5 points, common to all peaks.
            
            [NumMtx(:,1),NumMtx(:,2)] = findpeaks(-totNum,indVar);%Center Position Guesses are next few rows
            NumMtx(:,1)   = -NumMtx(:,1);
            incNum        = sortrows(NumMtx); % Increasing in number
            ordered       = sortrows(incNum(1:2,:),2); % Picks off expected number of peaks in order in increasing mj
            initCenter    = ordered(:,1);
            initCenterPos = ordered(:,2);
            
            FWHM          = 0.03e6; % Intuitive width, guessed by eye. The initial guess is applied to all peaks.
            initWidth     = FWHM/(2*sqrt(2*log(2))); % Fit parameter sigma
            
            intuitiveAmp  = fliplr([(initOffset - initCenter(1)), (initOffset - initCenter(2))]); % Amplitude as defined as the size of the peaks
            initAmp       = -2*pi*initWidth/2.*log(intuitiveAmp./initOffset); % Amplitude as defined in the fit function, differernt from the above
            
            % Fitting routine
            specFitModel = NonLinearModel.fit(indVar',totNum,specFit,[initAmp(1) initAmp(2) initCenterPos(1) initCenterPos(2) [1 1]*initWidth initOffset],...
                'CoefficientNames',{'AmplitudeN1','AmplitudeP1','Line CenterN1','Line CenterP1','SigmaN1','SigmaP1' 'Offset'});
            
            % Calculate output quantities
            lineCenterN1           = double(specFitModel.Coefficients{'Line CenterN1',{'Estimate', 'SE'}});
            lineCenterP1           = double(specFitModel.Coefficients{'Line CenterP1',{'Estimate', 'SE'}});
            lineCenters            = horzcat(lineCenterN1,lineCenterP1);
            
            fullwidthFunc         = @(sigma) sigma.*(2*sqrt(2*log(2))); %Converts fit parameter sigma to meaningfull widths (FWHM)
            fullWidthN1           = [fullwidthFunc(double(specFitModel.Coefficients{'SigmaN1','Estimate'})), fullwidthFunc(double(specFitModel.Coefficients{'SigmaN1','SE'}))];
            fullWidthP1           = [fullwidthFunc(double(specFitModel.Coefficients{'SigmaP1','Estimate'})), fullwidthFunc(double(specFitModel.Coefficients{'SigmaP1','SE'}))];
            fullWidths            = horzcat(fullWidthN1,fullWidthP1);
            
            fitIndVar = linspace(min(indVar),max(indVar),1e4)'; % Independent variable to plot fit guesses agianst.
            plotcheck = specFit([initAmp(1) initAmp(2) initCenterPos(1) initCenterPos(2) [1 1]*initWidth initOffset],fitIndVar);
            
            % Plot number vs. fit for inspection
            figure(specFitFig); subplot(subPlotRows,subPlotCols,iterVar);
            dataHan   = plot(indVar*1e-6,totNum,fitIndVar*1e-6,specFitModel.predict(fitIndVar));
            
            %%% Plot axis details
            title(num2str(labelVec(iterVar)));
            ylabel('Atom Number','FontSize',15,'FontWeight','bold');
            xlabel('Spec Slave Frequency (MHz)','FontSize',15,'FontWeight','bold');
            set(gca,'FontSize',23)
            set(gca,'LineWidth',1)

            grid on; axis tight;
            curYLim = get(gca,'YLim'); set(gca,'YLim',[0 max(curYLim)])
            set(dataHan(1)        ,...
                'LineStyle'       ,   'none'  ,...
                'Marker'          ,   'o'     ,...
                'MarkerSize'      ,   7      ,...
                'MarkerFaceColor' ,   'b'     );
            set(dataHan(2),...
                'LineWidth' ,   2.5   )
            if iterVar == length(indVarCell);
                set(gcf,'Name','Spectra Fits');
            end
            
            hold on
            scatter(initCenterPos*10^-6,initCenter,'filled','d','MarkerFaceColor','red') % Plot the datapoints used as peak initial guesses.
            hold on
            plot(fitIndVar*1e-6,plotcheck) % Plot the number function with the initial guesses
            sprintf('\nData Columns Arranged as m_j=-1|SE|m_j=+1|SE')
            legend('Data','Fit', 'Center Init. Guess', 'Fitfunc plot w/init. guesses')
        end
        
    case 1
        % Define Physical Functions Used in Calculations
        % coeffs has elements coeffs = [amplitude, line center, halfwidth, num. offset]
        %specFit = @(coeffs,x) coeffs(4)*(1 - coeffs(1)/coeffs(3)*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2))); % Small time approx.
        %specFit  = @(coeffs,x)
        %initOffset*exp(-coeffs(1)/(2*pi*coeffs(3))*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2)))
        specFit   = @(coeffs,x) coeffs(4)*(exp(-coeffs(1)/(2*pi*coeffs(3))*exp(-(x-coeffs(2)).^2/(2*coeffs(3)^2))));
        
        % Loop through each batch file or average scan value
        for iterVar = 1:length(indVarCell) % Indexed by iterVar
            % Reference variables in structure by shorter names for convenience
            indVar = nan(1,length(indVarCell{iterVar}));  indVar(:) = indVarCell{iterVar}*1e6; %The nan function initializes the array with nan values.
            totNum = nan(1,length(partNumCell{iterVar})); totNum(:) = partNumCell{iterVar};
            
            % Initial Guesses
            initOffset = mean(totNum([1:5 end-5:end])); % Mean of first 5 and last 5 points, common to all peaks.
            
            [NumMtx(:,1),NumMtx(:,2)] = findpeaks(-totNum,indVar);%Center Position Guesses are next few rows
            NumMtx(:,1)   = -NumMtx(:,1);
            incNum        = sortrows(NumMtx); % Increasing in number
            ordered       = sortrows(incNum(1:1,:),2); % Picks off expected number of peaks in order in increasing mj
            initCenter    = ordered(:,1);
            initCenterPos = ordered(:,2);
            
            FWHM          = 0.03e6; % Intuitive width, guessed by eye. The initial guess is applied to all peaks.
            initWidth     = FWHM/(2*sqrt(2*log(2))); % Fit parameter sigma
            
            intuitiveAmp  = (initOffset - initCenter(1)); % Amplitude as defined as the size of the peaks
            initAmp       = -2*pi*initWidth/2.*log(intuitiveAmp./initOffset); % Amplitude as defined in the fit function, differernt from the above
            
            % Fitting routine
            specFitModel = NonLinearModel.fit(indVar',totNum,specFit,[initAmp initCenterPos initWidth initOffset],...
                'CoefficientNames',{'Amplitude','Line Center','Sigma' 'Offset'});
            
            % Calculate output quantities
            lineCenters           = double(specFitModel.Coefficients{'Line Center',{'Estimate', 'SE'}});
            
            fullwidthFunc         = @(sigma) sigma.*(2*sqrt(2*log(2))); %Converts fit parameter sigma to meaningfull widths (FWHM)
            fullWidths            = [fullwidthFunc(double(specFitModel.Coefficients{'Sigma','Estimate'})), fullwidthFunc(double(specFitModel.Coefficients{'Sigma','SE'}))];
            
            fitIndVar = linspace(min(indVar),max(indVar),1e4)'; % Independent variable to plot fit guesses agianst.
            plotcheck = specFit([initAmp initCenterPos initWidth initOffset],fitIndVar);
           

            
            % Plot number vs. fit for inspection
            figure(specFitFig); subplot(subPlotRows,subPlotCols,iterVar);
            dataHan   = plot(indVar*1e-6,totNum,fitIndVar*1e-6,specFitModel.predict(fitIndVar));
            
                        %%% Plot axis details
            title(num2str(labelVec(iterVar)));
            stringTitle = sprintf('%u (A)',labelVec(iterVar));
            ylabel('Atom Number','FontSize',15,'FontWeight','bold');
            xlabel('Spec Slave Frequency (MHz)','FontSize',15,'FontWeight','bold');
            set(gca,'FontSize',23);
            set(gca,'LineWidth',1);
            grid on; axis tight;
            curYLim = get(gca,'YLim'); set(gca,'YLim',[0 max(curYLim)])
            set(dataHan(1)        ,...
                'LineStyle'       ,   'none'  ,...
                'Marker'          ,   'o'     ,...
                'MarkerSize'      ,   7      ,...
                'MarkerFaceColor' ,   'b'     );
            set(dataHan(2),...
                'LineWidth' ,   2.5   )
            if iterVar == length(indVarCell);
                set(gcf,'Name','Spectra Fits');
            end
            legend('Data','Fit', 'Center Init. Guess', 'Fitfunc plot w/init. guesses')
                        
            hold on
            scatter(initCenterPos*10^-6,initCenter,'filled','d','MarkerFaceColor','red') % Plot the datapoints used as peak initial guesses.
            hold on
            plot(fitIndVar*1e-6,plotcheck) % Plot the number function with the initial guesses
            legend('Data','Fit', 'Center Init. Guess', 'Fitfunc plot w/init. guesses')
            sprintf('\nData Columns Arranged as m_j=0')
        end
        
    otherwise
        warning('Number of peaks is not 1,2,or 3.')
        
end



lineCenters./10^6
fullWidths./10^3

end

