function [plotTime,fitGndState,fitExtState,rabiFreq,decohFreq,delta,totNum] = fitRabiOsc(measTime,atomPop,rabiParams,fitConstr)
%Function to fit rabi oscillations by solving thw two level optical bloch equations with damping
%
% INPUTS:
%   time        - vector of time points where Rabi oscillations measurements were obtained
%   atomPop     - vector of Rabi oscillation ground state population data
%                 normalized to total population (i.e. max population = 1)
%
% OUTPUTS:
%   fitTime     - 
%   fitGndState -
%   rabiFreq    -
%   decohFreq   -

%% Find inital guess for Rabi frequency
% This analysis finds the power spectrum and assumes a constant (or nearly constant) 
% sampling rate
[Pxx,f]   = periodogram(atomPop - mean(atomPop),[],[],1/(mean(diff(measTime))));
[~,loc]   = max(Pxx);
omegaRabi = 2*pi*f(loc);

%% Adjust data set for fitting
fitTime = [0; measTime + rabiParams.pulseOffset];

switch rabiParams.state
    case 'Ground' % If the ground state is measured
        fitAtomPop = [rabiParams.initRho(1); atomPop];
        rhoInd     = 1;
    case 'Excited' % If the excited state is measured
        fitAtomPop = [rabiParams.initRho(2); atomPop];
        rhoInd     = 2;
    otherwise
        error('Not a valid state selection.');
end

%% Fitting Routine
% Compile guess parameters into vector and make anonymous function for
% nonlinear fitting (this makes the fit caller into the normal form Matlab expects)
switch fitConstr
    case 'fitNum' %this case assumes you know the detuning of the laser
        guessParams = [omegaRabi, rabiParams.gammaDec, rabiParams.totNumMeas];
        fitCaller   = @(coeffs,time) fitOBEFunc(time,rabiParams,coeffs(1),coeffs(2),rabiParams.delta,coeffs(3),rhoInd);
        
    case 'fitDelta' %this case assumes you know the background atom number
        guessParams = [omegaRabi, rabiParams.gammaDec, rabiParams.delta];
        fitCaller   = @(coeffs,time) fitOBEFunc(time,rabiParams,coeffs(1),coeffs(2),coeffs(3),rabiParams.totNumMeas,rhoInd);
end

%Do the non-linear fitting
%oscCoeffsModel = NonLinearModel.fit(fitTime,fitAtomPop,fitCaller,guessParams,'CoefficientNames',coeffNames)
[oscCoeffsModel] = lsqcurvefit(fitCaller,guessParams,fitTime,fitAtomPop);
                                
% Save coefficients into variables
rabiFreq  = oscCoeffsModel(1); % [s^-1] Rabi rate
decohFreq = oscCoeffsModel(2); % [s^-1] decoherence rate
switch fitConstr
    case 'fitNum' %this case assumes you know the detuning of the laser
        totNum = oscCoeffsModel(3);
        delta  = rabiParams.delta; % [s^-1] Detuning       
    case 'fitDelta' %this case assumes you know the background atom number
        totNum = rabiParams.totNumMeas;
        delta  = oscCoeffsModel(3); % [s^-1] Detuning
end

%% Get fit oscillations from OBE
plotTime    = linspace(0,max(fitTime),1e3)';
fitGndState = fitOBEFunc(plotTime,rabiParams,rabiFreq,decohFreq,delta,totNum,1);
fitExtState = fitOBEFunc(plotTime,rabiParams,rabiFreq,decohFreq,delta,totNum,2);

end

function output = fitOBEFunc(tRange,rabiParams,omegaRabi,gammaCoh,delta,num,rhoInd)    
    [~,rho] = ode45(@(t,rho) funcOBE(t,rho,omegaRabi,rabiParams.gammaLife,gammaCoh,delta),tRange,rabiParams.initRho);
    output  = num*rho(:,rhoInd); % Solve for ground state atom population
end