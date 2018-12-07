function [indivCoeffMtx,indivFitStruct,indivWaveNum] = multiLorzFit(waveNum, atomNum, numLorz)
% Fits multiple Lorentzian loss peaks, extracts fit coefficients, for one
% scan.
% 
%INPUTS: Cell arrays of the data to fit; wavenumber readings and the
%corresponding atom number. Each element in the cell array is a seperate
%scan. numLorz is the number of expected loss peaks to be fit.
%
%
%OUTPUTS: Two matrices, one for fit coefficients, and one for the fits
%themselves. Rows correspond to an individual Lorentzian within the scan,
%columns are the fit coefficients or fir data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scanSpan = abs(waveNum(1)-waveNum(length(waveNum)));
lorzRng  = scanSpan/numLorz*0.5;

%% Primary User Inputs. Update row of init guess and bounds when number of expected peaks changes
GammaGuess = -0.003;     % Guess for the FWHM of the peaks
p2 = [-0.660 -0.626 -0.586 -0.545];% Guesses for individual peak center positions from left to right.

for i=1:numLorz % Split data to points around each individual Lorentzian
indivWaveNum(i) = {waveNum(and(waveNum>(p2(i)-lorzRng),waveNum<(p2(i)+lorzRng)))};
indivAtomNum(i) = {atomNum(and(waveNum>(p2(i)-lorzRng),waveNum<(p2(i)+lorzRng)))};
end

c1  = mean(atomNum(and(waveNum>(p2(1)-lorzRng),waveNum<(p2(1)+lorzRng)))); % Take mean of points around each expected Lorz to guess BG.
c2  = mean(atomNum(and(waveNum>(p2(2)-lorzRng),waveNum<(p2(2)+lorzRng))));
c3  = mean(atomNum(and(waveNum>(p2(3)-lorzRng),waveNum<(p2(3)+lorzRng))));
c4  = mean(atomNum(and(waveNum>(p2(4)-lorzRng),waveNum<(p2(4)+lorzRng))));

% Upper and lower bounds for each parameter = [LB1 LB2 LB3 LB4; UB1 UB2 UB3 UB4].
usebnds = 1; % Flag for choosing whether or not to use parameter bounds. 
if usebnds
    bnds(1:2,:) =[-inf -inf -inf -inf;inf inf inf inf];
    bnds(3:4,:) =[-inf -inf -inf -inf;inf inf inf inf];
    bnds(5:6,:) =[-inf -inf -inf -inf;inf inf inf inf];
    bnds(7:8,:) =[-inf -inf -inf -inf;inf inf inf inf];
else
    for i=1:numLorz
        bnds((2*i-1):(2*i),:) = [-inf,-inf,-inf,-inf;inf,inf,inf,inf];
    end
end

% Initial guesses for the Lorentzians. [P1 P2 P3 C]
p1 =  GammaGuess*2*10^2;   % Guess a common amplitude parameter
p3 = (GammaGuess/2)^2;    % Guess a common width parameter
p0(1,:) = [p1 p2(1) p3 c1];
p0(2,:) = [p1 p2(2) p3 c2];
p0(3,:) = [p1 p2(3) p3 c3];
p0(4,:) = [p1 p2(4) p3 c4];

%% Fitting
for i=1:numLorz
    [fit, params, resnorm, residual] = lorentzfit(indivWaveNum{i},indivAtomNum{i},p0(i,:),bnds((2*i-1):(2*i),:));
    output(i).fit      = fit;
    output(i).params   = [params sqrt(params(3))*2];% Calculate a fit Gamma as 5th coeff, based on the width parameter P3.
    output(i).resnorm  = resnorm;
    output(i).residual = residual;
    indivCoeffMtx(i,:) = output(i).params;
    indivFitStruct{i}  = output(i).fit;
end
end



