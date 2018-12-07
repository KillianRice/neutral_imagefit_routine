function accelFit = Pixel_Size_Fit(time,pos,accel,pixelSize)
% Function to fit and extract the pixel size from drop data
% NOTE: This function will assume that the entire image is used for analysis and that the cloud center was
% returned as an absolute position. (code written to analyze data from 2013_Rydberg_Dressing\2015.03.20)

%% Fitting function
accelFunc = @(coeffs,t) coeffs(2) + coeffs(1).*(1/2*accel*t.^2);

%% Guessing
initPos        = pos(1);
pixelPerMicron = 1/pixelSize;

%% Fitting
accelFit = NonLinearModel.fit(time,pos,accelFunc,[pixelPerMicron initPos],'CoefficientNames',{'Pixel_Size' 'Initial_Position'});

%% Print and Plot results
fprintf('\nFit pixel size is %g microns\n\n',1/accelFit.Coefficients.Estimate(1)*1e6)