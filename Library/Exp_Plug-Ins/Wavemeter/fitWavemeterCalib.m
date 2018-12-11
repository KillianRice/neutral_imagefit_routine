function indivBatch = fitWavemeterCalib(analyVar, indivBatch)
%PARAM_EXTRACT_IND_VAR_TO_WAVEMETER: Flagged on by "useWavemeter" variable
%in AnalysisVariables, and called at line ~119 in get_indiv_batch_data.
%
% Reads in wavemeter measurements from the last column of each batch file,
% interpolates between the points and generages a set of coefficeients per
% scan. Plots the results and changes the indipendent vatiable of the batch
% from voltage to wavenumbers via the interpolation coefficients for that
% scan. 
%
% To change the label for the horizontal axis, change "TimeOrDetune" to
% "Wavemeter", in analysis variables.
%
%STATUS: 2018.08.22: Wavemeter plugin does not work with averaging or normalization of atom number, so turn off the
%"normMeanNum" variable in AnalysisVariables when using. 
%
%
ind_var = indivBatch.imagevcoAtom;
wavemeter = indivBatch.sigParamAtom;

%Remove any readings flagged as 1 (or spuriously low)
correctWavemeterReadings = wavemeter  > 21690;
ind_var   = ind_var(correctWavemeterReadings);
wavemeter = wavemeter(correctWavemeterReadings);

% Calculate and display the fit coefficients
format long
indivBatch.indCalib = polyfit(ind_var,wavemeter,2);
indivBatch.indCalib
format short

% Save old x-axis into a new variable so we don't lose it entirely
indivBatch.rawimagevcoAtom = indivBatch.imagevcoAtom;

% Convert independent variable to wavenumbers 
indivBatch.imagevcoAtom = polyval(indivBatch.indCalib, indivBatch.imagevcoAtom);
end