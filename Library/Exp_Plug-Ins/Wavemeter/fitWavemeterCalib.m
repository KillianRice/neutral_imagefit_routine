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

%Set any readings flagged as 1 (or spuriously low), to whatever number came before it.
for i=1:length(wavemeter)
    if wavemeter(i) < 21690
        wavemeter(i) = wavemeter(i-1);
    end
end

% Calculate and display the fit coefficients
format long
indivBatch.indCalib = polyfit(ind_var,wavemeter,2);
indivBatch.indCalib
format short

% Save old x-axis into a new variable so we don't lose it entirely
indivBatch.rawimagevcoAtom = indivBatch.imagevcoAtom;

% Convert indipendent variable to wavenumbers 
indivBatch.imagevcoAtom = polyval(indivBatch.indCalib, ind_var);
end