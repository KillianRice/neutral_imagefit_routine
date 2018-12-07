function intStruct = loadIntenData(lcl_analyVar,lcl_indivDataset)
% Function to load intensity data from Picoscope 

%% Loop through each batch file and all images
%%%%%%%-----------------------------------%%%%%%%%%%

% Preallocate nested loop variables
intStruct   = cell(lcl_analyVar.numBasenamesAtom,1);

for basenameNum = 1:lcl_analyVar.numBasenamesAtom

% Preallocate nested loop variables   
intenData.intVec = cell(lcl_indivDataset{basenameNum}.CounterAtom,1);
intenData.p2p    = zeros(lcl_indivDataset{basenameNum}.CounterAtom,1);

    % Processes all the image files in the current batch
    for k = 1:lcl_indivDataset{basenameNum}.CounterAtom;
       % Construct filename and read in data
       intenFilename = [lcl_indivDataset{basenameNum}.fileAtom{k} '_' lcl_analyVar.molePicoChan '.pico'];
       intenDataFile = [lcl_analyVar.dataDir intenFilename];
       tmpIntVec     = dlmread(intenDataFile);
       
       % If the scanned variable is detuning then read in frequency from
       % the ind. variable. Otherwise, perform an FFT to estimate the freq.
       switch lcl_analyVar.TimeOrDetune
           case 'Detuning'
               % Retrieve the oscillatino frequency of the beat note
                freq = lcl_analyVar.funcDataScale(lcl_indivDataset{basenameNum}.imagevcoAtom(k))*1e3;
                
           otherwise
               % Compute power spectral density to estimate the frequency of oscillation
               [Pxx,f] = periodogram(tmpIntVec(2,:) - mean(tmpIntVec(2,:)),[],[],lcl_analyVar.fSamp);
               [~,loc] = max(Pxx);
               freq    = f(loc);
       end
       
       % Calculate the mean peak to peak oscillation amplitude       
       [~,maxPntsLoc] = findpeaks(tmpIntVec(2,:),'MinPeakDistance',1/freq*lcl_analyVar.fSamp);
       [~,minPntsLoc] = findpeaks(-tmpIntVec(2,:),'MinPeakDistance',1/freq*lcl_analyVar.fSamp);
       p2p            = mean(tmpIntVec(2,maxPntsLoc)) - mean(tmpIntVec(2,minPntsLoc));
       
       
       % Organize batch data into structure
       intenData.intVec{k} = tmpIntVec(2,:);
       intenData.p2p(k)    = p2p;
    end
    
    % Save intensity data into output structure
    intStruct{basenameNum} = intenData;
end
end