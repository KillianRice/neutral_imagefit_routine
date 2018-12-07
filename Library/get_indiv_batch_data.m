function indivDataset = get_indiv_batch_data(analyVar)
% Reads the master file list provided from AnalysisVariables (in analyVar) and
% opens each dataset batch file to create the dataset variables and the
% BackgroundAll matrix
%
% INPUTs:
%   analyVar - Structure from AnalysisVariables which enumerates all the
%              variables needed for the analysis
%
% OUTPUTS:
%   indivDataset - a cell of structures containing the individual dataset variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through each batch file listed in analyVar.basenamevectorAtom
indivDataset = cell(analyVar.numBasenamesAtom,1);
for basenameNum = 1:analyVar.numBasenamesAtom
    % Find basename for background
    whichBack = strcmpi(analyVar.basenamevectorBack,char(analyVar.basenamevectorAtom(basenameNum)));
    f = analyVar.basenamevectorBack{whichBack};
    %%%need to catch error about reference non-cell array when no names found in background file
    
    batchfileAtom = [analyVar.dataDir char(analyVar.basenamevectorAtom(basenameNum)) '.batch']; % current atom batchfile name
    batchfileBack = [analyVar.dataDir f '.batch'];                                           % current background batchfile name
    
    %read in all atom files, no limit
    try
        indivBatchAtomData = textscan(fopen(batchfileAtom), getVarList(analyVar),'commentstyle','%');
    catch atomBatchFileErr
        error(['Invalid file identifier. Unable to locate \n%s\n',...
               'Check that the background master batch (Files_%s.txt) is properly configured.'],...
               batchfileAtom,strrep(analyVar.dataDirName,'.',''))
    end
    %read in all background files, no limit
    
    try
        indivBatchBackData = textscan(fopen(batchfileBack), getVarList(analyVar),'commentstyle','%');
    catch backBatchFileErr
        error(['Invalid file identifier. Unable to locate \n%s\n',...
               'Check that the background master batch (Files_%s_Bg.txt) is properly configured.'],...
               batchfileBack,strrep(analyVar.dataDirName,'.',''))
    end
    % create structure with the variables from the filename
    indivBatch = cell2struct(cat(2,indivBatchAtomData,indivBatchBackData),cat(2,analyVar.indivBatchAtomVar,analyVar.indivBatchBackVar),2);
    
    indivBatch.CounterAtom = size(indivBatch.fileAtom,1); % determine how many atom files in batch associated with this basename
    indivBatch.CounterBack = size(indivBatch.fileBack,1); % determine how many background files in batch associated with this basename
    
    % Subplot number for plotting program
    [indivBatch.SubPlotRows, indivBatch.SubPlotCols] = optiSubPlotNum(indivBatch.CounterAtom);
    
    % Retrieve indexing matricies that outline the different regions of the image
    [indivBatch.image_Index,indivBatch.roiWin_Index] = get_image_regions(analyVar,basenameNum); 
    
    % Find number of elements inside the cloud window and around the cloud within the ROI window
    elCloud    = numel(nonzeros(indivBatch.image_Index ~= 0)); %% Enumerates all elements within the roiWinRadAtom (cloud and around)
    elNotCloud = numel(nonzeros(indivBatch.image_Index == 1)); %% Enumerates only elements not within the cloud
    
    %% Open raw images from camera and parse
    % The raw data is only needed before saving OD image (only OD image needed for fits 
    % and plotting) so this can be skipped for time. 
    noRawFile = {'imagefit_NumDistFit', 'imagefit_ParamEval'};
    %curStack = struct2cell(dbstack); 
    %callParent = curStack(2,length(curStack(2,:)) - 1);
    
    % Trying a new method to identify if the primary caller is one in noRawFile
    curStack = (dbstack); callParent = curStack(end).name; 
    % This will skip reading the raw file if the caller is one of the names in noRawFile
    if ~sum(strcmpi(noRawFile,callParent))
        % Aggregate background files into 2 matricies, 1 of area around cloud and the other including the cloud
        % Initialize matricies for speed
        AtomsCloud     = zeros(elCloud,indivBatch.CounterAtom);
        AtomsNotCloud  = zeros(elNotCloud,indivBatch.CounterAtom);
        rawAtomImages  = zeros(prod(analyVar.matrixSize), indivBatch.CounterAtom);
        for k = 1:indivBatch.CounterAtom
            s = [analyVar.dataDir char(indivBatch.fileAtom(k)) analyVar.dataAtom]; sFID = fopen(s,'r','ieee-be');
            %Read in binary file created in LabView
            %   LabView saves binary in "big endian" ('be') format: most significant bit in
            %   lowest memory address. Matlab needs this info to import binaryfile correctly.
            %   LabView and Matlab treat matrix coordinates differently.
    
            try % Bad image will cause this to fail so added an error handling case
                
                %   (0,0) for LabView is lower left corner; for Matlab is upper left corner.
                rawImageAtom_2D = double(fread(sFID,analyVar.matrixSize,'*int16')); fclose(sFID);
                               
                % Separate Atoms into cloud part and around cloud part
                AtomsCloud(:,k)    = rawImageAtom_2D(indivBatch.image_Index ~= 0);
                AtomsNotCloud(:,k) = rawImageAtom_2D(indivBatch.image_Index == 1);
                rawAtomImages(:,k) = rawImageAtom_2D(:);
                
            catch
                error('There was an error loading the image from batch number %g, image %g. Please disable this image and rerun',basenameNum,k)
            end
        end %%%% end of k = 1:CounterAtom

        % Aggregate background files into 2 matricies, 1 of background behind cloud and the other around the cloud
        % Initialize matricies for speed
        BackgroundCloud    = zeros(elCloud,indivBatch.CounterBack);
        BackgroundNotCloud = zeros(elNotCloud,indivBatch.CounterBack);
        rawBackImages      = zeros(prod(analyVar.matrixSize), indivBatch.CounterAtom);
        for j = 1:indivBatch.CounterBack
            t = [analyVar.dataDir char(indivBatch.fileBack(j)) analyVar.dataBack]; tFID = fopen(t,'r','ieee-be');
            rawImageBack_2D = double(fread(tFID,analyVar.matrixSize,'*int16')); fclose(tFID);

            % Separate Background into part behind cloud and part around cloud
            BackgroundCloud(:,j)    = rawImageBack_2D(indivBatch.image_Index ~= 0);
            BackgroundNotCloud(:,j) = rawImageBack_2D(indivBatch.image_Index == 1);
            rawBackImages(:,j)      = rawImageBack_2D(:);
        end

        % Save Atoms and Background into the indivBatch structure
        indivBatch.AtomsCloud         = AtomsCloud;
        indivBatch.AtomsNotCloud      = AtomsNotCloud;
        indivBatch.rawAtomImages      = rawAtomImages;
        indivBatch.BackgroundCloud    = BackgroundCloud;
        indivBatch.BackgroundNotCloud = BackgroundNotCloud;
        indivBatch.rawBackImages      = rawBackImages;
    end
   
    % In the rare case that you need to apply a calibration function to the
    % x-axis before proceeding on with fitting, the function defined in
    % AnalysisVariables for indFuncCalib will be called here. See the
    % discussion in the "Control the independent variable prensentation"
    % section of AnalysisVariables for further information
    if analyVar.applyIndCalib
        indivBatch = analyVar.funcIndCalib(analyVar, indivBatch);
    end
    
    %%%% Save the variables for each dataset into a variable containing data from all datasets listed in analyVar.basenamevectorAtom
    indivDataset{basenameNum} = orderfields(indivBatch);
end

%% Clean Workspace
fclose all;
end

function lcl_varFormatStr = getVarList(analyVar)
% function to return the string needed for reading in the indivudal image files
    % Read in list of variable names and determine the format string for textscan
    lcl_varFormatStr = repmat({'%f'},1,length(analyVar.indivBatchAtomVar)); % String is as long as number of variables
    lcl_varFormatStr(1) = {'%q'};                                       % Define fixed variables
    lcl_varFormatStr = horzcat(lcl_varFormatStr{:});                    % Concatenate cells into single string for textscan
end
