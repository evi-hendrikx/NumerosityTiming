function NumTim_data = NumTim_load_data(NumTim_path,new_subjNames, paths, minVEtime, minVEnum)
%% creates a structure called NumTim_data that contains coordinate ids, variance explained, rss, rawrss,and model parameters for all, even, and odd runs for both models in all areas for each participant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% NumTim_path: where you want the resulting structure to be saved
% new_subjNames: subject names you want in the structure
% paths: where the information about each participant is stored (should be
%       as long as and corresponding with new_subjNames)
% minVEtime: voxels are included as responding to time above this threshold
% minVEnum: voxels are included as responding to numerosity above this threshold
%
%
% Output
% NumTim_data: structure with
% subjects-->map-->hemisphere-->variable-->DTrun-->info to run analyses%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see if it already exists
if (exist(NumTim_path, 'file') == 2)
    load(NumTim_path);
    return
end


%% general info
TimingMapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI", "TFS", "GrayAll"];
NumerosityMapNames=["NLO", "NTO", "NPO", "NPCI", "NPCM", "NPCS", "NFI", "NFS", "GrayAll"];
mapNames=[TimingMapNames, NumerosityMapNames];

modelFolders=["NumTimingModels"];
modelNames{1}=["Timing", "Numerosity"]; %Names for models in output structure (for each folder)
modelFileNames{1}=["Oval", "1DGaussian-Log"]; %Unique parts of corresponding model file names (for each folder)

TimingDTs=["TimingSweeps (L1)", "OddScans (L1)", "EvenScans (L1)"];
NumerosityDTs=["NumerosityAll (L1)", "NumerosityAllOdd (L1)", "NumerosityAllEven (L1)"];
BothDTs=[TimingDTs,NumerosityDTs];

%% load in data with this script that already existed
% could definitely be prettier/ more efficient, but this already existed
% and I didn't want to touch it too much

for subj=1:length(new_subjNames)
    try
        close(1); mrvCleanWorkspace;
    end
    NumTim_data.(new_subjNames{subj}) = getTimingNumerosityModelData(paths{subj}, modelFolders, modelNames, modelFileNames, BothDTs, mapNames);
    % NOTE TO SELF: how indices were computed: [tmp1, roiIndices]=intersectCols(view.coords, view.ROIs(sourceRoi).coords);
end


%% adjust the structure to indicate which voxels are included after thresholding

Hemispheres=["Left", "Right"];
comparisonType = {'Timing','Numerosity'};

for comparison = 1:length(comparisonType)
    for subj=1:length(new_subjNames)
        for Hemisphere=1:length(Hemispheres)
            if string(comparisonType{comparison}) == "Timing"
                start_id = 1; end_id = length(TimingMapNames);
            elseif string(comparisonType{comparison}) == "Numerosity"
                start_id = length(TimingMapNames)+1; end_id = length(TimingMapNames)+length(NumerosityMapNames);
            end

            DTruns = strcat(comparisonType{comparison},{'All' 'Odd' 'Even'});

            for Map=start_id:end_id
                for run = 1:length(DTruns)
                    try
                        x0s = NumTim_data.(new_subjNames{subj}).(mapNames{Map}).(Hemispheres{Hemisphere}).(comparisonType{comparison}).(DTruns{run}).x0s;
                    catch % if the map doesn't exist for this participant
                        sprintf(" NO %s %s for %s for %s", Hemispheres{Hemisphere}, mapNames{Map}, new_subjNames{subj}, DTruns{run})
                        continue
                    end

                    if string(comparisonType{comparison}) == "Timing"
                        thresh_id_time = NumTim_data.(new_subjNames{subj}).(mapNames{Map}).(Hemispheres{Hemisphere}).(comparisonType{comparison}).(DTruns{run}).ves > minVEtime;
                        y0s = NumTim_data.(new_subjNames{subj}).(mapNames{Map}).(Hemispheres{Hemisphere}).(comparisonType{comparison}).(DTruns{run}).y0s;
                        max_x_id_time = x0s<1;
                        min_x_id_time = x0s>0.05;
                        max_y_id_time = y0s<1;
                        min_y_id_time = y0s>0.05;
                        MapThresh = thresh_id_time & max_x_id_time & min_x_id_time & max_y_id_time & min_y_id_time;
                    elseif string(comparisonType{comparison}) == "Numerosity"
                        thresh_id_num = NumTim_data.(new_subjNames{subj}).(mapNames{Map}).(Hemispheres{Hemisphere}).(comparisonType{comparison}).(DTruns{run}).ves > minVEnum;
                        max_x_id_num = x0s<log(7);
                        min_x_id_num = x0s>log(1.01);
                        MapThresh = thresh_id_num & max_x_id_num & min_x_id_num;
                    end

                    VOLUME_ids = NumTim_data.(new_subjNames{subj}).(mapNames{Map}).(Hemispheres{Hemisphere}).(comparisonType{comparison}).(DTruns{run}).roiIndices;
                    NumTim_data.(new_subjNames{subj}).(mapNames{Map}).(Hemispheres{Hemisphere}).(comparisonType{comparison}).(DTruns{run}).IndicesAboveThreshold = find(MapThresh==1);
                    NumTim_data.(new_subjNames{subj}).(mapNames{Map}).(Hemispheres{Hemisphere}).(comparisonType{comparison}).(DTruns{run}).VOLUMEIndicesAboveThreshold = VOLUME_ids(MapThresh==1);
                end
            end
        end
    end
end

save(NumTim_path, 'NumTim_data');
end



