function time_series = NumTim_time_series(save_path, new_subjNames, paths)
%% creates a structure called time_series that contains coordinate ids, and their time series for all datatypes areas for each participant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% time_series_path: where you want the resulting structure to be saved
% new_subjNames: subject names you want in the structure
% paths: where the information about each participant is stored (should be
%       as long as and corresponding with new_subjNames)
%
%
% Output
% time_series: structure with coordinate ids and their time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% general info
TimingMapNames=["TLO", "TTOP", "TTOA", "TPO", "TLS", "TPCI", "TPCM", "TPCS", "TFI", "TFS", "GrayAll"];
TimingMapFiles = [strcat("Left",TimingMapNames(1:end-1),"Map.mat"), strcat("Right",TimingMapNames(1:end-1),"Map.mat"), "gray-Layer1.mat"];
TimingMapNames=[TimingMapNames(1:end-1), TimingMapNames(1:end-1), "GrayAll"];

NumerosityMapNames=["NLO", "NTO", "NPO", "NPCI", "NPCM", "NPCS", "NFI", "NFS", "GrayAll"];
NumerosityMapFiles = [strcat("Left",NumerosityMapNames(1:end-1),"Map.mat"), strcat("Right",NumerosityMapNames(1:end-1),"Map.mat"), "gray-Layer1.mat"];
NumerosityMapNames=[NumerosityMapNames(1:end-1),NumerosityMapNames(1:end-1), "GrayAll"];

mapNames=[TimingMapNames, NumerosityMapNames];
mapFiles = [TimingMapFiles, NumerosityMapFiles];

TimingDTs=["TimingSweeps (L1)", "OddScans (L1)", "EvenScans (L1)"];
NumerosityDTs=["NumerosityAll (L1)", "NumerosityAllOdd (L1)", "NumerosityAllEven (L1)"];
BothDTs=[TimingDTs,NumerosityDTs];


%% get data
for subj = 1:length(new_subjNames)

    cd(paths{subj})
    load('mrSESSION.mat')
    VOLUME{1} = initHiddenGray;

    %get model info voxels
    for dt = 1:length(BothDTs)
        dt_id = find({dataTYPES.name} == BothDTs(dt));
        VOLUME{1}.curDataType = dt_id;
        stimName = regexprep(BothDTs(dt), '[^\w'']','');
        
        % get info for specific roi
        for Map = 1:length(mapNames(1:end-1))

            ROI_path = strcat(paths{subj},'/Gray/ROIs/',mapFiles(Map));
            
            % need to check if file exists, because not all
            % subjects have all maps
            if exist(ROI_path, 'file') == 2
                load(fullfile(ROI_path))                         %Load ROIs
                [tmp, iCrds] = intersectCols(VOLUME{1}.coords, ROI.coords);         %Isolate ROI locations

                if contains(mapFiles(Map),"Right")
                    Hemisphere = "Right";
                    time_series.(new_subjNames{subj}).(mapNames(Map)).(Hemisphere).(stimName).iCoords = iCrds;
                elseif contains(mapFiles(Map),"Left")
                    Hemisphere = "Left";
                    time_series.(new_subjNames{subj}).(mapNames(Map)).(Hemisphere).(stimName).iCoords = iCrds;
                else
                    time_series.(new_subjNames{subj}).(mapNames(Map)).(stimName).iCoords = iCrds;
                end
                    
                scanFolders = dir(strcat(paths{subj},'/Gray/', BothDTs{dt}, '/TSeries/Scan*'));
                
                for scan = 1:length(scanFolders)
                    cd(strcat(paths{subj},'/Gray/', BothDTs(dt), '/TSeries/', scanFolders(scan).name));
                    load('tSeries1');

                    if contains(mapFiles(Map),"Right") ||contains(mapFiles(Map),"Left")
                        time_series.(new_subjNames{subj}).(mapNames(Map)).(Hemisphere).(stimName).(scanFolders(scan).name) = tSeries(:,iCrds);
                    else
                        time_series.(new_subjNames{subj}).(mapNames(Map)).(stimName).(scanFolders(scan).name) = tSeries(:,iCrds);
                    end

                end
                
                clear tSeries iCrds
                
                
            else
                fprintf('Map does not exist: %s for %s \n',mapNames(Map), new_subjNames{subj})
            end
            
            
        end
    end
    mrvCleanWorkspace
    close all
end

cd(save_path);
savename = 'time_series.mat';
save(savename, 'time_series','-v7.3');

