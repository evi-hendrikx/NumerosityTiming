function data = getTimingNumerosityModelData(path, modelFolders, modelNames, modelFileNames,  DTnamesIn, mapNames)
%% Get data from various candidate timing models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% path: where the information about the participant is stored 
% modelFolders: folder in which the retinotopy models are saved gray -->
%       DT_run --> FOLDER
% modelNames: names of quantities you want to compare
% modelFileNames: unique identifier in retinotopy models corresponding to
%       modelNames
% DTnamesIn: DT_runs names (e.g., OddScans (L1) )
% mapNames: names of the maps you want to include
%
%
% Output
% data: structure with
% map-->hemisphere-->variable-->DTrun-->info to run analyses%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(path)
mrGlobals;
mrVista 3;
path=char(path);

% get indices in mrVista from requested datatypes
for m=1:length(DTnamesIn)
    for n=1:length(dataTYPES)
        if dataTYPES(n).name==DTnamesIn(m)
            DTs(m)=n;
        end
    end
end

% what you want the datatypes to be renamed as in your new structure
DTnames=["TimingAll", "TimingOdd", "TimingEven", "NumerosityAll", "NumerosityOdd", "NumerosityEven"];

%% Load all ROIs
VOLUME{1}=deleteAllROIs(VOLUME{1}); VOLUME{1}=refreshScreen(VOLUME{1},0);
okROIs=zeros(length(mapNames), 2);
for whichMap=1:length(mapNames)

    [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Left', mapNames(whichMap), 'Map')), [],[],[],1);
    if ok
        okROIs(whichMap,1)=1;
    else
         % still load an ROI so the indexing doesn't go crazy
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 1]);
    end

    [VOLUME{1}, ok] = loadROI(VOLUME{1}, char(strcat('Right', mapNames(whichMap), 'Map')), [],[],[],1);
    if ok
        okROIs(whichMap,2)=1;
    else
         % still load an ROI so the indexing doesn't go crazy
        VOLUME{1} = makePointROI(VOLUME{1},[1 whichMap 2]);
    end
end
VOLUME{1} = refreshScreen(VOLUME{1}, 0);


%% get info
for whichDT=1:length(DTs)
    %     if DTs(whichDT)>0
    VOLUME{1}=viewSet(VOLUME{1}, 'curdt', DTs(whichDT));
    for whichFolder=1:length(modelFolders)

        %% select and load model 
        folder=[path filesep 'Gray' filesep dataTYPES(DTs(whichDT)).name filesep char(modelFolders(whichFolder))];
        whichModelNames=modelNames{whichFolder};
        whichModelFiles=modelFileNames{whichFolder};

        if whichDT<4 % timing 
            whichModel=1;
        else % numerosity
            whichModel=2;
        end
        modelFile=dir([char(strcat(folder,filesep, '*', whichModelFiles(whichModel), '*'))]);
        disp(modelFile)

        VOLUME{1}=rmSelect(VOLUME{1}, 1,[folder, filesep, modelFile.name]);

        try
            VOLUME{1} = rmLoadDefault(VOLUME{1});
        catch
            VOLUME{1}
        end
        VOLUME{1}=refreshScreen(VOLUME{1},0);


        %% get data for specific model
        for whichMap=1:length(mapNames)
            dataName=char(strcat('data.', mapNames(whichMap), '.Left.', whichModelNames(whichModel), '.', DTnames(whichDT)));

            %Left hemisphere first entries
            if okROIs(whichMap,1)
                if whichDT==1 && whichFolder==1 && whichModel==1
                    dataTmp = RoiDistanceRatio(VOLUME{1}, [], [], 1+(whichMap-1)*2);

                % I think this is similar to the previous script that's
                % called, but the indices are already available so that
                % doesn't need to happen anymore
                else
                    sourceDataName=char(strcat('data.', mapNames(whichMap), '.Left.', modelNames{1}(1), '.', DTnames(1)));
                    try
                        eval(['dataTmp=', sourceDataName, ';'])
                    catch
                        sourceDataName
                        %
                    end
                    dataTmp.x0s=VOLUME{1}.rm.retinotopyModels{1}.x0(dataTmp.roiIndices);
                    dataTmp.y0s=VOLUME{1}.rm.retinotopyModels{1}.y0(dataTmp.roiIndices);
                    dataTmp.sigmas=VOLUME{1}.rm.retinotopyModels{1}.sigma.major(dataTmp.roiIndices);
                    dataTmp.sigmaMinor=VOLUME{1}.rm.retinotopyModels{1}.sigma.minor(dataTmp.roiIndices);
                    dataTmp.sigmaTheta=VOLUME{1}.rm.retinotopyModels{1}.sigma.theta(dataTmp.roiIndices);
                    if isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exp')
                        dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exp(dataTmp.roiIndices);
                    elseif isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exponent')
                        dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exponent(dataTmp.roiIndices);
                    end

                    ves=1-(VOLUME{1}.rm.retinotopyModels{1}.rss(dataTmp.roiIndices)./VOLUME{1}.rm.retinotopyModels{1}.rawrss(dataTmp.roiIndices));
                    ves(~isfinite(ves)) = 0;
                    ves = max(ves, 0);
                    ves = min(ves, 1);
                    dataTmp.ves=ves;

                end
                eval([dataName, '=dataTmp;'])

            end

            dataName=char(strcat('data.', mapNames(whichMap), '.Right.', whichModelNames(whichModel), '.', DTnames(whichDT)));
            %Right hemisphere

            if okROIs(whichMap,2)
                if whichDT==1 && whichFolder==1 && whichModel==1
                    dataTmp = RoiDistanceRatio(VOLUME{1}, [], [], 2+(whichMap-1)*2);
                else
                    sourceDataName=char(strcat('data.', mapNames(whichMap), '.Right.', modelNames{1}(1), '.', DTnames(1)));
                    try
                        eval(['dataTmp=', sourceDataName, ';'])
                    catch
                        sourceDataName
                    end
                    dataTmp.x0s=VOLUME{1}.rm.retinotopyModels{1}.x0(dataTmp.roiIndices);
                    dataTmp.y0s=VOLUME{1}.rm.retinotopyModels{1}.y0(dataTmp.roiIndices);
                    dataTmp.sigmas=VOLUME{1}.rm.retinotopyModels{1}.sigma.major(dataTmp.roiIndices);
                    dataTmp.sigmaMinor=VOLUME{1}.rm.retinotopyModels{1}.sigma.minor(dataTmp.roiIndices);
                    dataTmp.sigmaTheta=VOLUME{1}.rm.retinotopyModels{1}.sigma.theta(dataTmp.roiIndices);
                    if isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exp')
                        dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exp(dataTmp.roiIndices);
                    elseif isfield(VOLUME{1}.rm.retinotopyModels{1}, 'exponent')
                        dataTmp.exp=VOLUME{1}.rm.retinotopyModels{1}.exponent(dataTmp.roiIndices);
                    end

                    ves=1-(VOLUME{1}.rm.retinotopyModels{1}.rss(dataTmp.roiIndices)./VOLUME{1}.rm.retinotopyModels{1}.rawrss(dataTmp.roiIndices));
                    ves(~isfinite(ves)) = 0;
                    ves = max(ves, 0);
                    ves = min(ves, 1);
                    dataTmp.ves=ves;
                end
                eval([dataName, '=dataTmp;'])
            end

        end
    end
end


close(1); mrvCleanWorkspace;
end

