%% Makes the surface renderings for the timing response model data (the numerosity model fits were done manually)

% Paths to find subject data
paths{1}=[anonymized]
paths{2}=[anonymized]
paths{3}=[anonymized]
paths{4}=[anonymized]
paths{5}=[anonymized]
paths{6}=[anonymized]
paths{7}=[anonymized]
paths{8}=[anonymized]


meshNames = {'meshLeftInflated','meshRightInflated'};

for subj = 1:length(paths)
    cd(paths{subj})
    
    for time_id = 1:length({'duration','period'})
        if time_id == 1
            parameter = 'map';
        else
            parameter = 'amp';
        end
        
        
        for meshes=1:length(meshNames) % back left right
            
            % set data type to timing (all the data for timing, not just
            % odd or even)
            mrVista 3
            dt_id = find({dataTYPES.name} == "TimingSweeps (L1)");
            VOLUME{1}.curDataType = dt_id;
            
            %Get correct folder for stimulus configuration
            stimFolder = fullfile(paths{subj},'Gray',dataTYPES(1,dt_id).name,'NumTimingModels');
            rmName = char(strcat([stimFolder,'/', ls(stimFolder)]));
            
            %load model
            VOLUME{1} = rmSelect(VOLUME{1}, 1, rmName); VOLUME{1} = rmLoadDefault(VOLUME{1});
            
            % set duration and period as "map"  and "amp"  within mrVista (DON'T SET
            % AS PHASE), so they can be visualized
            VOLUME{1} = rmLoad(VOLUME{1}, 1, 'x0', 'map');
            VOLUME{1} = rmLoad(VOLUME{1}, 1, 'y0', 'amp');
            
            % set color bar
            VOLUME{1} = setClipMode(VOLUME{1}, 'map', [0 1.1]);
            VOLUME{1} = setClipMode(VOLUME{1}, 'amp', [0 1.1]);
            VOLUME{1} = refreshScreen(VOLUME{1}, 1);
            VOLUME{1}.ui.mapMode=setColormap(VOLUME{1}.ui.mapMode, 'hsvCmap');
            VOLUME{1}.ui.ampMode=setColormap(VOLUME{1}.ui.ampMode, 'hsvCmap');
            
            % only parameters inside de range
            VOLUME{1}.co{1}(VOLUME{1}.amp{1}>1 | VOLUME{1}.amp{1}<0.05)=0;
            VOLUME{1}.co{1}(VOLUME{1}.map{1}>1 | VOLUME{1}.map{1}<0.05)=0;
            
            VOLUME{1}=refreshScreen(VOLUME{1}, 1);
            
            % set coherence threshold (0.2 fig 2, 0.1 supplementary)
            VOLUME{1} = setCothresh(VOLUME{1}, 0.2);VOLUME{1} = refreshScreen(VOLUME{1}, 1);
            
            VOLUME{1} = setDisplayMode(VOLUME{1},parameter); VOLUME{1} = refreshScreen(VOLUME{1});
            
            
            % load 3d mesh brain
            H = open3DWindow;
            
            % use previously saved settings for lighting etc.
            meshSettings = load('MeshSettings.mat');
            subjSettings = length(meshSettings.settings);
            currentSettings = [subjSettings-1,subjSettings]; %left,right
            
            % load mesh
            meshPath = strcat([pwd,'/',meshNames{meshes}]);
            [VOLUME{1},OK] = meshLoad(VOLUME{1},meshPath,1);
            msh = viewGet(VOLUME{1}, 'currentmesh'); % refresh happens there
            if ~isempty(msh)
                VOLUME{1} = viewSet(VOLUME{1},'currentmesh',msh);
            end
            mrmSet(msh, 'cursoroff');
            
            MSH = viewGet(VOLUME{1}, 'Mesh');
            vertexGrayMap = mrmMapVerticesToGray( meshGet(MSH, 'initialvertices'), viewGet(VOLUME{1}, 'nodes'),...
                viewGet(VOLUME{1}, 'mmPerVox'), viewGet(VOLUME{1}, 'edges') );
            MSH = meshSet(MSH, 'vertexgraymap', vertexGrayMap); VOLUME{1} = viewSet(VOLUME{1}, 'Mesh', MSH);
            clear MSH vertexGrayMap
            meshRetrieveSettings(viewGet(VOLUME{1}, 'CurMesh'), currentSettings(meshes));
            
            VOLUME{1} = viewSet(VOLUME{1}, 'mesh', msh);
            VOLUME{1} = meshColorOverlay(VOLUME{1});
            
            % screenshot currently done manually to match sizes ROI and
            % numerosity files
% %             wSize=[1056,1855]; %@1920x1080
% %             mrmSet(msh,'windowSize',wSize(1),wSize(2));
% %             
% %             
% %             fname = fullfile(save_mesh,[saveNames{meshes}, '_', subjNames{subj},'.png']);
% %             udata.rgb = mrmGet(msh,'screenshot')/255;
% %             imwrite(udata.rgb, fname);
% %             set(gcf,'userdata',udata);
% %             
            % nicely close everything
            allMeshes = viewGet(VOLUME{1}, 'allmeshes');
            mrmSet(allMeshes, 'closeall');
            h = findobj('Name', '3DWindow (mrMesh)');
            close(h);
            VOLUME{1}=rmfield(VOLUME{1},'mesh');
            
            close all
            mrvCleanWorkspace
            
        end
    end
end

