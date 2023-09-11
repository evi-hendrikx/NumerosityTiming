function topo_vector = NumTim_get_topo(path, NumTim_data,subjName, HemiName, MapName, DTrun)
%% Calculated topographic vectors for all voxels in the requested map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% path: where the information about the participant is stored 
% NumTim_data: made with the NumTim_load_data script. Contains all values
%       that you could need for all following analyses
% subjName: subject for which you want to get the topographic vectors
% HemiName: hemisphere for which you want to get the topographic vectors
% MapName: map for which you want to get the topographic vectors
% DTrun: the first data run you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
%
%
% Output
% topo_vector: struct with normals, cordinates ids and parameter vectors of the
%       requested map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% general info
cd(path)

%% get NumTim coords here
VOLUME{1} = initHiddenGray;
all_coords = VOLUME{1}.coords;
load('mrSESSION.mat');
mmPerPix = readVolAnatHeader(vANATOMYPATH);

%% Only layer 1 nodes
% %     % Only needs to be done once for each subject
% %     segFile = dir([char('*class*')]);
% %     segPaths{1}=segFile.name;
% %     segPaths{2}=segFile.name;
% %     buildGrayCoords(VOLUME{1}, 'L1coords', 0, segPaths, 0);

load('L1coords.mat');
if string(HemiName) == "Left"
    nodes = double(L1coords.allLeftNodes);
    edges = double(L1coords.allLeftEdges);
else
    nodes = double(L1coords.allRightNodes);
    edges = double(L1coords.allRightEdges);
end

%% load mesh to find normals
% have to work with hemispheres separately, otherwise mrVista's distance formula doesn't
% work well
if (string(HemiName) == "Left" && ~isfield(L1coords, 'vertexGrayMapLeft')) || (string(HemiName) == "Right" && ~isfield(L1coords, 'vertexGrayMapRight'))
    meshPath = strcat([pwd,'/folded_layer0.mat']);
    [VOLUME{1},OK] = meshLoad(VOLUME{1},meshPath,1);
    MSH = viewGet(VOLUME{1}, 'Mesh');
    vertexGrayMap = mrmMapVerticesToGray(meshGet(MSH, 'initialvertices'), nodes, mmPerPix, edges);
    
    allMeshes = viewGet(VOLUME{1}, 'allmeshes');
    mrmSet(allMeshes, 'closeall');
    h = findobj('Name', '3DWindow (mrMesh)');
    close(h);
    
    if string(HemiName) == "Left"
        L1coords.vertexGrayMapLeft = vertexGrayMap;
    else
        L1coords.vertexGrayMapRight = vertexGrayMap;
    end
    L1coords.all_normals = VOLUME{1}.mesh{1}.normals;
    
    save('L1coords.mat','L1coords')
else
    if string(HemiName) == "Left"
        vertexGrayMap = L1coords.vertexGrayMapLeft;
    else
        vertexGrayMap = L1coords.vertexGrayMapRight;
    end
end


%% load ROI into VOLUME so it can be used for computations
VOLUME{1} = loadROI(VOLUME{1}, strcat(HemiName,MapName,'Map'), [],[],[],1);
map_bool = arrayfun(@(x) string(x.name)==string(strcat(HemiName,MapName,'Map')),VOLUME{1}.ROIs);
map_id = find(map_bool == 1);
if isempty(map_id)
    return
end
targetROI = VOLUME{1}.ROIs(map_id);
sourceROI = VOLUME{1}.ROIs(map_id);

%% make sure all info is available for included nodes
%  'nodes' is an 8xn array. The first 3 rows correspond to the y, x, and z coordinates of
%   the nodes %% remappedNodes has same order as coords
remappedNodes = [nodes(2,:); nodes(1,:); nodes(3,:)];

% find ROI's coords in remapped nodes
[~, sourceNodeIndices_remap, ~]=intersectCols(remappedNodes, VOLUME{1}.ROIs(map_id).coords); %% importantly this is the same order we use as the creation of NumTim_data

% some ROIs exist from a sligtly other segmentations
[~, sourceNodeIndices_remap_id, ~] = intersectCols(remappedNodes(:,sourceNodeIndices_remap), L1coords.coords);

% sourceNodeIndices indices of L1 ROI in remappedNodes
sourceNodeIndices = sourceNodeIndices_remap(sourceNodeIndices_remap_id);
lineNodeIndices = sourceNodeIndices;
targetROI.coords = remappedNodes(:,sourceNodeIndices);
sourceROI.coords = remappedNodes(:,sourceNodeIndices);

%% find normals
% find ids of the coords in the mesh
% if there is multiple, take the mean
normals = NaN(3,length(sourceNodeIndices));
for sourceNodeID = 1:length(sourceNodeIndices)
    normal_ids = find(vertexGrayMap == sourceNodeIndices(sourceNodeID));
    if isempty(normal_ids)
            error('Not all VOLUME coordinates can be found in mesh')
    end
    normals(:,sourceNodeID) = mean(L1coords.all_normals(:,normal_ids),2);
end

topo_vector.norm_surf = normals;

%% get distances to other voxels
allDist = zeros(numel(sourceNodeIndices), numel(lineNodeIndices));
last_points = zeros(numel(sourceNodeIndices), length(remappedNodes));
for g=1:numel(sourceNodeIndices)
    [tmp,~,last_point] = mrManDist(nodes, edges, double(sourceNodeIndices(g)), mmPerPix, 100000, 0);
    allDist(g,:) = tmp(lineNodeIndices);
    last_points(g,:) = last_point;
end

%% Make topo's timing and numerosity separately
% look at indices above the threshold for the specific variable, later look at only the shared voxels
% I already removed out of bounds/ VE threshold before to make IndicesAboveThreshold when loading in the data

if startsWith(MapName,"T") && ~(contains(MapName,"All"))
    parameter_maps = {'durationMap','periodMap'};

    % get parameter values for the voxels in the map
    above_thresh = NumTim_data.(subjName).(MapName).(HemiName).Timing.(DTrun).IndicesAboveThreshold;
    durationMap = NaN(size(NumTim_data.(subjName).(MapName).(HemiName).Timing.(DTrun).x0s));
    periodMap = durationMap;
    durationMap(above_thresh) = NumTim_data.(subjName).(MapName).(HemiName).Timing.(DTrun).x0s(above_thresh);
    periodMap(above_thresh) = NumTim_data.(subjName).(MapName).(HemiName).Timing.(DTrun).y0s(above_thresh);
    
    % voxels witin map (or within roiIndices of the entire brain volume)
    % are structured differently from voxels within only one of the
    % hemispheres. We need the value data and the normal (etc.) data to be
    % structured so it's clear what belongs to the same coordinates
    roi_coord_id = NumTim_data.(subjName).(MapName).(HemiName).Timing.(DTrun).roiIndices; % order of data in NumTim
    coords_time = all_coords(:,roi_coord_id); % get coords by selecting from all voxels in VOLUME of the entire brain (this was how roiIndices was originally made)
    % match the order of the values of NumTim to the order of the remappedNodes
    [~,~,Tim_ids] = intersectCols(remappedNodes(:,sourceNodeIndices), coords_time); 
    durationMap = durationMap(Tim_ids);
    periodMap = periodMap(Tim_ids);
    topo_vector.included_coord_ids = roi_coord_id(Tim_ids);
 
    % sanity check
    if length(Tim_ids) ~= length(coords_time) || length(Tim_ids) ~= length(sourceNodeIndices)
        error('Map in NumTim dataset and newly loaded in map do not have the same amount of voxels')
    end
    
elseif contains(MapName,"All")
    parameter_maps = {};

elseif startsWith(MapName,"N") && ~(contains(MapName,"All"))
    parameter_maps = {'numerosityMap'};

    above_thresh = NumTim_data.(subjName).(MapName).(HemiName).Numerosity.(DTrun).IndicesAboveThreshold;
    numerosityMap = NaN(size(NumTim_data.(subjName).(MapName).(HemiName).Numerosity.(DTrun).x0s));
    numerosityMap(above_thresh) = NumTim_data.(subjName).(MapName).(HemiName).Numerosity.(DTrun).x0s(above_thresh);
    
    roi_coord_id = NumTim_data.(subjName).(MapName).(HemiName).Numerosity.(DTrun).roiIndices;
    coords_num = all_coords(:,roi_coord_id);
    [~,~,Num_ids] = intersectCols(remappedNodes(:,sourceNodeIndices), coords_num);
    numerosityMap = numerosityMap(Num_ids);
    topo_vector.included_coord_ids = roi_coord_id(Num_ids);
    
    if length(Num_ids) ~= length(coords_num) || length(Num_ids) ~= length(sourceNodeIndices)
        error('Map in NumTim dataset and newly loaded in map do not have the same amount of voxels')
    end
end


%% calculate topographic vectors
for par_map = 1:size(parameter_maps,2)
    map = eval(parameter_maps{par_map});
    
    % it would be faster if made a matrix and average that outside of the loop
    mean_vector = nan([3,length(sourceNodeIndices)]);
    cross_mean_vector = mean_vector;
    cross_cross_mean_vector = mean_vector;
    if string(parameter_maps{par_map}) == "numerosityMap"
        mean_vector_exp = mean_vector;
        cross_mean_vector_exp = mean_vector;
        cross_cross_mean_vector_exp = mean_vector;
    end
    
    for gg = 1:length(sourceNodeIndices)
        
        % surface normal gg
        norm_surf = normals(:,gg);

        vectors_from_g = nan([3,numel(lineNodeIndices)]);
        if string(parameter_maps{par_map}) == "numerosityMap"
            vectors_from_g_exp = nan([3,numel(lineNodeIndices)]);
        end
        
        startPoint = sourceNodeIndices(gg);
        lastPoint = last_points(gg,:);
        
        for n=1:numel(lineNodeIndices)
            
            if lineNodeIndices(n) ~= sourceNodeIndices(gg)
                %% make geodesic to capture all steps (!) from n to gg (!) and includes gg and n
                endPoint = lineNodeIndices(n);
                nextPoint = endPoint;
                geodesic=endPoint;
                while(nextPoint~=startPoint)
                    geodesic(end+1) = lastPoint(nextPoint); %#ok<SAGROW>
                    nextPoint = lastPoint(nextPoint);
                end
                
                angle_cross=NaN(length(geodesic)-1,1);
                distances=zeros(length(geodesic)-1,1);
                for step = length(geodesic)-1 :-1:1
                    % here x and y are switched again, so use the unaltered normals?
                    ai=remappedNodes(2,geodesic(step))-...
                        remappedNodes(2,geodesic(end));
                    bj=remappedNodes(1,geodesic(step))-...
                        remappedNodes(1,geodesic(end));
                    ck=remappedNodes(3,geodesic(step))-...
                        remappedNodes(3,geodesic(end));
                    distances(step)=sqrt(ai^2+bj^2+ck^2);
                    step_vec(:,step) = [ai/ distances(step), bj/ distances(step), ck/ distances(step)];
                    
                    % since both use unit vectors, length of
                    % AxB is same, just the angle is different
                    % for different steps
                    if step == length(geodesic)-1
                        first_cross = cross(step_vec(:,step), norm_surf);
                    end
                    crossproduct=cross(step_vec(:,step), norm_surf);
                    
                    % this angle is between 0 and 180, which is fine for what we're doing here (checking it does not exceed 90)
                    cos_angle_timing_numerosity = max(min(dot(first_cross,crossproduct)./(norm(first_cross).*norm(crossproduct)),1),-1);
                    angle_cross(step) = real(acosd(cos_angle_timing_numerosity));
                end
                
                %% max distance where cross product where the vector is <90 deg first step cross-product
                under_90 = find(angle_cross < 90);
                [~,I] = max(distances(under_90));
                unit_vec_surface = step_vec(:,under_90(I));
                % normalize change in preferred value with cortical
                % distance
                quantity_gradient = (map(n) - map(gg))/allDist(gg,n);
                % here vector lengths will change, which is
                % important for the calculation of the mean
                vectors_from_g(:,n) = quantity_gradient .* unit_vec_surface;
                
                if string(parameter_maps{par_map}) == "numerosityMap"
                    exp_quantity_gradient = ((exp(map(n)) - exp(map(gg)))/allDist(gg,n));
                    vectors_from_g_exp(:,n) = exp_quantity_gradient .* unit_vec_surface;
                end
                
            end
        end
        
        mean_vector(:,gg) = nanmean(vectors_from_g,2);
        %% for later calculations and visualizations
        cross_mean_vector(:,gg)= cross(mean_vector(:,gg),normals(:,gg));
        cross_cross_mean_vector(:,gg)= cross(normals(:,gg),cross_mean_vector(:,gg));
        
        if string(parameter_maps{par_map}) == "numerosityMap"
            mean_vector_exp(:,gg) = nanmean(vectors_from_g_exp,2);
            
            %% for later calculations and visualizations
            cross_mean_vector_exp(:,gg)= cross(mean_vector_exp(:,gg),normals(:,gg));
            cross_cross_mean_vector_exp(:,gg)= cross(normals(:,gg),cross_mean_vector_exp(:,gg));
        end
    end
    
    %% save vector field info
    topo_vector.(parameter_maps{par_map}).mean_vector = mean_vector;
    topo_vector.(parameter_maps{par_map}).cross_mean_vector = cross_mean_vector;
    topo_vector.(parameter_maps{par_map}).cross_cross_mean_vector = cross_cross_mean_vector;
    if string(parameter_maps{par_map}) == "numerosityMap"
        topo_vector.(parameter_maps{par_map}).mean_vector_exp = mean_vector_exp;
        topo_vector.(parameter_maps{par_map}).cross_mean_vector_exp = cross_mean_vector_exp;
        topo_vector.(parameter_maps{par_map}).cross_cross_mean_vector_exp = cross_cross_mean_vector_exp;
    end
end

mrvCleanWorkspace;
close all;

end