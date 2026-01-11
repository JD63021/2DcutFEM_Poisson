function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh( ...
    gmshVelFile, order, excludedVelFlags)
% mesh5_gmsh reads a Gmsh file for a triangular approach and builds:
%   nodeInfo.velocity (with node coordinates),
%   elemInfo.velElements (element connectivity),
%   boundaryInfo with boundary node sets.
%
% The mesh file read depends on the input "order":
%   order = 1  => P1 elements: expects 'TRIANGLES' (3 nodes/element) and 'LINES' (2 nodes/line)
%   order = 2  => P2 elements: expects 'TRIANGLES6' (6 nodes/element) and 'LINES3' (3 nodes/line)
%   order = 3  => P3 elements: expects 'TRIANGLES10' (10 nodes/element) and 'LINES4' (4 nodes/line)
%
% If 'excludedVelFlags' is specified, those boundary segments are ignored.

if ~exist('excludedVelFlags','var') || isempty(excludedVelFlags)
    excludedVelFlags = [];
end

%% (A) Read Velocity Mesh
if ischar(gmshVelFile)
    run(gmshVelFile);  % loads struct 'msh' into workspace
    msh_vel = msh;
else
    error('Velocity mesh file not provided as a string.');
end

switch order
    case 1
        if ~isfield(msh_vel, 'TRIANGLES')
            error('Expected TRIANGLES in the velocity mesh (P1).');
        end
        triElementsVel = msh_vel.TRIANGLES(:,1:3);
    case 2
        if ~isfield(msh_vel, 'TRIANGLES6')
            error('Expected TRIANGLES6 in the velocity mesh (P2).');
        end
        triElementsVel = msh_vel.TRIANGLES6(:,1:6);
    case 3
        if ~isfield(msh_vel, 'TRIANGLES10')
            error('Expected TRIANGLES10 in the velocity mesh (P3).');
        end
        triElementsVel = msh_vel.TRIANGLES10(:,1:10);
    otherwise
        error('Order must be 1, 2, or 3.');
end
numElemsVel = size(triElementsVel,1);

allVelCoords = msh_vel.POS;    % Nx3 or Nx2
numVelNodes  = size(allVelCoords,1);

vel_nodes.id = (1:numVelNodes).';
vel_nodes.x  = allVelCoords(:,1);
vel_nodes.y  = allVelCoords(:,2);

%% (B) (Optional Pressure Mesh reading can be added similarly)

%% (C) Store element connectivity
elemInfo.velElements = triElementsVel;  % node count depends on order

%% (D) Fill nodeInfo
nodeInfo.velocity = vel_nodes;

%% (E) Build boundaryInfo from velocity boundary lines.
boundaryInfo = struct();

switch order
    case 1
        % For P1 elements, expect LINES (Nx3: [nodeA, nodeB, flag])
        if ~isfield(msh_vel, 'LINES') || isempty(msh_vel.LINES)
            warning('No LINES found in velocity mesh => boundaryInfo.velLineElements = empty.');
            boundaryInfo.velLineElements = struct();
        else
            lines_vel = msh_vel.LINES;  % Nx3: [nA, nB, flag]
            uniqueVFlags = unique(lines_vel(:,3));
            velLineStruct = struct();
            for iF = 1:numel(uniqueVFlags)
                gFlag = uniqueVFlags(iF);
                if ismember(gFlag, excludedVelFlags)
                    continue;
                end
                mask = (lines_vel(:,3) == gFlag);
                theseLines = lines_vel(mask,1:2);  % columns 1-2: node IDs
                fieldName  = sprintf('flag_%d', gFlag);
                velLineStruct.(fieldName) = theseLines;
            end
            boundaryInfo.velLineElements = velLineStruct;
        end
    case 2
        % For P2 elements, expect LINES3 (Nx4: [nA, nM, nB, flag])
        if ~isfield(msh_vel, 'LINES3') || isempty(msh_vel.LINES3)
            warning('No LINES3 found in velocity mesh => boundaryInfo.velLine3Elements = empty.');
            boundaryInfo.velLine3Elements = struct();
        else
            lines3_vel = msh_vel.LINES3;  % Nx4: [nA, nM, nB, flag]
            uniqueVFlags = unique(lines3_vel(:,4));
            velLine3Struct = struct();
            for iF = 1:numel(uniqueVFlags)
                gFlag = uniqueVFlags(iF);
                if ismember(gFlag, excludedVelFlags)
                    continue;
                end
                mask = (lines3_vel(:,4) == gFlag);
                theseLines = lines3_vel(mask,1:3);  % columns 1-3: node IDs
                fieldName  = sprintf('flag_%d', gFlag);
                velLine3Struct.(fieldName) = theseLines;
            end
            boundaryInfo.velLine3Elements = velLine3Struct;
        end
    case 3
        % For P3 elements, expect LINES4 (Nx5: [nA, interior nodes, nB, flag])
        if ~isfield(msh_vel, 'LINES4') || isempty(msh_vel.LINES4)
            warning('No LINES4 found in velocity mesh => boundaryInfo.velLine4Elements = empty.');
            boundaryInfo.velLine4Elements = struct();
        else
            lines4_vel = msh_vel.LINES4;  % Nx5: [nA, n?, n?, nB, flag]
            uniqueVFlags = unique(lines4_vel(:,5));
            velLine4Struct = struct();
            for iF = 1:numel(uniqueVFlags)
                gFlag = uniqueVFlags(iF);
                if ismember(gFlag, excludedVelFlags)
                    continue;
                end
                mask = (lines4_vel(:,5) == gFlag);
                theseLines = lines4_vel(mask,1:4);  % first 4 columns: node IDs
                fieldName = sprintf('flag_%d', gFlag);
                velLine4Struct.(fieldName) = theseLines;
            end
            boundaryInfo.velLine4Elements = velLine4Struct;
        end
    otherwise
        error('Order must be 1,2, or 3.');
end

% Build a combined list of all boundary nodes.
boundaryInfo.allVelNodes = [];
if order == 1
    if isfield(boundaryInfo, 'velLineElements')
        fnames = fieldnames(boundaryInfo.velLineElements);
        for iF = 1:numel(fnames)
            thisField = fnames{iF};
            lineEls   = boundaryInfo.velLineElements.(thisField);
            if isempty(lineEls), continue; end
            boundaryNodes = unique(lineEls(:));
            boundaryInfo.(thisField) = boundaryNodes;
            boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
        end
    end
elseif order == 2
    if isfield(boundaryInfo, 'velLine3Elements')
        fnames = fieldnames(boundaryInfo.velLine3Elements);
        for iF = 1:numel(fnames)
            thisField = fnames{iF};
            lineEls   = boundaryInfo.velLine3Elements.(thisField);
            if isempty(lineEls), continue; end
            boundaryNodes = unique(lineEls(:));
            boundaryInfo.(thisField) = boundaryNodes;
            boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
        end
    end
elseif order == 3
    if isfield(boundaryInfo, 'velLine4Elements')
        fnames = fieldnames(boundaryInfo.velLine4Elements);
        for iF = 1:numel(fnames)
            thisField = fnames{iF};
            lineEls   = boundaryInfo.velLine4Elements.(thisField);
            if isempty(lineEls), continue; end
            boundaryNodes = unique(lineEls(:));
            boundaryInfo.(thisField) = boundaryNodes;
            boundaryInfo.allVelNodes = [boundaryInfo.allVelNodes; boundaryNodes];
        end
    end
end
boundaryInfo.allVelNodes = unique(boundaryInfo.allVelNodes);

%% (G) (Optional) Pressure boundary conditions can be added similarly.
end
