%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 7/25/2024
    Last update date: 12/21/2025

    ~~~~~~ adapted from meshCleanup.m ~~~~~~

this script is for cleaning the mesh data, delete problematic cells after
oufti automatic cell outline detection

this version is specifically for phase images taken with SPT lens
(loci tracking experiments with phase images)

pixel size:     64.5 nm (phase) vs 160 nm (SPT) 
---------------------------------------------------------------------------
%}

clear, clc, close all

% load oufti mesh output (cell meshes)
load( 'mesh.mat')


% elimination criteria
minArea = 100;  % minimal cell area (pixel)
maxArea = 200;  % maximal cell area (pixel)
minWid = 6;     % minimal cell width (pixel)
maxWid = 8.3;   % maximal cell width (pixel)

    
% set up cell mesh and remove bad cells
nImg = size( cellListN, 2);
cellAll = cell( nImg, 1);

for imgNum = 1: nImg % number of images

    cellMesh = struct();
    meshData = cellList.meshData{ imgNum}';
    cellid = double( cellList.cellId{ imgNum});

    nCell = size( meshData, 1); % number of cells

    if nCell == 0  % in case no cell are detected
        warning( '~~~~~~ Image %2d have no cell mesh ~~~~~~\n', imgNum)
        continue
    end
    
    cellMesh( nCell, 1) = struct();
    badCells = false( 1, nCell);
    
    for cellNum = 1: nCell
        
        mesh = meshData{ cellNum}.mesh; % cell outline: [n,4] - (x1, y1, x2, y2), (right, left)
        
        % find empty cell mesh
        if length( mesh) < 10
            badCells( cellNum) = true;
            fprintf( '~~~~~~ Image %2d Cell #%-3d have no mesh ~~~~~~\n', imgNum, cellid( cellNum))
            continue
        end
            
        % calculate cell mesh properties for later use
        meshOut = [ mesh(:, 1:2); flipud( mesh(:, 3:4))]; % reshape the mesh matrix to form a circle [2n, 2]
        cellMesh( cellNum).area = double( polyarea( meshOut(:,1), meshOut(:,2)));

        meshMid = [ mean( mesh( :, [1 3]), 2), mean( mesh( :, [2 4]), 2)]; % midline along the long axis
        gridLen = vecnorm( diff( meshMid), 2, 2); % length of each grid (L direction)
        cellMesh( cellNum).length = double( sum( gridLen));
            
        % find cells that are too small
        if ( cellMesh( cellNum).area < minArea)
            badCells( cellNum) = true;
            fprintf( '~~~~~~ Image %2d Cell #%-3d is too small ~~~~~~\n', imgNum, cellid( cellNum))
            continue
        end

        % find cells that are too big
        if ( cellMesh( cellNum).area > maxArea)
            badCells( cellNum) = true;
            fprintf( '~~~~~~ Image %2d Cell #%-3d is too big ~~~~~~\n', imgNum, cellid( cellNum))
            continue
        end

        cellWid = cellMesh( cellNum).area/ cellMesh( cellNum).length;
        if cellWid < minWid || cellWid > maxWid
            badCells( cellNum) = true;
            fprintf( '~~~~~~ Image %2d Cell #%-3d has abnormal width ~~~\n', imgNum, cellid( cellNum))
            continue
        end
        
        % if imgNum == 2 && cellList.cellId{ imgNum}(cellNum) == 73
        %     a = 1;
        % end

        % find curled cells
        movNum = 3;
        meshTest = movmean( meshOut( [1:end 1:movNum], :), movNum, 'Endpoints', 'discard'); % smooth cell mesh
        for k = 1: length( meshTest)
            pt = meshTest( k, :);
            % find distance to each segment (short axis)
            [~, dist] = findPerpFoot( pt, mesh(:, 3:4), mesh(:, 1:2)); % -: below, +: above, should be - to +
            bra = find( abs( diff( sign( dist))) == 2); % find the grid idx where distance changes sign, first & end = nan
            if length( bra) > 1 % bad cell, outline curls back
                badCells( cellNum) = true;
                fprintf( '~~~~~~ Image %2d Cell #%-3d is curled ~~~~~~\n', imgNum, cellid( cellNum))
                break
            end
        end

        if cellNum > length( cellList.cellId{ imgNum}) % sometimes oufti has this error
            cellMesh( cellNum).cellId = nan;
            fprintf( '~~~~~~ Image %2d Cell #%-3d has problematic cellID ! ~~~~~~\n', imgNum, cellid( cellNum))
            continue
        end
        cellMesh( cellNum).cellId = [ imgNum cellid( cellNum)];            
    end
    
    cellList.meshData{ imgNum}( badCells) = [];
    cellList.cellId{ imgNum}( badCells) = [];
    
    cellMesh( badCells) = [];
    cellAll{ imgNum} = cellMesh;
end

% recalculate the cell number in each image
cellListN = cellfun( @length, cellList.cellId);         
    
% save as mesha.mat file
save( 'mesha', 'p', 'rawPhaseFolder', 'cellList', 'cellListN', 'paramString', 'cellAll')

fprintf( '\n~~~~~~ all images cleanup done ~~~~~~\n\n')

cellAll = vertcat( cellAll{:});

% display the cell length & wid
cellArea = vertcat( cellAll.area);
cellLength = vertcat( cellAll.length);
cellWid = cellArea./ cellLength;
fprintf( '  ~~~ length: %.2f pix,  wid: %.2f pix,  area: %.2f pix,  %d cells  ~~~ \n',...
    mean( cellLength), mean( cellWid), mean( cellArea), sum( cellListN)) 


%% Function

function [D, dist] = findPerpFoot( pt, B, C)
% this function returns the coordinate of the perpendicular foot D so that 
% AD perpendicular to BC (everything in 2D) and the distance of pt to line BC
% dist > 0 if pt is on the right side of line BC (pt, B, C: counterclockwise)
% dist < 0 if pt is on the  left side of line BC (pt, B, C: clockwise)

    AB = B - pt; % vector
    BC = C - B;  % vector

    area = AB(:,1).*BC(:,2) - AB(:,2).*BC(:,1); % cross product
    side = sign( area); % -1: left side, +1: right side

    normVec = [ BC(:,2) -BC(:,1)]; % normal vector of BC in 2D 
    unitNormVec = normVec./ vecnorm( normVec, 2, 2); % unit normal vector of BC
    AD = dot(unitNormVec, AB, 2).* unitNormVec; % AD is perpendicular to BC, dot product
    D = pt + AD; % D point of intersection, Perpendicular Foot
    dist = side.* vecnorm( AD, 2, 2);
end    
