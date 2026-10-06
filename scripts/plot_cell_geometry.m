%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 7/31/2026
    Last update date: 7/31/2026

~~~~~~~~~~~ adapted from plot_spotNorm.m ~~~~~~~~~~~

Description: this script plots and compares the distributions of cell width
and cell length in selected Loci SpotNorm analysis files. 
 
Every retained cell in cellInfo is included in horizontal raincloud plots
so files with different numbers of cells can be compared directly.
---------------------------------------------------------------------------
%}

clear, clc, close all

% Set the path that stores the analysis-result files.
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm', 'comb');
lociPath = fullfile( varPath, 'Loci SpotNorm', 'rif');

% Plotting settings.
showPoints = false;       % show individual cells as lightly jittered points
maxPointsPerFile = 500;  % plot a random subset when a file has more cells
widthLimits = [];        % e.g. [0.5 1.2] to keep this axis fixed across runs
lengthLimits = [1 6];       % e.g. [1 8] to keep this axis fixed across runs

% Choose the files to compare
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);

    % loadStrainName

strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'LacZ'...'12tetO@LacZ' ... 'LacZ'...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};


% Load geometry from each selected file before plotting.
numFiles = numel( plotNum);
cellWidth = cell( numFiles, 1);
cellLength = cell( numFiles, 1);
labelName = cell( numFiles, 1);
labelName2 = cell( numFiles, 1);

for k = 1: numFiles

    fileIndex = plotNum( k);

    % load lociPos file
    load( fullfile( lociList( fileIndex).folder, lociList( fileIndex).name))

    % find strain name
    index = find( strcmp( string( strainList), strain(3:end)));    
    if isempty( index), strainName = strain;
    else, strainName = nameList{ index}; end

        tInt = findInt( extraName, expDate, strain);
        extraName = erase( extraName, [ ' ' tInt]);
        expDate = erase( expDate, ' comb');

    % data = load( fullfile( lociList( fileIndex).folder, lociList( fileIndex).name), 'cellInfo', 'strain', 'extraName', 'folderName');
    % if ~isfield( data, 'cellInfo') || isempty( data.cellInfo)
    %     error( '"%s" does not contain cellInfo.', lociList( fileIndex).name)
    % end

    % cellInfo dimensions are saved in metres; convert to micrometres.
    cellWidth{ k} = [ cellInfo.width]'* 1e6;
    cellLength{ k} = [ cellInfo.length]'* 1e6;

    % % cellInfo dimensions are saved in metres; convert to micrometres.
    % cellWidth{ k} = [ data.cellInfo.width]' * 1e6;
    % cellLength{ k} = [ data.cellInfo.length]' * 1e6;

    validCells = isfinite( cellWidth{ k}) & isfinite( cellLength{ k}) & ...
        cellWidth{ k} > 0 & cellLength{ k} > 0;
    cellWidth{ k} = cellWidth{ k}( validCells);
    cellLength{ k} = cellLength{ k}( validCells);


    legtxt = sprintf( '%s%s', strain, extraName);
    legtxt2 = sprintf( '%s%s', strainName, extraName);

    labelName{ k} = legtxt;
    labelName2{ k} = legtxt2; % makeLabel( data, lociList( fileIndex).name);
    
    fprintf( '  %s: %d cells, width = %.3f +/- %.3f um, length = %.3f +/- %.3f um\n', ...
        labelName{ k}, numel( cellWidth{ k}), mean( cellWidth{ k}), std( cellWidth{ k}), ...
        mean( cellLength{ k}), std( cellLength{ k}))
end

colorList = lines( numFiles);
%%
% Each row is one file: half-violin = density, box = median/IQR, and dots
% = individual cells. The horizontal measurement axis never depends on the
% number of selected files.
figure( 'Position', [400 400 800 110+75*numFiles])
tiledlayout( 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
plotGeometryRaincloud( cellWidth, labelName, colorList, showPoints, ...
    maxPointsPerFile, widthLimits, 'Cell Width (μm)', 'Cell Wid (tracking)')

nexttile
plotGeometryRaincloud( cellLength, labelName2, colorList, showPoints, ...
    maxPointsPerFile, lengthLimits, 'Cell Length (μm)', 'Cell Length (tracking)')


% figure( 'Position', [400 400 420 max(320, 110+ 75*numFiles)])
% plotGeometryRaincloud( cellWidth, labelName, colorList, showPoints, ...
%     maxPointsPerFile, widthLimits, 'Cell Width (μm)', 'Cell Width')
% 
% figure( 'Position', [830 400 420 max(320, 110+ 75*numFiles)])
% plotGeometryRaincloud( cellLength, labelName2, colorList, showPoints, ...
%     maxPointsPerFile, lengthLimits, 'Cell Length (μm)', 'Cell Length')



function plotGeometryRaincloud( geometry, labelName, colorList, showPoints, ...
    maxPoints, xLimits, xLabel, plotTitle)
    numFiles = numel( geometry);
    hold on

    for k = 1: numFiles
        % With YDir reversed, row 1 is displayed at the top.
        drawRaincloud( geometry{ k}, k, colorList( k,:), showPoints, maxPoints)
    end

    set( gca, 'FontName', 'Arial', 'FontSize', 14, 'LineWidth', 1, ...
        'YTick', 1:numFiles, 'YTickLabel', labelName, 'YDir', 'reverse')
    ylim( [0.4 numFiles+ 0.25])
    if ~isempty( xLimits)
        xlim( xLimits)
    end
    xlabel( xLabel)
    title( plotTitle)
    grid on, set(gca, 'GridLineWidth', 0.5)
    box off
end


function drawRaincloud( values, row, color, showPoints, maxPoints)

    % Draw a half violin above the row. Its height encodes density only;
    % horizontal position always encodes the cell measurement.
    if numel( values) > 1
        [density, x] = ksdensity( values);
        violinHeight = 0.45 * density / max( density);
        fill( [x fliplr( x)], [row - violinHeight, row + zeros( size( x))], color, ...
            'FaceAlpha', 0.25, 'EdgeColor', color, 'LineWidth', 1)
    end

    % Draw a compact horizontal Tukey box plot (median, IQR, and whiskers).
    quartiles = quantile( values, [0.25 0.5 0.75]);
    iqrValue = quartiles(3) - quartiles(1);
    inlierValues = values( values >= quartiles(1) - 1.5*iqrValue & ...
        values <= quartiles(3) + 1.5*iqrValue);
    whiskers = [min( inlierValues), max( inlierValues)];
    boxHeight = 0.2;
    boxBottom = row - boxHeight/2;

    line( whiskers, [row row], 'Color', color, 'LineWidth', 1.2)
    line( whiskers(1)*[1 1], [row-0.07 row+0.07], 'Color', color, 'LineWidth', 1.2)
    line( whiskers(2)*[1 1], [row-0.07 row+0.07], 'Color', color, 'LineWidth', 1.2)
    rectangle( 'Position', [quartiles(1), boxBottom, iqrValue, boxHeight], ...
        'EdgeColor', color, 'LineWidth', 1, 'FaceColor', color, 'FaceAlpha', 0.18)
    line( quartiles(2)*[1 1], [boxBottom boxBottom+boxHeight], 'Color', color, 'LineWidth', 1.5)

    if showPoints
        if numel( values) > maxPoints
            pointIndex = randperm( numel( values), maxPoints);
            values = values( pointIndex);
        end
        yJitter = row + 0.12 + 0.12*rand( size( values));
        scatter( values, yJitter, 9, color, 'filled', ...
            'MarkerFaceAlpha', 0.18, 'MarkerEdgeAlpha', 0.18)
    end
end
