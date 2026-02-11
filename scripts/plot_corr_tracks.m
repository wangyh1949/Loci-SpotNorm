%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 1/21/2026
    Last update date: 2/11/2026

Description: This script plot the correlation of any two track quantities
(such as amplitude, spotNorm, cellLength), each track as one data point

It makes 1 plots with scatter points for each datasets 

this script is specifically for loci tracking experiments with cell images
---------------------------------------------------------------------------
%}


clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);

% set up strain name info
strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'ori' 'right' 'left' 'lacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};

minTL = 90; % minimal trackLength for plotting
c = 0;

for j = plotNum
        
    c = c + 1;
    
    % load loci file
    load( fullfile( lociList( j).folder, lociList( j).name))
    
    % find strain name
    index = find( strcmp( string( strainList), strain(3:end)));
    if isempty( index), strainName = strain;
    else, strainName = nameList{ index}; end
    
    % set up figures
    figure, set( gcf, "Position", [355*c-330 550 350 300])

    % count spot number for each cell
    cellSpots = accumarray( cellNum(:), 1, [totalCells, 1]);

    % extract spotNorm variables
    xxNorm = tracksxNorm(:,1);   xNorm40 = tracksxNorm40;
    llNorm = tracksLNorm(:,1);   lNorm40 = tracksLNorm40;
    mid = tracksMid(:,1);       mid40 = tracksMid40;
    % tracksMid: [first spot, whole track]

    % extract MSD variables
    msd1 = EnsMSD(:, 1); % ensemble MSD(1)
    tamsd1 = EnsTAMSD(:, 1); % time-averaged MSD(1)

    % extract amplitude variables
    amp = tracksAmp(:,1); % take amplitude at first frame
    % amp = mean( tracksAmp(:,1:40), 2); % average amplitude of first 40 frames

    % ~~~~~~~~ Condition ~~~~~~~~~

        % 1. condition for track length
        tl = tracksLength;    frame = tracksFrame;
        condTracks = ( tl>= minTL) & frame(:,1) == 1;
        % condTrack = ( tl>= minTL);
        
        % 2. cell spot number condition 
        goodCells = find( cellSpots > 0);   extra = '1+ spot'; 
        % goodCells = find( cellSpots == 1);   extra = '1 spot';
        % goodCells = find( cellSpots == 2);   extra = '2 spots';
        condSpots = ismember( cellNum, goodCells); % flag for spots in selected cells
    
        % 3. cap region condition
        condxNorm = all( mid40, 2); % exclude tracks with any cap points

        % ~~~~ final condition ~~~~
        % cond = condTracks & condSpots & condxNorm;
        cond = true( numel( tl), 1); % include all tracks

        fprintf( '   %d/%d cells,   %d/%d tracks    %s\n', ...
            numel( goodCells), totalCells, sum( cond), nTracks, lociName)
        % fprintf( '     condTracks:  %d/%d tracks\n', sum( condTracks), nTracks)
        % fprintf( '     condSpots:   %d/%d tracks\n', sum( condSpots), nTracks)
        % fprintf( '     condxNorm:   %d/%d tracks (no cap points)\n\n', sum( condxNorm), nTracks)

    xPos = xxNorm.* [cellInfo( cellNum).width]'*1e6; % absolute x position, unit: um


        % legtxt = sprintf( '%s (%s)', strain, expDate);
        legtxt = sprintf( '%d tracks', sum( cond));
        legtxt2 = sprintf( '%s %s', strain, strainName);

    % select x and y for plotting
           
       
        % % 1. MSD(1) vs. spotNorm
        % y = tamsd1( cond);          ylabelTxt = 'MSD(1) (µm^2)';
        % x = abs( xxNorm( cond));    xlabelTxt = '|xNorm|';
        % % x = abs( llNorm( cond)-0.5);    xlabelTxt = '|LNorm-0.5|';

        % % 2. MSD(1) vs. amp
        % y = tamsd1( cond); ylabelTxt = 'MSD(1) (µm^2)';
        % x = amp( cond);    xlabelTxt = 'amplitude (a.u.)';

        % % 3. spotNorm vs cellLength
        % x = [ cellInfo( cellNum).length]'*1e6;  xlabelTxt = 'cell length (µm)';
        % y = abs( xxNorm( cond));    ylabelTxt = '|xNorm|';
        % y = abs( llNorm( cond)-0.5);    ylabelTxt = '|LNorm-0.5|';

        % % 4. amp vs. spotNorm (no correlation)
        % y = amp( cond);    ylabelTxt = 'amplitude (a.u.)';
        % x = abs( xxNorm( cond));    xlabelTxt = '|xNorm|';
        % % x = abs( llNorm( cond)-0.5);    xlabelTxt = '|LNorm-0.5|';

        % 5. amp vs. trackLength
        y = amp( cond);    ylabelTxt = 'amplitude (a.u.)';
        x = tl( cond);    xlabelTxt = 'trackLength';


    % plot scatter 
    scatter( x, y, 20, 'filled', 'MarkerFaceAlpha', 0.1, 'DisplayName', legtxt), hold on

    % plot binned mean ± std
    plotBinnedMeanStd( x, y)


    % figure setting
    figure( gcf), set( gca, 'FontSize', 14)
    xlabel( xlabelTxt), ylabel( ylabelTxt)
    legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 12)
    title( legtxt2)
    % xlim( [2 5]), ylim( [0 0.5])
end




%% Function

function plotNum = getPlotNum( list)

    plotNum = 1; % for multiple file of the same strain    
    
    if length( list) > 1
        fprintf( '\n~~~~~ There are %d files for this strain ~~~~~\n', length(list));
        for k = 1:length( list)
            disp( ['  ' num2str(k) '. ' list(k).name])
        end
        tmp = string( input( '\nWhich files do you want to plot (like: 1 3 4): ', 's')); % just input: 1 3 4 5
        % plotNum = double( regexp( tmp, '\d+', 'match')); % find all numbers but not extra space
        plotNum = sscanf( tmp, '%d')'; % parses any valid number format        
    end
    
fprintf( '\n')
end


function plotBinnedMeanStd( x, y)

    % Cap extreme values to suppress outliers
    xCap = min(x, prctile(x, 99));
    xCap = max(xCap, prctile(x, 1)); % also cap low values

    xLim = [min(xCap) max(xCap)];  nBins = 12;
    [~, edges, bin] = histcounts(x, nBins, 'BinLimits', xLim);

    % Compute bin centers
    xCtr = edges(1:end-1) + diff(edges)/2;

    % Select valid data points: inside bins and with non-NaN y values
    valid = ~isnan(y) & bin > 0;
    yMean = accumarray(bin(valid), y(valid), [nBins 1], @mean, NaN);
    yStd  = accumarray(bin(valid), y(valid), [nBins 1], @std,  NaN);

    errorbar(xCtr, yMean, yStd, 'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    hold off
end