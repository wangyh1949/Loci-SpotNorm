%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 1/8/2026
    Last update date: 1/8/2026

    ~~~~~~ adapted from plotLociMSD_ampBin.m ~~~~~~

Description: This script plot loci subpopulation's MSD binned by cell
geometry (length or width).

It makes 2 plots with binned data
    1. EA-MSD       (no fit)
    2. EATA-MSD     (no fit)

this script is specifically for loci tracking experiments with cell images
---------------------------------------------------------------------------
%}


clear, clc, close all

% set the path to store analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);

% set strain name
strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'LacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};


% set up figures
f1 = figure( 'Position', [300 400 400 380]);
f2 = figure( 'Position', [710 400 400 380]);

minTL = 90; % minimal trackLength for plotting

for j = plotNum
        
    % load loci file
    load( fullfile( lociList( j).folder, lociList( j).name))
    
    % find strain name
    index = find( strcmp( string( strainList), strain(3:end)));
    if isempty( index), strainName = strain;
    else, strainName = nameList{ index}; end
    

    % ~~~~~~~~ Condition ~~~~~~~~~
    tl = tracksLength;    frame = tracksFrame;
    cond = ( tl>= minTL) & frame(:,1) == 1;
    % cond = ( tl>= minTL);
    
    fprintf( '   plotting tracks with %d+ frames,  %d/%d tracks    %s\n\n', ...
        minTL, sum( cond), nTracks, lociName)
    
        maxTau = size( EnsMSD, 2);
        time = (1: maxTau)* timeStep;

    % divide data by cellular information
    binData = [ cellInfo.length]';  binName = 'length';
    % binData = [ cellInfo.width]';   binName = 'width';

    binPer = [0 0.3 0.5 0.7 1]; % binning percentage
    binList = quantile( binData, binPer); 
    
    colorList = flip( winter( numel( binPer))); c = 0; 

    for k = 1: numel( binList)- 1
    
        c = c + 1;

        % ~~~~~~~~ Condition ~~~~~~~~~
        goodCell = find( binData >= binList(k) & binData < binList(k+1)); % cell in this bin
        condCell = ismember( cellNum, goodCell); % flag for spots in selected cells
    
        condTracks = cond & condCell; % overall flag for good tracks

        % legtxt = sprintf( '%s (%s)', strain, expDate);
        legtxt = sprintf( '%s: %g%%-%g%%', binName, binPer(k)*100, binPer(k+1)*100);

        fprintf( '      %d/%d tracks,   %s\n', sum( condTracks), sum( cond), legtxt)
        
            eaMSD = mean( EnsMSD( condTracks,:), 1, 'omitnan'); % unit: um^2
            eataMSD = mean( EnsTAMSD( condTracks, :), 1, 'omitnan');    
            

        % plotRange = 1: maxTau;
        plotRange = 1: minTL;
        
        % 1. plot EA-MSD
        figure( f1)
        scatter( time( plotRange), eaMSD( plotRange), 30, 'filled', 'MarkerFaceColor', colorList(c,:), ...
            'DisplayName', legtxt), hold on

        % 2. plot EATA-MSD
        figure( f2)
        scatter( time( plotRange), eataMSD( plotRange), 20, 'filled', 'MarkerFaceColor', colorList(c,:), ...
            'DisplayName', legtxt), hold on
    end                
end


%% Figure Setting

[limX, limY] = findLim( timeStep, strain);

% 1. EA-MSD
figure( f1), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( 'Time (s)'), ylabel( 'EA-MSD (µm^2)')
legend( 'Location', 'southeast', 'box', 'off', 'FontSize', 10)
title( sprintf( '%s-%s %s (%d+f)', expDate, strain, strainName, minTL), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
xlim( limX), ylim( limY)


% 2. EATA-MSD
figure( f2), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( 'Time (s)'), ylabel( 'EATA-MSD (µm^2)')
legend( 'Location', 'southeast', 'box', 'off', 'FontSize', 10)
title( sprintf( '%s-%s %s (EATA)', expDate, strain, strainName), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
xlim( limX), ylim( limY)



%% Function

function plotNum = getPlotNum( list)

    plotNum = 1; % for multiple file of the same strain    
    
    if length( list) > 1
        fprintf( '\n~~~~~ There are %d files for this strain ~~~~~\n', length(list));
        for k = 1:length( list)
            disp( ['  ' num2str(k) '. ' list(k).name])
        end
        tmp = string( input( '\nWhich file do you want to plot:  ', 's')); % just input: 1 3 4 5
        % plotNum = double( regexp( tmp, '\d+', 'match')); % find all numbers but not extra space
        plotNum = sscanf( tmp, '%d')'; % parses any valid number format
        
        plotNum = plotNum(1); % only plot one file
    end
    
fprintf( '\n')
end


function [limX, limY] = findLim( timeStep, strain)

    if timeStep == 1
        limX = [0.8 100]; limY = [5e-3 0.1]; % 100-1000ms
    elseif timeStep == 0.2
        limX = [0.1 20]; limY = [4e-3 0.06]; % 20-200ms
        if strcmp( strain, 'SK311')
            limX = [0.1 30]; limY = [1e-3 0.1]; % 20-200ms, with SK311
        end
    elseif timeStep == 0.02
        limX = [0.01 3]; limY = [3e-4 0.03]; % 20-200ms, with SK311
    else
        limX = 'auto'; limY = 'auto';
    end
end