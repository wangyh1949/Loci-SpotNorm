%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 1/12/2026
    Last update date: 1/12/2026

    ~~~~~~~ adapted from plot_MSDxy_cell.m ~~~~~~~

Description: This code calculates the MSD along cellular long & short axes
    for loci data that have been spot-normalized.

It makes one figure for each dataset
    1. MSD along cell short axis (x)
    2. MSD along cell long axis (L/y)

---------------------------------------------------------------------------
%}

%% 1. 

clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';

lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
tfPath   = fullfile( lociPath, 'tracksFinal'); % subfolder under lociPath

% find which files to plot
tfFile = dir( fullfile( tfPath, ['tf oufti*']));
plotNum = getPlotNum( tfFile);


% set up strain name info
strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'ori' 'right' 'left' 'lacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};


minTL = 90; % minimal trackLength for MSD calculation
c = 0;

for i = plotNum

    c = c + 1;
    % load the trackFinal file
    load( fullfile( tfPath, tfFile(i).name)),  fprintf( '    ~~~ %s loaded ~~~', tfName)
    
        % find strain name
        index = find( strcmp( string( strainList), strain(3:end)));    
        if isempty( index), strainName = strain;
        else, strainName = nameList{ index}; end
        
    nTracks = numel( tracksFinal);

    tl = cellfun( @length, {tracksFinal.amp}');
    cond = ( tl >= minTL);

    maxT = max( tl); plotRange = 1: minTL;

    % preallocate for MSD
    EnsMSDx = nan( nTracks, maxT-1);    EnsTAMSDx = nan( nTracks, maxT-1);
    EnsMSDy = nan( nTracks, maxT-1);    EnsTAMSDy = nan( nTracks, maxT-1);

    for j = find( cond)'
        
        % get spot positions in cell coordinates
        xPos = tracksFinal(j).spotPosInCell(:,2)* pixelSize* 1e6; % unit: um
        LPos = tracksFinal(j).spotPosInCell(:,1)* pixelSize* 1e6; % unit: um
        
        [ EnsMSDx( j,:), EnsTAMSDx( j,:)] = getMSD( xPos, maxT);
        [ EnsMSDy( j,:), EnsTAMSDy( j,:)] = getMSD( LPos, maxT);

    end

    eataMSDx = mean( EnsTAMSDx( cond, :), 1, 'omitnan'); % unit: um^2
    eataMSDy = mean( EnsTAMSDy( cond, :), 1, 'omitnan'); % unit: um^2
    eaMSDx = mean( EnsMSDx( cond,:), 1, 'omitnan');
    eaMSDy = mean( EnsMSDy( cond,:), 1, 'omitnan');


    % determine frame time (timeStep) for the movie
    frameT = str2double( erase( regexp( extraName, '\d+ms', 'match', 'once'), 'ms'));
    timeStep = frameT* 1e-3; % unit: s    
    fprintf( '        ~~~~ frame time is  %d ms ~~~~\n', timeStep*1e3)
    time = (1: maxT-1)* timeStep;


    % set legend text
    legtxt = sprintf( '%s, %s', strain, expDate); 
    legtxt2 = sprintf( '%s,%s', strainName, extraName);

    % set up figures
    figure, set( gcf, "Position", [410*c-110 550 400 270])

    % plot EATA-MSD for x Pos
    scatter( time( plotRange), eataMSDx( plotRange), 20, 'filled', 'DisplayName', 'transverse x'), hold on

    % plot EATA-MSD for L/y Pos
    scatter( time( plotRange), eataMSDy( plotRange), 20, 'filled', 'DisplayName', 'longitudinal L')

    % figure setting
    [limX, limY] = findLim( timeStep, strain);

    % 1. EATA-MSDx
    figure( gcf), set( gca, 'LineWidth', 1, 'FontSize', 14)
    xlabel( 'Time (s)'), ylabel( 'EATA-MSD (µm^2)')
    legend( 'Location', 'southeast', 'box', 'off', 'FontSize', 12)
    title( sprintf( '%s %s (%s, %d+f)', strain, strainName, expDate, minTL), 'FontSize', 14)
    set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
    xlim( limX), ylim( limY)

end



%% function

function plotNum = getPlotNum( list)

    plotNum = 1; % for multiple file of the same strain    
    
    if length( list) > 1
        fprintf( '\n~~~~~ There are %d files for this strain ~~~~~\n', length(list));
        for k = 1:length( list)
            disp( ['  ' num2str(k) '. ' list(k).name])
        end
        tmp = string( input( '\nWhich files do you want to plot (like: 1 3 4): ', 's')); % just input: 1 3 4 5
        plotNum = sscanf( tmp, '%d')'; % parses any valid number format
    end
    
fprintf( '\n')
end


function [ EnsMSD, EnsTAMSD] = getMSD( traj, maxT)

    nFrames = numel( traj);
    
    EnsMSD = nan(1, maxT-1);
    EnsTAMSD = nan(1, maxT-1);

    % ~~~~ Ensemble-Averaged MSD ~~~~~
    dr = traj( 2: min( nFrames, maxT), :) - traj( 1, :);
    EnsMSD( 1: numel( dr)) = sum( dr.^2, 2); % ensemble MSD for EA-MSD plotting        

    % ~~~~ Time-Averaged MSD ~~~~~
    for tau = 1: min( nFrames, maxT)-1
        dr = traj( tau+1:end, :) - traj( 1:end-tau, :);
        dr2 = sum( dr.^2, 2);
        EnsTAMSD(tau) = mean( dr2, 'omitnan'); % ensemble TA-MSD for EATA-MSD plotting
    end    
end


function [limX, limY] = findLim( timeStep, strain)

    if timeStep == 1
        limX = [0.8 100]; limY = [5e-3 0.1]; % 100-1000ms
    elseif timeStep == 0.2
        limX = [0.1 20]; limY = [1e-3 0.03]; % 20-200ms
        if strcmp( strain, 'SK311')
            limX = [0.1 30]; limY = [1e-3 0.1]; % 20-200ms, with SK311
        end
    elseif timeStep == 0.02
        limX = [0.01 3]; limY = [3e-4 0.03]; % 20-200ms, with SK311
    else
        limX = 'auto'; limY = 'auto';
    end
end