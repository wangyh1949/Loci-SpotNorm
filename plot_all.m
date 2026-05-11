%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 1/29/2026
    Last update date: 1/29/2026

~~~~~~~~~~~ adapted from plotLoci.m ~~~~~~~~~~~

Description: This code make plots for loci tracking data. It plots basic
properties:
        1. track length distribution (1-CDF)
        2. EA/EATA MSD & (optional) linear fitting
        3. Photobleaching curve
        4. D_alpha histogram
        5. alpha histogram
        6. locErr histogram (posStd from uTrack spot detection)

(optional) redo the linear fitting of individual TA-MSD with a different
fitting range, default range is 25%
---------------------------------------------------------------------------
%}

clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
% lociPath = fullfile( varPath, 'Loci SpotNorm', 'single day'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);

    strainList = { 727 728 729 730 731 734 701 662 311 725};
    nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'LacZ'... %'12tetO@LacZ' ...
        '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};


%% Load data and make plots

% plotting condition
minTL = 90; % minimal trackLength for MSD calculation

% plotting parameters
dAlphaBin = 4e-4;
fitFlag = true; % flag: show MSD fitting

% fitting range
fitR = 1:40;    fitR_nl = fitR;
fitTxt1 = sprintf( '%d:%d fit', min( fitR), max( fitR));
fitTxt_nl = sprintf( '%d:%d fit', min( fitR_nl), max( fitR_nl));

% set up figures
f1 = figure( 'Position', [200 650 400 380]);
f2 = figure( 'Position', [610 650 400 380]);
f3 = figure( 'Position', [1020 650 400 380]);
f4 = figure( 'Position', [200 185 400 380]);
f5 = figure( 'Position', [610 185 400 380]);
f6 = figure( 'Position', [1020 185 400 380]);
colorList = get( gca,'colororder');
colorList = repmat( colorList, [2, 1]); c = 0;
% colorList = flip( winter( numel( plotNum)));


for j = plotNum
    
    c = c + 1;
    % load loci file
    load( [ lociList( j).folder '\' lociList( j).name])
    
        % find strain name
        index = find( strcmp( string( strainList), strain(3:end)));    
        if isempty( index), strainName = strain;
        else, strainName = nameList{ index}; end
    
        tInt = findInt( extraName, expDate, strain);    
        extraName = erase( extraName, [ ' ' tInt]);
        expDate = erase( expDate, ' comb');
        expT = str2double( regexp( tInt, '\d+', 'match', 'once'))/1e3; % exposure time: s 


    % set up legend text
    legtxt = sprintf( '%s%s', strainName, extraName);
    % legtxt = sprintf( '%s%s', strain, extraName);
    % legtxt = sprintf( '%s, %s', strain, expDate);
    legtxt3 = sprintf( '%s %s%s', strain, strainName, extraName);

    % ~~~~~~~~ Condition ~~~~~~~~~
    tl = tracksLength;
    frame = tracksFrame;
    cond = ( tl>= minTL) & frame(:,1) == 1;
    % cond = ( tl>= minTL);
    fprintf( ' plotting tracks with %d+ frames,  %d tracks    %s\n', minTL, sum( cond), lociName)  
   

    % 1. tracksLength 1-CDF
    figure( f1)

    % tl_good = tl;
    tl_good = tl( frame(:,1) == 1); % only consider tracks start from frame 1

        [f,x] = ecdf( tl_good); x = x(2:end); f = f(1:end-1);
        plot( x, 1-f, 'LineWidth', 2, 'Color', colorList(c,:), 'DisplayName',...
            sprintf( '%d tracks [%s]', length(tl_good), legtxt)), hold on
    
    
    % 2. EA-MSD & fitting
    figure( f2)

    % calculate EA-MSD & EATA-MSD
    eaMSD = mean( EnsMSD( cond,:), 1, 'omitnan'); % unit: um^2
    eataMSD = mean( EnsTAMSD( cond, :), 1, 'omitnan');
    msd = eaMSD;    msdName = 'EA-MSD';
    % msd = eataMSD;    msdName = 'EATA-MSD';
    
        % determine time step
        maxTau = size( EnsMSD, 2);
        time = (1: maxTau)* timeStep;
            
        if ~logical( fitFlag) % no fitting
            scatter( time, msd, 20, 'filled', 'MarkerFaceColor', colorList(c,:),...
                'MarkerFaceAlpha', 0.5, 'DisplayName',...
                sprintf( '%s, %d tracks', legtxt, sum( cond))), hold on
            fitTxt1 = 'no fit';
        else % with fitting
            scatter( time, msd, 20, 'filled', 'MarkerFaceColor', colorList(c,:),...
                'MarkerFaceAlpha', 0.5, 'HandleVisibility','off'), hold on
    
            % linear fitting,  MSD = 4*D*t^alpha        
            f = polyfit( log( time( fitR)), log( msd( fitR)), 1);
            alphaFit = f(1);    DFit = exp( f(2))/4;
            t = time( fitR);    MSDFit =  4* DFit* t.^ alphaFit;
            plot( t, MSDFit, 'LineWidth', 1, 'color', colorList(c,:), ...
                'DisplayName', sprintf( 'D\\alpha=%.1e, \\alpha=%.2f [%s]', DFit, alphaFit, legtxt))
                % 'DisplayName', sprintf( '\\alpha=%.2f  [%s]', alphaFit, legtxt))
        end
    

    % 3. amplitude over frame
    figure( f3)

    % calculate average amplitude of all tracks
    ampAvg = mean( tracksAmp( cond, :), 'omitnan');%* sigToPhoton;
    amp = ampAvg;
    % amp = movmean( ampAvg, 5); % make it more smooth

        tFrame = 1: minTL; % frame N
        plot( tFrame, amp(tFrame), '-', 'LineWidth', 2, 'color', colorList(c,:), ...
            'DisplayName', legtxt), hold on
    

        
    % (optional) fit Individual TA_MSD to get D & alpha
    fitRange = 1:40;
    fitTxt = sprintf( '%d:%d fit', min( fitRange), max( fitRange));
    alpha = nan( nTracks, 1);   Dalpha = nan( nTracks, 1);
    time = fitRange* timeStep;

    for k = 1: size( EnsTAMSD, 1)
        f = polyfit( log( time), log( EnsTAMSD( k, fitRange)), 1);
        alpha(k) = f(1);    Dalpha(k) = exp( f(2))/4;
    end
    

    % 4. D_alpha histogram
    figure( f4)
    dAlpha = Dalpha( cond);  meanD = mean( dAlpha);
    histogram( dAlpha, 'BinWidth', dAlphaBin, 'LineWidth', 2, 'EdgeColor', colorList(c,:), ...
        'Normalization', 'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
        sprintf( '\\langleD\\rangle=%.1e  [%s]', meanD, legtxt)), hold on
    
    
    % 5. alpha histogram
    figure( f5)
    aAlpha = alpha( cond);  meanAlpha = mean( aAlpha);
    histogram( aAlpha, 'BinWidth', 0.04, 'LineWidth', 2, 'EdgeColor', colorList(c,:),...
        'Normalization', 'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
        sprintf( '\\langle\\alpha\\rangle=%.2f\\pm%.2f  [%s]', meanAlpha, std(aAlpha), legtxt)), hold on
    % xline( alphaFit, '--', 'color', colorList(c,:), 'LineWidth', 1.5, 'HandleVisibility', 'off')
    
    
    % 6. Localization Error Histogram
    figure( f6)
    lLocErr = abs( posStd( cond));     meanLocErr = mean( lLocErr, 'omitnan');
    histogram( lLocErr, 'BinWidth', 2, 'LineWidth', 2,  'EdgeColor', colorList(c,:), ...
        'Normalization', 'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
        sprintf( '\\langle\\sigma\\rangle=%.0f nm  [%s]', meanLocErr, legtxt))
    hold on        
end

fprintf( '\n~~~~~~ !!!  All Complete  !!! ~~~~~~\n')


%% Figure Setting

% trackLength distribution
figure( f1), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( 'Track Length (frames)'), ylabel( '1-CDF (portion)')
legend( 'Location', 'best', 'box', 'off', 'FontSize', 11)
title( 'Track Length Distribution', 'FontSize', 14)
% xline(12, '--', 'LineWidth', 1.5, 'FontSize', 13, 'HandleVisibility', 'off')
ylim( [0 1])


% MSD & fitting
figure( f2), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( 'Time (s)'), ylabel( 'MSD (µm^2)')
title( sprintf( '%s (%d+ frames, %s)', msdName, minTL, fitTxt1), 'FontSize', 14)
legend( 'Location', 'northwest', 'box', 'off', 'FontSize', 10)
set( gca, 'Xscale', 'log', 'YScale', 'log')

[limX, limY] = findLim( timeStep, strain);
xlim( limX), ylim( limY), box on


% Amplitude over time
figure( f3), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( 'Frame'), ylabel( 'Signal Intensity (a.u.)')
title( sprintf( 'Photobleaching  (%d+ frames)', minTL), 'FontSize', 14)
legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 12)
   

% D_alpha
figure( f4), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( 'D\alpha (\mum^2/s)'), ylabel( 'Probability')
title( sprintf( 'D\\alpha (%d+ frames, %s)', minTL, fitTxt), 'FontSize', 14)
legend('Location', 'northeast', 'box', 'off', 'FontSize', 12)
% xlim([-0.001 0.01])


% alpha
figure( f5), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( sprintf( 'alpha')), ylabel( 'Probability')
title( sprintf( '\\alpha (%d+ frames, %s)', minTL, fitTxt), 'FontSize', 14)
legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 12)
ylim( [0 0.25])


% locErr
figure( f6), set( gca, 'LineWidth', 1, 'FontSize', 14)
xlabel( 'Localization Error \sigma (nm)'), ylabel( 'Probability')
title( sprintf( 'locErr (from uTrack)'), 'FontSize', 14)
legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 12)



%% Function

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