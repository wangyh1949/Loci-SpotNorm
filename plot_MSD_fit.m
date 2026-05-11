%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 1/12/2026
    Last update date: 1/12/2026

~~~~~~~~~~~ adapted from plotLociMSD.m ~~~~~~~~~~~

Description: This script plots MSD of loci tracking data and does fitting.

        1. EA-MSD & linear fitting
        2. EA-MSD & non-linear fitting
        3. EATA-MSD & linear fitting
        4. EATA-MSD & non-linear fitting

---------------------------------------------------------------------------
%}

clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));

% lociPath = fullfile( varPath, 'Loci', '20ms'); % subfolder under varPath
% % find which files to plot
% lociList = dir( fullfile( lociPath, 'Loci*'));

plotNum = getPlotNum( lociList);

    strainList = { 727 728 729 730 731 734 701 662 311 725};
    nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'LacZ'... %'12tetO@LacZ' ...
        '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};

        
%% Load data and Plot MSD

% plotting condition
minTL = 90; % minimal trackLength for MSD calculation

% fitting range
fitR = 1:90;    fitR_nl = fitR;
fitTxt1 = sprintf( '%d:%d', min( fitR), max( fitR));
fitTxt_nl = sprintf( '%d:%d', min( fitR_nl), max( fitR_nl));

% set up figures
f1 = figure( 'Position', [300 400 400 380]);
f2 = figure( 'Position', [710 400 400 380]);
f3 = figure( 'Position', [1120 400 400 380]);
f4 = figure( 'Position', [1530 400 400 380]);
colorList = get( gca,'colororder');     
colorList = repmat( colorList, [2, 1]); c = 0;

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


    % set up legend texts
    legtxt = sprintf( '%s%s', strainName, extraName);
    legtxt2 = sprintf( '%s, %s', strain, expDate);
    

    % ~~~~~~~~ Condition ~~~~~~~~~
    tl = tracksLength;
    frame = tracksFrame;
    cond = ( tl>= minTL) & frame(:,1) == 1;
    % cond = ( tl>= minTL);
    fprintf( ' plotting tracks with %d+ frames,  %d tracks    %s\n', minTL, sum( cond), lociName)  


    % calculate EA-MSD & EATA-MSD
    eaMSD = mean( EnsMSD( cond,:), 1, 'omitnan'); % unit: um^2
    eataMSD = mean( EnsTAMSD( cond, :), 1, 'omitnan');

        % determine time step 
        time = (1: maxT-1)* timeStep;
        plotRange = 1: minTL;
        % plotRange = 1: maxT-1;

    
    % 1. plot EA-MSD & linear fit
    figure( f1)
    scatter( time( plotRange), eaMSD( plotRange), 20, 'filled', 'MarkerFaceColor', colorList(c,:),...
        'MarkerFaceAlpha', 0.7, 'HandleVisibility','off'), hold on
        
        % linear fit
        f = polyfit( log( time( fitR)), log( eaMSD( fitR)), 1);
        alphaFit2 = f(1);    Dapp2 = exp( f(2))/4;
        t = time( fitR);  MSDFit =  4 * Dapp2* t.^ alphaFit2;
        plot( t, MSDFit, 'LineWidth', 1, 'color', colorList(c,:), ...
            'DisplayName', sprintf( 'D\\alpha=%.1e, \\alpha=%.2f [%s]', Dapp2, alphaFit2, legtxt2))


    % 2. plot EA-MSD & non-linear fit
    figure( f2)
    scatter( time( plotRange), eaMSD( plotRange), 20, 'filled', 'MarkerFaceColor', colorList(c,:),...
        'MarkerFaceAlpha', 0.7, 'HandleVisibility','off'), hold on

        % non-linear fit
        [DFit, alphaFit, locErrFit, f] = nonLinearMSDFit( time, eaMSD, fitR_nl, expT);
        t =  linspace( time(1), time( max( minTL, fitR_nl(end))), 100);
        MSDFit =  4* DFit* t.^ alphaFit + 4*f.c;
        plot( t, MSDFit, 'LineWidth', 1, 'color', [ colorList(c,:) ], 'DisplayName', ...
            sprintf( 'D=%.1e, \\alpha=%.2f, \\sigma=%.0f [%s]', DFit, alphaFit, locErrFit*1e3, legtxt))
        fprintf( '     EA fit:  f.c = %f,   locE = %g\n', f.c, locErrFit*1e3)
    

    % 3. plot EATA-MSD & linear fit
    figure( f3)
    scatter( time( plotRange), eataMSD( plotRange), 20, 'filled', ...
        'MarkerFaceColor', colorList(c,:), 'HandleVisibility', 'off'), hold on
        
        % linear fit
        f = polyfit( log( time( fitR)), log( eataMSD( fitR)), 1);
        alphaFit = f(1);    Dapp = exp( f(2))/4; 
        t = time( fitR);   MSDFit =  4* Dapp* t.^ alphaFit;
        plot( t, MSDFit, 'LineWidth', 1, 'color', [ colorList(c,:) ], ...
            'DisplayName', sprintf( '\\alpha=%.2f [%s]', alphaFit, legtxt))


    % 4. plot EATA-MSD & non-linear fit
    figure( f4)
    scatter( time( plotRange), eataMSD( plotRange), 20, 'filled', ...
        'MarkerFaceColor', colorList(c,:), 'HandleVisibility', 'off'), hold on

        % non-linear fit
        [DFit, alphaFit, locErrFit, f] = nonLinearMSDFit( time, eataMSD, fitR_nl, expT);
        t =  linspace( time(1), time( max( minTL, fitR_nl(end))), 100); % extend to 100f
        t =  linspace( time(1), time( fitR_nl(end))); % fit region
        MSDFit =  4* DFit* t.^ alphaFit + 4*f.c;
        plot( t, MSDFit, 'LineWidth', 1, 'color', [ colorList(c,:) ], 'DisplayName', ...
            sprintf( '\\alpha=%.2f, \\sigma=%.0f [%s]', alphaFit, locErrFit*1e3, legtxt))
        fprintf( '   EATA fit:  f.c = %f,   locE = %g\n\n', f.c, locErrFit*1e3)


    % % display the cell numbers & spot numbers
    % fprintf( '  ~~~ %3d images, %4d cells, %5d/%5d tracks,    %s-%s%s\n\n',...
    %     size( cellRecord, 1), totalCells, sum( cond), nTracks, expDate, strain, extraName)    
end

%% figure setting

[limX, limY] = findLim( timeStep, strain);

% 1. EA-MSD & linear fitting
figure( f1), set( gca, 'LineWidth', 1, 'FontSize', 13)
xlabel( 'Time (s)'), ylabel( 'EA-MSD (µm^2)')
legend( 'Location', 'northwest', 'box', 'off', 'FontSize', 10)
title( sprintf( 'EA-MSD (%d+frame, %s fit)', minTL, fitTxt1), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
xlim( limX), ylim( limY)

% 3. EATA-MSD & linear fitting
figure( f3), set( gca, 'LineWidth', 1, 'FontSize', 13)
% set( gca, 'Xtick', [1 10 100])
xlabel( 'Time (s)'), ylabel( 'EATA-MSD (µm^2)') % ylabel( 'EATA-MSD (µm^2)')
legend( 'Location', 'northwest', 'box', 'off', 'FontSize', 10)
title( sprintf( '%s linear fit (%d+frame)', fitTxt_nl, minTL), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
xlim( limX), ylim( limY)

% 2. EA-MSD & non-linear fitting
figure( f2), set( gca, 'LineWidth', 1, 'FontSize', 13)
xlabel( 'Time (s)'), ylabel( 'EA-MSD (µm^2)')
legend( 'Location', 'northwest', 'box', 'off', 'FontSize', 10)
title( sprintf( 'EA-MSD (%d+frame, %s)', minTL, fitTxt_nl), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
xlim( limX), ylim( limY)

% 4. EATA-MSD & non-inear fitting
figure( f4), set( gca, 'LineWidth', 1, 'FontSize', 13)
xlabel( 'Time (s)'), ylabel( 'EATA-MSD (µm^2)')
legend( 'Location', 'northwest', 'box', 'off', 'FontSize', 10)
title( sprintf( '%s non-linear fit (%d+frame)', fitTxt_nl, minTL), 'FontSize', 14)
set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
xlim( limX), ylim( limY)



%% Function

function [ DFit, alphaFit, locErrFit, f] = nonLinearMSDFit( time, MSD, fitR, expT)

    % Linear fit for initial estimation
    f = polyfit( log( time( fitR)), log( MSD( fitR)), 1);
    alphaFit = f(1);    DFit = exp( f(2))/4;    locErrFit = 20; % nm

    % non-linear fit using a comprehensive form, MSD = 4Dt^a + b
    fun = fittype( 'log( 4*a*(x)^b+4*c)');
    x0 = [ DFit, alphaFit, (locErrFit/1e3)^2];
    xmin = [ 0, 0, -2e-3];   xmax = [ 0.01, 1, 2e-3];
    % xmin = [ 0, 0, -inf];   xmax = [ inf, 1, inf];

    [ f, ~] = fit( ( time( fitR)'), log( MSD( fitR)'), fun, 'StartPoint', x0, 'Lower', xmin, 'Upper', xmax);
    DFit = f.a;  alphaFit = f.b;  %locErrFit = sign( f.c)* sqrt( abs( f.c)); % unit: nm
    
    motionBlur = 4*2* DFit* expT^alphaFit/ (( 1+alphaFit)* (2+alphaFit));
    loc2 = ( 4*f.c + motionBlur)/ 4;
    locErrFit = sign( loc2)* sqrt( abs( 4*f.c + motionBlur)/ 4); % unit: um

    if abs( f.c) >= abs( xmax(3))- 1e-6
        warning(' Parameter locErr has hit its bound %.4g. Bound may be too tight.', xmax(3))
    end
end


function [limX, limY] = findLim( timeStep, strain)

    if timeStep == 1
        limX = [0.8 100]; limY = [5e-3 0.1]; % 100-1000ms
    elseif timeStep == 0.2
        limX = [0.1 20]; limY = [4e-3 0.06]; % 20-200ms
        if strcmp( strain, 'SK311')
            limX = [0.1 30]; limY = [1e-3 0.1]; % 20-200ms, with SK311
        end
    elseif timeStep - 0.02 < 2e-3
        limX = [0.01 3]; limY = [1e-3 0.05]; % 20-200ms, with SK311
    else
        limX = 'auto'; limY = 'auto';
    end
end