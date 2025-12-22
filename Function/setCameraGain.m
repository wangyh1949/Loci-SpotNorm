%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 12/21/2025
    Last update date: 12/21/2025

~~~~~~~~~~ adapted from lociAnalysis.m ~~~~~~~~~~~

Description: this script calculate conversion factors (from readout to
photon) and set pixelSize based on camera setting

    ~~~~ photon # = signal (AD)* sensitivity/ EM Gain ~~~~

Default imaging setting for loci tracking experiments
    ~~~~~~~~ uB, Andor, Gain 3, 300x ~~~~~~~~
%}


% cameraFlag = input( '\n ~~~ What is the camera (1:andor, 2:hama):  ');
cameraFlag = 1; 

if cameraFlag == 1
    pixelSize = 160* 1e-9; % Andor camera
    % gainFlag = input( '\n ~~~ What is the gain setting (1 or 3, assumed 300x):  ');
    gainFlag = 3;

    % sigToPhoton = (2^16-1)* 4.88/  10; % uB, Andor, Gain 3, 10x (SK311, 250710)
    % sigToPhoton = (2^16-1)* 15.3/ 100; % uB, Andor, Gain 1, 100x (SK311, 250628)
    % sigToPhoton = (2^16-1)* 5.11/ 100; % uA, Andor, Gain 3, 100x 

    if gainFlag == 1
        sigToPhoton = (2^16-1)* 15.3/ 300; % uB, Andor, Gain 1, 300x   (~3e3)
    elseif gainFlag == 3    
        sigToPhoton = (2^16-1)* 4.88/ 300; % uB, Andor, Gain 3, 300x   (~1e3, default)
    elseif gainFlag == 10   
        sigToPhoton = (2^16-1)* 4.88/  10; % uB, Andor, Gain 3, 10x (SK311, 250710)
    else
        sigToPhoton = input( '\n Please input the value manually as (2^16-1)* 4.88/ 300:  ');   
    end
elseif cameraFlag == 2
    pixelSize = 64.5*10^-9; % Hamamatsu camera
    fprintf( '\n ~~~ Hama, assumed gain 100x ~~~\n')
    sigToPhoton = (2^16-1)* 0.275/ (10^( 100/255));  % Hama, Gain 100  (~7e3)
else        
    error( '\n ~~~ Problematic Camera Info ~~~\n')
end