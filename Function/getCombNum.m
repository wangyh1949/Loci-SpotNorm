%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 12/8/2025
    Last updated at 12/8/2025

Description: this script ask you to choose files from a list for combining

input: 
    0 - combine all available files
    1 3 4 - combine file #1, 3, 4

---------------------------------------------------------------------------
%}

function combNum = getCombNum( tfList)

    listN = length( tfList); 
    
    if listN > 1
        fprintf( '\n~~~~~ There are %d tracking files ~~~~~\n', length( tfList))
        
        % list all tracking folders in the directory and ask which to combine
        for k = 1: listN
            fprintf( '%3d. %s\n', k, tfList( k).name)
        end

        tmp = string( input( '\nWhich files do you want to combine (0: all, or like: 1 3 4): ', 's'));
        combNum = double( regexp( tmp, '\d+', 'match')); % find all numbers but not extra space
        
        if combNum == 0, combNum = 1: listN; end
    end

    fprintf( '\n')
end