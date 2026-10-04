function tInt = findInt( extraName, expDate, strain)
    % this function find the tInt text, containing exposure & frame
    % interval information, for example
    % 20-200ms --> 200 (ms)
    %    200ms --> 200 (ms)

    % Created at 6/11/2024
    
    tInt = regexp( extraName, '\d+-\d+ms', 'match', 'once');
    if isempty( tInt)
        tInt = regexp( extraName, '\d+ms', 'match', 'once'); 
    end