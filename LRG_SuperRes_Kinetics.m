function [GroupLocat, Ensemble]=LRG_SuperRes_Kinetics(GroupLocat,e)

%% Kinetics_LK102113
% Kinetics code for application in super-resolution suite produced by Bo,
% Jixin, Joey, and Lydia October 2013

% INPUT: GroupLocat-locations determined by Jixin's super-resolution or Joey's FindFinalLocat 
%        e - global variable; use e.dataSpace, e.dataTime, e.datafreq
% OUTPUT:   GroupLocat.Dwell - single location dwell times
%           GroupLocat.Assoc - single location association times
%           Ensemble.Dwell - ensemble dwell times
%           Ensemble.Assoc - ensemble association times
% DEFINITIONS: 
% dwell time - the time a molecule is present for consecutive frames
% assoc time - the time at an individual site from when one molecule leaves
%              till the next molecule arrives

% Calculate single site kinetics
for i=1:size(GroupLocat,2) %cycle through all individual adsorption sites
    TimeSort=sort(GroupLocat(1,i).RawSites(:,3));  %sort in events in order of which frame they appeared
    dwell=[]; %initialize storage for all dwell/assoc times at the site
    assoc=[];
    dwellflag=1; %intitialize flags for counting between events
    for j=2:numel(TimeSort)
        if TimeSort(j-1)==TimeSort(j)-1 %calculate dwell times
            dwellflag=dwellflag+1;   
        else
            if TimeSort(j-dwellflag)~=1 %don't count events that started in first frame
                dwell=[dwell, dwellflag];
            end
            dwellflag=1;
        end
        
        if TimeSort(j-1)~=TimeSort(j)-1 %calculate assoc times
            assoc=[assoc, (TimeSort(j)-TimeSort(j-1))];
        end
    end
    
    % convert dwell and assoc times to ms instead of frame nos.
    GroupLocat(1,i).Dwell=dwell.*((e.dataSpace/2)+e.dataTime); %convert dwell time to ms
    % Based on Jixin's explination for "average time" since we don't know what
    % is occuring during the dark time (explination based on 30 ms int/32 ms down):
    % "The first frame of the dwell time is set to 46.5 ms which is the average 
    % of 30+32 ms dark+light+dark time, meaning that when an event shows up 
    % only once, the average dwell time should be 46.5 ms instead of 62 ms."
    GroupLocat(1,i).Assoc=assoc.*e.datafreq; %frames * frame rate in ms
end
    
% Calculate ensemble kinetics
Ensemble.Dwell=[]; %initialize storage variable
Ensemble.Assoc=[];

for i=1:size(GroupLocat,2) %cycle through all individual adsorption sites
    Ensemble.Dwell=[Ensemble.Dwell, GroupLocat(1,i).Dwell]; %store all dwell times
    Ensemble.Assoc=[Ensemble.Assoc, GroupLocat(1,i).Assoc]; %store all association times
end

end