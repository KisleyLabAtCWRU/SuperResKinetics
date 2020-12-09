%% Group events to get a list of the locations of interest
% This function first makes a list of the coordinates of all identified
% events and then groups them based on a distance threshold specified in
% e.sameLocation. This function is for use with diffraction limited data. 
% For use with the super-resolution suite produced by Bo, Jixin, Joey, and
% Lydia -October 2013
%----------------------------------Inputs----------------------------------
% LocatStore - The structure containing the identified event information.
% Specifically the field PSFfinal is used by this function
% e - The structure containing the user defined parameters. Specifically
% this function uses e.sameLocation
%----------------------------------Outputs----------------------------------
% GroupLocat(1:n) - Is a structure where each element contains the fields RawSites and Centroid. 
% n is the number of identified unique sites. GroupLocat(1:n).RawSites 
% has all of the events grouped to that site and is a m by 4 matrix where m
% is the number of events in the group and each row contains the event info formated as [y x frame# I].
% GroupLocat(1:n).Centroid contains the final locations calculated by
% averaging the x and y positions of all events grouped to a unique
% location. It is a 1 by 4 matirx and has the format [y x stdy stdx]
%                                                                                           ---Joey Tauzin     Oct. 2013
function [GroupLocat] = LRG_SuperRes_FindFinalLocation(LocatStore,e)
tic % Initialize timer
sameLocation=e.sameLocation; % Create a copy of e.sameLocation in the function workspace
for f = 1:size(LocatStore,2)   % add the frame number to the list event parameters for use later in this function
     if isempty(LocatStore(1,f).PSFfinal)==0
         LocatStore(1,f).PSFfinal(:,7)=f;
     end
end

LocatStoreAll=[]; % put all molecule locations in a single matrix (LocatStoreAll) for ease of use
for f=1:size(LocatStore,2) % loop over all frames
   if isempty(LocatStore(1,f).PSFfinal)==0 
       LocatStoreAll=[LocatStoreAll; LocatStore(1,f).PSFfinal];
   end
end
 
GroupLocat=struct('Centroid',[],'RawSites',[]); % initialize storage
UniqueLocations=[];
GroupLocat(1,1).RawSites=LocatStoreAll(1,[1 2 7 5]); % start with first location
UniqueLocations(1,1:2)=LocatStoreAll(1,1:2); % Initialize the list of unique locations starting with the first location

%% Identify the unique sites
for i=2:size(LocatStoreAll,1); % loop over all events
    % now we will calculate the distance from event i to all of the unique
    % event sites identified so far
    LESqure=(LocatStoreAll(i,1)-UniqueLocations(:,1)).^2 +(LocatStoreAll(i,2)-UniqueLocations(:,2)).^2;
    % Find the minimum of the list of distances
    [mindist b]=min(LESqure);
    if mindist <= sameLocation^2 % Is the minimum less than or equal to our samelocation threshhold?
        % if so group it with other events at that site
        GroupLocat(1,b).RawSites=[GroupLocat(1,b).RawSites; LocatStoreAll(i,[1 2 7 5])];
    else
        % if not then create a new unique site
        GroupLocat(1,size(GroupLocat,2)+1).RawSites=LocatStoreAll(i,[1 2 7 5]);
        UniqueLocations(size(UniqueLocations,1)+1,1:2)=LocatStoreAll(i,1:2);
    end
end
%% Calculate final locations
% Now we will calculate the average final site location based on the
% grouped event locations
for i=1:size(GroupLocat,2) % loop over all sites
    if size(GroupLocat(1,i).RawSites,1)==1 
        GroupLocat(1,i).Centroid=[GroupLocat(1,i).RawSites(1,1:2) 0 0]; % If there is only 1 event in the group 
        % then there is no need to averge so we just add it to the final
        % location list
    else
        % calculate the location statistics
        meany=mean(GroupLocat(1,i).RawSites(:,1));
        stdy=std(GroupLocat(1,i).RawSites(:,1));
        meanx=mean(GroupLocat(1,i).RawSites(:,2));
        stdx=std(GroupLocat(1,i).RawSites(:,2));
        % and add them to the list of final sites
       GroupLocat(1,i).Centroid=[ meany, meanx stdy stdx];
    end
end           
time=toc; % Store the time it took this function to run. Useful for diagnostics. This can obviously be commented out
return
