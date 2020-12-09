function [GroupLocat, Ensemble]=LRG_SuperRes_Kinetics_CP(GroupLocat,e,Data1)

%% Kinetics_LK102113
% Kinetics code for application in super-resolution suite produced by Bo,
% Jixin, Joey, and Lydia October 2013

% LK021914 - update LRG_SuperRes_Kinetics to include change point analysis 
%            to determine kinetics instead of molecule ID. This is
%            important for data with low S/N or high particle density. New
%            e variables added - e.chngpt determines if the change point
%            analysis is being used.  e.CPstate lists the number of user 
%            defined states or have the program determine it.

% INPUT: GroupLocat-locations determined by Jixin's super-resolution or Joey's FindFinalLocat 
%        e - global variable; use e.dataSpace, e.dataTime, e.datafreq
%                                 e.chngpt,e.CPstate,e.nframes
%        Data1 - the raw data 
% OUTPUT:   GroupLocat.Dwell - single location dwell times
%           GroupLocat.Assoc - single location association times
%  (for CP  GroupLocat.IntCP - Intensity (row1) and change point (row2) traces
%   only):  GroupLocat.nst_localthd - # states ID'd and local threshold
%                                     used for site
%           Ensemble.Dwell - ensemble dwell times
%           Ensemble.Assoc - ensemble association times
% DEFINITIONS: 
% dwell time - the time a molecule is present for consecutive frames
% assoc time - the time at an individual site from when one molecule leaves
%              till the next molecule arrives


if e.chngpt==0 % Calculate single site kinetics by molecule ID
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
else % Calculate single site kinetics by change point of intensity trace
    
    for i=1:size(Data1,3) %calculate average local intensity threshold map
        local_thd(:,:,i)=LRG_SuperRes_LocalThrMap(Data1(:,:,i));
    end
    thd_map=mean(local_thd,3); %calculated intensity for e.wide2 pixels

    for i=1:size(GroupLocat,2) %cycle through all individual adsorption sites
        ycenter=round(GroupLocat(1,i).Centroid(1,1));
        xcenter=round(GroupLocat(1,i).Centroid(1,2));

        if ycenter+5<size(Data1,1) && ycenter-5>0 && xcenter+5<size(Data1,2) && xcenter-5>0 %make sure not at edge of image

            for g=1:e.nframes %get intensity of 5x5 pixels around center
                k(g)=sum(sum(Data1(ycenter-5:ycenter+5,xcenter-5:xcenter+5,g)));
            end
            
            % calculate change point analysis
            addpath(strcat((e.codedir),'\change-point'));
            if e.CPstate~=0 
                % user defined no. of states
                [G, BIC, states, breaks, eff, output, records, excluded] = main_wavelet_CP_userstates(k,e.CPstate);
            else
                % analysis determines no. of states
                [G, BIC, states, breaks, eff, output, records, excluded] = main_wavelet_CP(k);
            end

            dwell1=[]; %initialize storage for all dwell/assoc times at the site
            assoc1=[];
            dwell=[];
            assoc=[];
            siteon=[];
            dwellflag=1; %intitialize flags for counting between events
            assocflag=1;
            
            [minstate indexmin]=min(states); %calculate minimum state
            thd_center=sum(sum((thd_map(ycenter-5:ycenter+5,xcenter-5:xcenter+5)))); %calculate local threshold

%             if minstate <= thd_center% make sure 'min' is a true 'off state' based on comparison to local threshold
%                 minindexmin=min(indexmin);
%                 for j=minindexmin:numel(states)-1 %store change points as dwell/assoc times                    
%                     if states(j)~=minstate %calc dwell times only as changes from min. (background) state
%                         dwellflag=dwellflag+1;
%                         assoc1=[assoc1,assocflag];
%                         assocflag=1;
%                         siteon1=[siteon, j]; %add storage of the 'on' frames for figure making
%                     else
%                         dwell1=[dwell1,dwellflag];
%                         dwellflag=1;
%                         assocflag=assocflag+1;
%                     end
%                 end
% 
%                 %calc dwell times (dwell1=1 not real/fillers)
%                 for j=1:numel(dwell1)
%                     if dwell1(j)~=1
%                         dwell=[dwell, dwell1(j)];
%                     end
%                 end
% 
% 
%                 % calc assoc times (assoc1=1 not real/fillers)
%                 for j=1:numel(assoc1)
%                     if assoc1(j)~=1
%                         assoc=[assoc, assoc1(j)];
%                     end
%                 end
% 
%                 % trim assoc if first frame is at the min state
%                 
%                     if minindexmin==1
%                         if numel(assoc)>1
%                             assoc=assoc(2):assoc(numel(assoc));
%                         else
%                             assoc=[];
%                         end
%                     end
% 
% %                 figure
% %                 hold on
% %                 plot(eff);
% %                 plot(states, '-r');
% %                 title(strcat('molec at (x,y) (', num2str(xcenter),',',num2str(ycenter),') no states=',num2str(size(unique(states)))));
%             end
            if minstate <= thd_center
                minindexmin=min(indexmin);
                for j=minindexmin:numel(states)-1 %store change points as dwell/assoc times                    
                    if states(j)~=minstate
                        dwellflag=dwellflag+1;%increase dwell time by 1 frame
                        if assocflag~=1
                            assoc=[assoc assocflag]; %record assoc time if first change
                        end
                        assocflag=1;
                        siteon=[siteon j]; %record frame # of 'on' event
                    else
                        assocflag=assocflag+1; %increase assoc time by 1 frame
                        if dwellflag~=1
                            dwell=[dwell dwellflag]; %record dwell time if change
                        end
                        dwellflag=1;
                    end
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
            
            GroupLocat(1,i).IntCP=[eff;states];%store intensity and cp states
            GroupLocat(1,i).nst_locthd=[size(unique(states),2),thd_center];
            GroupLocat(1,i).frameson=siteon;%store the frame #s where a molecule is ID'd

        end
    end
end
    
% Calculate ensemble kinetics
Ensemble.Dwell=[]; %initialize storage variable
Ensemble.Assoc=[];

for i=1:size(GroupLocat,2) %cycle through all individual adsorption sites
    Ensemble.Dwell=[Ensemble.Dwell, GroupLocat(1,i).Dwell]; %store all dwell times
    Ensemble.Assoc=[Ensemble.Assoc, GroupLocat(1,i).Assoc]; %store all association times
end

end