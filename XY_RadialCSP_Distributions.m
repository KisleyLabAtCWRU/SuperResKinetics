%% Radial Distributions at set Z, in XY
%main distribution analysis for CSP single molecule -set slice -data
%creates Analyzed output structure with all data, and grouped radial
%distributions
clear;close all;clc

%%
startloc='Z:\RicardoMongeN\IX-73 collections\2023_05_08\';
numFiles=3;
multiSelect = 'false';

cd(startloc);

fileNames = "";
filePaths = "";

if strcmp(multiSelect,'true')==0
    for i=1:numFiles
        [file, path] = uigetfile('.mat')
        fileNames = [fileNames, file];
        filePaths = [filePaths, path];
    end
else
   [file,path] = uigetfile(startloc,...
   'Select One or More Files', ...
   'MultiSelect', 'on');
    numFiles = size(file,2);
    for i=1:numFiles
        fileNames = [fileNames,file(i)];
        filePaths = [filePaths, path];
    end
end


disp('All Files Selected')

%% Find Radial Distribution
close all
Analyzed = struct([]);

opts.Interpreter = 'tex';
% Include the desired Default answer
opts.Default = 'Yes';
% Use the TeX interpreter to format the question
quest = 'Use Contour Fitting?';
methodFit = questdlg(quest,'Fitting Method',...
                  'Yes','No (Manual Fit)',opts)


for i=2:numel(fileNames); 
   
    
    e.filename=fileNames(i);
    
    e.codedir= filePaths(i);
    e.path= filePaths(i);
    addpath(e.path); 
    
    disp('Loading All Analyzed File')
    %load(strcat(e.filename,'_analyzed.mat'));
    %load(e.filename);
    load(strcat(e.path,e.filename));
    i = i-1;
    Analyzed(i).RadDist = [];
    Analyzed(i).SiteRadDist = [];
    Analyzed(i).numSites = [];
    Analyzed(i).ParticlesPerFrame = [];
    Analyzed(i).Dwell = [];
    Analyzed(i).xDwell = [];
    Analyzed(i).yDwell = [];

    Analyzed(i).Asc = [];
    Analyzed(i).xAsc = [];
    Analyzed(i).yAsc = [];

    Analyzed(i).Rcenter = [];
    %%%%%%%%% ONLY for contour fitted data. Empty otherwise (or not well
    %%%%%%%%% fitted)
    Analyzed(i).Elli = []; %recorded ellipse parameters
    Analyzed(i).acdArea_tot = []; %total area sampled by analyte in particle
    Analyzed(i).acdArea_In = []; %
    Analyzed(i).acdArea_Out = []; %
    Analyzed(i).acdArea_ratio = [];
    Analyzed(i).rIn_Out = [];
    %%%%%%%%%%%
    close all
            
            if strcmp(methodFit,'Yes')~=0; %Countour Fitting
                addpath('Z:\RicardoMongeN\MATLAB')
                addpath('Z:\RicardoMongeN\MATLAB\LRG_SuperRes_Kinetics_Final\LRG_SuperRes_Kinetics_Final')
                %e.nzoom = 20;
                if i<numFiles
                        e.sigmarad = 25; %nm
                    else
                        e.sigmarad = 200;
                end
                D2 = LRG_SuperRes_GenerateSR(LocatStore,e);
                %D4 = sum(sum(Data1,3),3);
                figure();imagesc(D2);
                D3 = imbinarize(D2); figure();imagesc(D3);
                                
                pctList = [0.75 0];
                M1 = max(max(D3)); %contour based on multiples of max intensity
                mLims =[];
                for j =1:size(pctList,2)
                    mLims = [mLims, pctList(j)*M1]; %number of levels at which contours are drawn (defualt two at 75% and 50%)
                end
                %figure();hold on;
                [I3,h3]  = imcontour(D3,mLims);

                if isempty(I3)==0 
                %close;close
                % Extract Contour Data
                % clc

                S = contourdata(I3);
                numLvls = numel(mLims);
                SnumEl = [];
                for j=1:size(S,2)
                    SnumEl = [SnumEl, S(j).numel];   
                end
                [mX,mY] = maxk(SnumEl,numLvls);

                % Fit the Contour Data to Ellipses
                clc 
                AS = [];%Area
                RS = [];%R^2 value of fit

                aX = figure();hold on;
                for j=1:numLvls
                    x = S(mY(j)).xdata;
                    y = S(mY(j)).ydata;
                    plot(x,y);
                %     plot(x,y);axis([0,512,0,512]);
                    Elli(j) = fit_ellipse(x,y,aX); %fitted ellipse parameters
                    if isempty(Elli(j).a)==0
                    Analyzed(i).Elli(j,:) = [Elli(j).a, Elli(j).b,...
                        Elli(j).X0_in, Elli(j).Y0_in].*((e.pixelSize/e.nzoom)/1000) ; %in um
                        %[a , b, x0, y0];
                    AS = [AS,pi*Elli(j).a*Elli(j).b*((e.pixelSize/(e.nzoom*1000))^2)];%area of ellipse, um^2
                    end
                % %    
                end
                end
                
                %Confirm Fit Quality
                opts.Interpreter = 'tex';
                % Include the desired Default answer
                opts.Default = 'Yes';
                % Use the TeX interpreter to format the question
                quest = 'Is this a good fit?';
                goodFit = questdlg(quest,'Contour Fitting',...
                                  'Yes','No (Cancel)',opts)
                if strcmp(goodFit,'Yes')~=0; %save the fitted ellipse x, y, and radius
                    fitXlim = ylim(gca);
                    fitYlim = xlim(gca);               
                   
                    mEll = [];
                    for j=1:numLvls;
                        
                        mEll= [mEll, Elli(j).a];
                    end
                    [maEll, iEll] = max(mEll); %inner circle
                    [mIll, iMll] = min(mEll); %outer circle
                    
                    Rad0 = mean([Elli(iEll).a,Elli(iEll).b])*(e.pixelSize/e.nzoom)/1000;%in um, radius
                    Rcenter = ([Elli(iEll).X0_in, Elli(iEll).Y0_in].*e.pixelSize/e.nzoom)/1000; %in um
                    Rcenter = round([Rcenter,Rad0],2); %center of circular ROI, [x0, y0, R0] in um
                    
                    Analyzed(i).acdArea_Out = pi().*(mean([Elli(iEll).a,Elli(iEll).b])*(e.pixelSize/e.nzoom)/1000)^2 ;%in um^2
                    Analyzed(i).acdArea_In = pi().*(mean([Elli(iMll).a,Elli(iMll).b])*(e.pixelSize/e.nzoom)/1000)^2 ;%in um^2
                    Analyzed(i).rIn_Out = abs((mean([Elli(iMll).a,Elli(iMll).b])...
                        -mean([Elli(iEll).a,Elli(iEll).b]))*(e.pixelSize/e.nzoom)/1000); %porous layer thickness
                    
                    close all
                    
                    figure();hold on;
                    %Binned = (sum(Data1,3));
                    %imagesc(Binned);
                        e.sigmarad = 25; %nm

                    D4 = LRG_SuperRes_GenerateSR(LocatStore,e);D4 (D4 < 0.5) = 0;
                    D4 = imbinarize(D4);imagesc(D4);
                    %imagesc(D2);
                                      
                   
                    axis image;colorbar;%colormap hot;
                     % the ellipse
                        theta_r         = linspace(0,2*pi);
                        ellipse_x_r     = [Elli(iEll).X0_in + Elli(iEll).a*cos( theta_r )];%./e.nzoom;
                        ellipse_y_r     = [Elli(iEll).Y0_in + Elli(iEll).b*sin( theta_r )];%./e.nzoom;
                        rotated_ellipse = 1* [ellipse_x_r;ellipse_y_r];

                        % draw
                        hold all; plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
                    
                    
                    D4 = D4(min(rotated_ellipse(2,:)):max(rotated_ellipse(2,:)),...
                        min(rotated_ellipse(1,:)):max(rotated_ellipse(1,:)));
                    figure();imagesc(D4);
                    
                    SpA2 = bwarea(D4)*((e.pixelSize/e.nzoom)/1000)^2;%in um^2
                    Analyzed(i).acdArea_tot = SpA2; 
                    Analyzed(i).acdArea_ratio = SpA2/Analyzed(i).acdArea_Out;
                    
                    disp('Press any button or click figure to proceed')
                    waitforbuttonpress
                else % Cancel and start over
                    if strcmp(goodFit,'No (Cancel)')~=0;
                       error('Canceled Operation')
                       close all
                       return
                    end
                end
                
            else %Manual Fitting
                
                close all
                 figure(1);hold on;
                    Binned = (sum(Data1,3));
                    imagesc(Binned);
                    axis image;colorbar;%colormap hot;
                    
                if i<3
                    ROI = imellipse(gca,[50 50 50 50]); %plot initial rectangle for choosing ROI
                else
                    ROI = imellipse(gca,tempROI);
                end
                hold on
                title({'Shape Circle and double-click to choose region of interest' ...
                    mat2str(ROI.getPosition)});
                addNewPositionCallback(ROI,@(p) title({...
                    'Double-click to choose region of interest,use Fix Aspect Ratio option' ...
                    ['[xmin ymin width height] = ' mat2str(p,3)]}));
                % make sure rectangle stays within image bounds
                fcn = makeConstrainToRectFcn('imellipse',get(gca,'XLim'),get(gca,'YLim'));
                setPositionConstraintFcn(ROI,fcn);
                tempROI = wait(ROI);
                tempROI = getPosition(ROI);
                close;
                Rcenter = round([tempROI(1)+(tempROI(3)/2) , tempROI(2)+(tempROI(4)/2)],3);
                Rad0 = (mean([tempROI(3),tempROI(4)])*e.pixelSize/1000)/2;%in um, radius
                Rcenter = (Rcenter.*e.pixelSize)/1000; %in um
                Rcenter = round([Rcenter,Rad0],2); %center of circular ROI, [x0, y0, R0] in um
                
            end
  
    
    %particle distribution
    kCount = 0;%particle counter, over all frames
    for j=1:size(LocatStore,2)
        
        for k=1:size(LocatStore(j).PSFfinal,1)
            
            delX = (LocatStore(j).PSFfinal(k,2).*e.pixelSize/1000) - Rcenter(1);%um
            delY = (LocatStore(j).PSFfinal(k,1).*e.pixelSize/1000) - Rcenter(2);   %um
            
            Rdist = sqrt(delX^2 + delY^2); %um
            if Rdist <= round(Rad0,2)+0.2;
                kCount = kCount +1;
                Analyzed(i).RadDist(kCount) = ((round(Rad0,2)-Rdist)/round(Rad0,2))*100; %um  -> percent of Radius
                
            end
            
        end

        
    end
    Analyzed(i).ParticlesPerFrame = (kCount)/size(LocatStore,2);
    Analyzed(i).particleDensity = ((kCount)/size(LocatStore,2))/Analyzed(i).acdArea_Out;
    
    %site distribution
    for j=1:size(GroupLocat,2)
            
            delX = (GroupLocat(j).Centroid(2).*e.pixelSize/1000) - Rcenter(1);%um
            delY = (GroupLocat(j).Centroid(1).*e.pixelSize/1000) - Rcenter(2);   %um
            
            Analyzed(i).SiteRadDist(j) = ((round(Rad0,2)-sqrt(delX^2 + delY^2))/round(Rad0,2))*100; %um -> percent of Radius
            
           
    end
    Analyzed(i).Dwell = Ensemble.Dwell(:);  
    [Analyzed(i).xDwell, Analyzed(i).yDwell]= cumuldist(transpose(Ensemble.Dwell(:)),...
                                                transpose(unique(Ensemble.Dwell(:))));

     Analyzed(i).Asc = Ensemble.Assoc(:);  
    [Analyzed(i).xAsc, Analyzed(i).yAsc]= cumuldist(transpose(Ensemble.Assoc(:)),...
                                                transpose(unique(Ensemble.Assoc(:))));
    
    Analyzed(i).numSites = size(GroupLocat,2);
    Analyzed(i).Rcenter = Rcenter;
    
    i=i+1;
end
close all
% save(strcat(filePaths(i),"Analyzed.mat"),'Analyzed');
disp('Done')
%% Export Variables
close all
pRD = [];% particle radial distribution
sRD = [];% sites radial distribution

fullDwell = [];
fullAsc = [];


for i=2:size(Analyzed,2)+1; 
    Analyzed(i-1).xDist = [];
    Analyzed(i-1).yDist = [];
    figure(1);hold on; histogram(Analyzed(i-1).RadDist(:));
    xlabel('Radial Distance from edge (%)');ylabel('Identified Particles');
    
    %Particle Dist
    [xd,yd]=cumuldist(transpose(Analyzed(i-1).RadDist(:)),transpose(unique(Analyzed(i-1).RadDist(:))));
    Analyzed(i-1).xDist = xd;
    Analyzed(i-1).yDist = yd;    
figure(2);hold on
    plot(xd,yd,'o-')
    set(gca,'YScale','linear');
    xlabel('Distance from edge(%)')
    ylabel("P(r'>r)")
    title('Radial Distance Distribution') 
    pRD = [pRD, transpose(Analyzed(i-1).RadDist(:))];
    
    %Sites Dist
    [xd,yd]=cumuldist(transpose(Analyzed(i-1).SiteRadDist(:)),transpose(unique(Analyzed(i-1).SiteRadDist(:))));
figure(3);hold on
    plot(xd,yd,'o-')
    set(gca,'YScale','log');
    xlabel('Distance from edge(%)')
    ylabel("P(r'>r)")
    title('Sites Radial Distance Distribution')
    sRD = [sRD, transpose(Analyzed(i-1).SiteRadDist(:))];
    
    %Dwell Times Dist
    fullDwell = [fullDwell, transpose(Analyzed(i-1).Dwell)];
    fullAsc = [fullAsc, transpose(Analyzed(i-1).Asc)];
end



[pXD,pYD]=cumuldist(pRD,unique(pRD));
figure(2);hold on;legend();
    plot(pXD,pYD,'o-')
     %set(gca,'YScale','log');
%     xlabel('Distance (um)')
     ylabel("P(r'>r)")
%     title('Radial Distance Distribution')

figure(4);histogram(pRD);xlabel('Radial Distance from edge (%)');ylabel('Particles');
figure(5);histogram(sRD);xlabel('Sites Radial Distance from edge (%)');ylabel('Sites');

[sDX,sDY]=cumuldist(sRD,unique(sRD));
figure(3);hold on;
    plot(sDX,sDY,'o-')
%     set(gca,'YScale','log');
%     xlabel('Distance (um)')
%     ylabel("P(r'>r)")
%     title('Radial Distance Distribution')

[xDwell, yDwell] = cumuldist(fullDwell, unique(fullDwell));
figure(6);
for i=2:size(Analyzed,2)+1; 
    figure(6);hold on;
    plot(Analyzed(i-1).xDwell,Analyzed(i-1).yDwell,'o-')
end
    plot(xDwell,yDwell,'o-')
    set(gca,'YScale','log');
    xlabel('Dwell (ms)')
    ylabel("P(t'>t)")
    title('Overall Dwell Time Distribution');legend();
    

%% Save Distribution Data

save(strcat(filePaths(i),"3RadiiDistributionAnalyzed3.mat"),'Analyzed',...
    'fullDwell','fullAsc','xDwell','yDwell',...
    'pRD','pXD','pYD',...
    'sRD','sDX','sDY');
saveas(figure(2),strcat(filePaths(i),"Dist.fig"));
disp('Done & Saved Distribution')

%%  Uncertainty Error calculation
% xLength = [];
% for i=2:numel(fileNames);
%     xLength = [xLength, size(Analyzed(i-1).yDist,1)];
% end
% [mA, mAid] = max(xLength);
% simErrX = zeros(1,mA);%[];
% simErrY = zeros(1,mA);%[];
% for i=2:numel(fileNames); 
%     %Analyzed(i-1).xDist
%     roundPop = round(pXD,3);%round to 3 dec - ensemble
%     roundA_i = round(Analyzed(i-1).xDist,3);% round current data set
%     roundA_max= round(Analyzed(mAid).xDist,3);%round data for largest data set
% 
%     [iCA, iCB] = ismember(roundA_i,roundPop);%get indeces where X values overlap
%     [iCA, iCBi] = ismember(roundPop,roundA_i);%get indeces where X values overlap
%     roundA_i = roundA_i(nonzeros(iCBi));
%     [iAA,iAB] = ismember(roundA_i,roundA_max); %where data sets overlap
%     iAB = nonzeros(iAB);
%     for j =1:length(iAB);
%         simErrX(iAB(j)) = Analyzed(i-1).xDist(iAB(j));%make list of indexes X 
%         simErrY(iAB(j)) =  simErrY(iAB(j))+transpose((pYD(iCB)-Analyzed(i-1).yDist));
%     %simErrY = [simErrY, transpose((pYD(iCB)-Analyzed(i-1).yDist))];
%     end
%     
% end
% [simErrX2, sortIDx] = sort(simErrX,'ascend');%sort X values
% simDatY2 = pYD(unique(simErrX2));
% simErrX2 = pXD(unique(simErrX2));
% simErrY2 = simErrY(sortIDx); %sort Y values to match X sorting
% 
% 
% globalErr = (simErrY2./(numel(fileNames)-1))./2; %get the associated error as the mean difference
% 
% widSize = 3;
% globalErr2 = zeros(1,size(globalErr,2));
% for i=1:(size(globalErr,2)-widSize);
%     globalErr2(1,i)= mean(globalErr(1,i:i+widSize));
% end
% 
% figure(7);hold on; h = errorbar(simErrX2,simDatY2,globalErr2,...
%     "--s","MarkerSize",8,'LineWidth',1);
% legend();
% % Set transparency level (0:1)
% alpha = 0.2;   
% % Set transparency (undocumented)
% set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha]);
