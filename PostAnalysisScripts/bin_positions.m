
%% group sites for stdev analysis

% other parameters
e.xmin=1;
e.xmax=512;
e.ymin=1;
e.ymax=512;
e.pixelSize=64; % nm
e.nframes=size(LocatStore,2);

% SuperResolution Parameters
e.nzoom=6; % Super resolution zoom value
e.sigmarad=25; %nm sigma radius used to generate superresolution data. 

% SuperResolution Event Grouping Parameters
e.FinalLocatSigma=0.5;
e.nevent=1;
e.SRCorrFactor=0.6;
e.FinalLocatThresh=6.25;

% run analysis
addpath('U:\Lydia\MATLAB\SuperPosition-20130215\LRG_SuperRes_Kinetics_Final\')
superdata=LRG_SuperRes_GenerateSR(LocatStore,e);
disp('Identifying Final Sites of Interest')
[GroupLocat, BSALocatStore]=LRG_SuperRes_GroupsfromSR(LocatStore,superdata,e);

stdpts=[];
for i=1:size(GroupLocat,2)
    x=mean(GroupLocat(1,i).Centroid(1,3:4));
    if x>0.001
        stdpts=[stdpts, x];
    end
end

figure
hist(stdpts,sqrt(numel(stdpts)))
set(gca,'xtick',[0.16,0.31,0.47,0.625,0.78,0.94,1.09,1.25],'xticklabel',[10,20,30,40,50,60,70,80],'FontSize',14)
xlabel('Stdev (nm)','FontSize',14)
ylabel('Occurances','FontSize',14)

%% bin data representation figure
bindata= zeros((e.ymax-e.ymin+1)*e.nzoom,(e.xmax-e.xmin+1)*e.nzoom);

alllocat=[];
for i = 1:size(BSALocatStore,2)
        alllocat = [alllocat;BSALocatStore(1,i).PSFfinal(:,1:6)];   
end

for i = 1:size(alllocat,1)
    newy0 = round(alllocat(i,1)*e.nzoom);  % new center y0
    newx0 = round(alllocat(i,2)*e.nzoom);  % new center x0 

    bindata(newy0,newx0) = bindata(newy0,newx0)+1;
end

%% figure of binned data, Bo's representation
figure
subplot(1,3,1)
imagesc(bindata)
axis square
set(gca,'xtick',[],'ytick',[])
title('bin data')

subplot(1,3,2)
imagesc(cnx, cnx, highres)
axis square
set(gca,'xtick',[],'ytick',[])
title('weighted-SNR sigma of PSF')

subplot(1,3,3)
imagesc(cnx, cnx, highres2)
axis square
set(gca,'xtick',[],'ytick',[])
title('weighted-Inten. amp. of PSF')

figure
imagesc(cnx, cnx, highres)
axis square
set(gca,'xtick',[],'ytick',[])
title('weighted-SNR sigma of PSF')
caxis([0 200])

