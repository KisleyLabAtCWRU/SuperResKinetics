h=figure;

movies=[];

for i=1:size(Data1,3) %plot frames and id'd events on frames
    imagesc(Data1(:,:,i))
    colormap(gray)
    axis image
    caxis([0 1000])
    hold on
    
    if e.chngpt==1
        for j=1:size(GroupLocat,2)
            if any(i==GroupLocat(1,j).frameson)
                plot(GroupLocat(1,j).Centroid(1,2),GroupLocat(1,j).Centroid(1,1),...
                    'co','MarkerSize',10,'LineWidth',2)
            end
        end
    else
        %add stuff w/ GroupLocat(1,j).RawSites(:,3)
    end      
        for j=1:size(GroupLocat,2)
            if any(i==GroupLocat(1,j).RawSites(:,3))
                plot(GroupLocat(1,j).Centroid(1,2),GroupLocat(1,j).Centroid(1,1),...
                    'co','MarkerSize',10,'LineWidth',2)
            end
        end
        
    hold off
    
    
    
    movies=[movies getframe(h)];
    clf(h)
end

[h, w, p] = size(movies(1).cdata);
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf,'Position', [150 150 w h]);
axis off
% Place frames at bottom left
movie(hf,movies,1,25,[0 0 0 0]);

% %Save movie
% movie_name='U:\Lydia\Data\030714\singlesitemovie'

% movie2avi(movies, movie_name, 'fps',10);
