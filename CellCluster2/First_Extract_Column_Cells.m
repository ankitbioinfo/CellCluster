

clear all

% 
% 
% dt_path_wt={    '../data/Nuclei_and_Cells_DT_S84_m3_wt/','../data/Nuclei_and_Cells_DT_S17_m2_wt/','../data/Nuclei_and_Cells_DT_S84_m4_wt/',...
% 		'../data/Nuclei_and_Cells_DT_S18_m6_wt/','../data/Nuclei_and_Cells_DT_S51_m2_wt/'};
% pt_path_wt = {  '../data/Nuclei_and_Cells_PT_S84_m3_wt/','../data/Nuclei_and_Cells_PT_S17_m2_wt/','../data/Nuclei_and_Cells_PT_S84_m4_wt/',...  
%                 '../data/Nuclei_and_Cells_PT_S18_m6_wt/','../data/Nuclei_and_Cells_PT_S51_m2_wt/'};
% dt_path_mut= {'../data/Nuclei_and_Cells_DT_S17_m1_mut/', '../data/Nuclei_and_Cells_DT_S18_m2_mut/' ,...
%               '../data/Nuclei_and_Cells_DT_S84_m1_mut/', '../data/Nuclei_and_Cells_DT_S84_m5_mut/'};
% pt_path_mut = {'../data/Nuclei_and_Cells_PT_S17_m1_mut/', '../data/Nuclei_and_Cells_PT_S18_m2_mut/',...
%                '../data/Nuclei_and_Cells_PT_S84_m1_mut/', '../data/Nuclei_and_Cells_PT_S84_m5_mut/', };
% du_path_wt={'../data/Nuclei_and_Cells_DU_S51_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m3_wt/'};


dt_path_wt={ '../data/Nuclei_and_Cells_DT_S18_m6_wt/', '../data/Nuclei_and_Cells_DT_S17_m2_wt/',...
             '../data/Nuclei_and_Cells_DT_S84_m3_wt/', '../data/Nuclei_and_Cells_DT_S51_m2_wt/',...
             '../data/Nuclei_and_Cells_DT_S84_m4_wt/'};
pt_path_wt = {  '../data/Nuclei_and_Cells_PT_S18_m6_wt/','../data/Nuclei_and_Cells_PT_S17_m2_wt/',...  
                '../data/Nuclei_and_Cells_PT_S84_m3_wt/','../data/Nuclei_and_Cells_PT_S51_m2_wt/',...
                '../data/Nuclei_and_Cells_PT_S84_m4_wt/'};
dt_path_mut= {'../data/Nuclei_and_Cells_DT_S17_m1_mut/', '../data/Nuclei_and_Cells_DT_S18_m2_mut/' ,...
              '../data/Nuclei_and_Cells_DT_S84_m1_mut/', '../data/Nuclei_and_Cells_DT_S84_m5_mut/'};    
pt_path_mut = {'../data/Nuclei_and_Cells_PT_S17_m1_mut/', '../data/Nuclei_and_Cells_PT_S18_m2_mut/',...
               '../data/Nuclei_and_Cells_PT_S84_m1_mut/', '../data/Nuclei_and_Cells_PT_S84_m5_mut/', };
du_path_wt={'../data/Nuclei_and_Cells_DU_S51_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m3_wt/'};
       
allpath={dt_path_wt; pt_path_wt; dt_path_mut; pt_path_mut; du_path_wt}; 
% 
% 
% RZ{1}={[1:10],[16:25],[1:10],[1:14],[1:7]};
% RZ{2}={[1:10],[17:34],[1:12],[20:39],[16:30]};
% RZ{3}={[1:8],[15:22],[1:8],[17:25]};
% RZ{4}={[16:29],[1:14],[11:20],[1:10]};
% RZ{5}={[15:20],[1:6],[1:6]};
% 
% PZ{1}={[11:16],[11:15],[11:16],[15:18],[8:14]};
% PZ{2}={[11:20],[12:16],[13:21],[10:19],[12:15]};
% PZ{3}={[9:14],[7:14],[9:14],[8:16]};
% PZ{4}={[8:15],[15:22],[7:10],[11:17]};
% PZ{5}={[10:14],[7:12],[7:12]};
% 
% PHZ{1}={[17:19],[8:10],[17:19],[19:22],[15:17]};
% PHZ{2}={[21:24],[8:11],[22:25],[6:9],[9:11]};
% PHZ{3}={[15:17],[4:6],[15:17],[6:7]};
% PHZ{4}={[4:7],[23:26],[4:6],[18:19]};
% PHZ{5}={[7:9],[13:14],[13:15]};
% 
% HZ{1}={[20:22],[1:7],[20:22],[23:30],[18:23]};
% HZ{2}={[25:30],[1:7],[26:31],[1:5],[1:8]};
% HZ{3}={[18:23],[1:3],[18:20],[1:5]};
% HZ{4}={[1:3],[27:29],[1:3],[20:24]};
% HZ{5}={[1:6],[15:21],[16:21]};


meanOfAllCelltemp=load('meanOfAllCell.dat');
count=1;
for gi=1:length(allpath)
	for gj=1:length(allpath{gi})
        GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,1:3);
        count=count+1;
    end
end
 

for gi=4:length(allpath)
	for gj=1:length(allpath{gi})
		path=allpath{gi}{gj};
		disp(path)
        s=strsplit(path,'Nuclei_and_Cells_');
        outputpath=strcat('MakeListColumnarStructurePrediction/',s{2});    
        if ~exist([outputpath],'dir')
              mkdir([outputpath]);
        end

            
      if exist([outputpath,'centroid_and_surface_cells.mat'], 'file') == 0  

          alignment=load(['./../../',path,'Alignment_matrix.dat']);
            vec=alignment(:,1:3);


             [numbers,txt,raw] = xlsread(['./../../',path,'Tile_coordinates.xlsx']);
               coordinates = zeros(size(txt,1)-3,5);
                for i = 4:size(txt,1),
                    temp =  char(txt(i,1));
                    res = strsplit(temp,'_POS');
                    coordinates(i-3,1) = str2num(char(res(2)));
                    coordinates(i-3,2:5) = numbers(i-3,:);
                end
                tile=coordinates(:,2:end);


    clear alltileid
    clear Repeat_surfaces
    %clear Repeat_ellipsoid
    clear Repeat_centroids
    clear Repeat_pixels
    %clear Repeat_triangulation
    clear Repeat_volume
    
scount=1;
for position = coordinates(:,1)'  % intersect(PZ{gi}{gj},coordinates(:,1)')        %setdiff(coordinates(:,1)',RZ{gi}{gj}) 
    load(strcat('./../../',path,'c_n_pos',num2str(position),' (Characteristics).mat'));
    indtemp = find(G.inter.volume_ratio>1);
    fi = find(coordinates(:,1) == position);
    for i=1:length(indtemp)
        j=indtemp(i);
    %for j=1:size(C.coords,1)
            value=G.cel.centroids(j,:)+ repmat([tile(fi,2), -tile(fi,1),tile(fi,3)],1,1);   
            Repeat_centroids(scount,:)=value*vec;

            value=C.surfaces(j).vertices + repmat([tile(fi,2), -tile(fi,1),tile(fi,3)],1,1);
            Repeat_surfaces{scount,1}=value*vec;
            %Repeat_triangulation{scount,1}=C.surfaces(j).faces;
            %value={G.cel.ellipsoid_center(j,:), G.cel.ellipsoid_radii(j,:), G.cel.ellipsoid_evecs{j}, G.cel.ellipsoid_v{j}};   
            %Repeat_ellipsoid{scount,1}=value;
            Repeat_volume(scount,:)=G.cel.volume(j,:);
            Repeat_pixels(scount,:) = calc_centroids(C.masks(j), C.origins(j,:));
            alltileid(scount,:)=[fi,j];
         
            scount=scount+1;
    end    
 
end


 ia=RemoveBadcell(Repeat_centroids);
 [size(Repeat_centroids,1), size(ia,1)]

centroid=Repeat_centroids(ia,:);
nuc=Repeat_surfaces(ia,:);
%fitellipsoid=Repeat_ellipsoid(ia,:);
unique_pixel=Repeat_pixels(ia,:);
unique_tileid=alltileid(ia,:);
%faces=Repeat_triangulation(ia,:);
celvolume= Repeat_volume(ia,:);

save([outputpath,'centroid_and_surface_cells.mat'],'nuc','centroid','celvolume','unique_pixel','unique_tileid','-v7.3');
      else
              load([outputpath,'centroid_and_surface_cells.mat']);
      end

        % meanCentroid=mean(centroid);
        % [meanCentroid, GlobalCenter{gi}{gj}];
        % [min(centroid(:,3)),max(centroid(:,3))]
        bonetype=gi;
        if (bonetype==3)|(bonetype==1)
                tempCentroid=[centroid(:,1:2),-centroid(:,3)];
        else
                tempCentroid=[centroid(:,1:2),centroid(:,3)];
        end

        newcentroid=tempCentroid-GlobalCenter{gi}{gj};
        PD_bins=linspace(min(newcentroid(:,3)),max(newcentroid(:,3)),11);
        
        if exist([outputpath,'threshold_along_PD_axis.mat'], 'file') == 0  
        
            disp('calculate minimum negihbors distance');
            [min_neigh_dist,clist,cel_normalizationFactor] = calculate_nuclei_density(newcentroid, [1, 1, 1], 2);
            for i=1:10
                index=find( (newcentroid(:,3)>PD_bins(i)) &  (newcentroid(:,3)<=PD_bins(i+1)));
                threshold(i)=mean(min_neigh_dist(index));
            end

            save([outputpath,'threshold_along_PD_axis.mat'],'threshold','-v7.3');
            dlmwrite([outputpath,'threshold_along_PD_axis.dat'],threshold,'\n');
        
        else
              load([outputpath,'threshold_along_PD_axis.mat']);
        end
        
        
        

%cellsInPZandPHZindex= find((newcentroid(:,3)>-100)&(newcentroid(:,3)<1000));
%cellsInPZandPHZindex=find(newcentroid(:,3)>PD_bins(3));
%RZtoHZrange=[min(newcentroid(:,3)),max(newcentroid(:,3))]
%CellsInPZandPHZrange=[size(centroid,1),length(cellsInPZandPHZindex)]
%disp('I am here')



d=20;
fid=fopen([outputpath,'NeighboringCell_in_20_micron_cube.dat'],'w');
for i=1:length(centroid)
    main=centroid(i,:);
    normalized_main=newcentroid(i,:);
    for j=1:10
          if( (normalized_main(3)>PD_bins(j)) &  (normalized_main(3)<=PD_bins(j+1)))
                cutoff=threshold(j);
          end
    end
    
    %ankit(i,:)=[normalized_main(:,3),cutoff];
    
    index=find(((main(1)-d)<=centroid(:,1)) & ((main(1)+d)>=centroid(:,1)) & ((main(2)-d)<=centroid(:,2)) & ((main(2)+d)>=centroid(:,2) ) & ((main(3)-d)<=centroid(:,3)) & ((main(3)+d)>=centroid(:,3) )    );
    if length(index)>1
        col=distance_between_neighbor_cell(centroid, index,cutoff);
        for j=1:size(col,1)              
              fprintf(fid,'%d\t%d\t%0.4f\t%0.4f\n',col(j,1),col(j,2),col(j,3),col(j,4));
        end         
    end
  
end
fclose(fid);







    end
end



function [column]=distance_between_neighbor_cell(Cent,ind,cutoff)
         C=Cent(ind,:);
         n=size(C,1);
         c2=1;
         column=[];
         for i=1:n
             for j=i+1:n
                 ed=pdist(C([i,j],:));  
                 zd=diff(C([i,j],3));
                 if ed<=cutoff
                     column(c2,:)=[ind(i),ind(j),ed,abs(zd)];
                     c2=c2+1;
                 end
             end
         end
        
end



function centroids = calc_centroids(masks, origins)
centroids = nan(length(masks),3);
for i = 1 : length(masks)
    [x,y,z] = ind2sub(size(masks{i}), find(masks{i}));
    centroids(i,:) = (mean([x,y,z],1) + double(origins(i,:) - 1));
end
end




function goodcellindex=RemoveBadcell(centroid)
        %goodOrBadId=zeros(size(centroid,1),1);
        n=size(centroid,1);
        d=10;
        count=1;
        
        edges=[];
        for i=1:n
            main=centroid(i,:);
            index=find(((main(1)-d)<=centroid(:,1)) & ((main(1)+d)>=centroid(:,1)) & ((main(2)-d)<=centroid(:,2)) & ((main(2)+d)>=centroid(:,2) ) & ((main(3)-d)<=centroid(:,3)) & ((main(3)+d)>=centroid(:,3) )    );
            ankit(i)=length(index);
            for j=1:length(index)
                for k=j+1:length(index)
                    dist=pdist(centroid([index(j),index(k)],:));
                    if dist<2
                        edges(count,:)=[index(k) index(j)];
                        count=count+1;
                    end
                end
            end    
        end

         if size(edges,1)>0
             [~,ia]=unique(edges,'rows');
             edges=edges(ia,:);
             badcell=unique(edges(:));
             clear oneBadIsGood
             LCCall=LargestConnectedComponents(edges);
             for i=1:length(LCCall)
                 oneBadIsGood(i)=LCCall{i}(1);
             end
             only_goodcellindex= setdiff([1:size(centroid,1)]', badcell);
             goodcellindex=union(only_goodcellindex,oneBadIsGood);
         else
             goodcellindex=[1:size(centroid,1)]';
         end

end



function LCC=LargestConnectedComponents(edges)
        cellIds=unique(edges(:));
        for j=1:length(cellIds)
            old2new(cellIds(j),1)=j;
            new2old(j,1)=cellIds(j);    
        end
        [length(old2new),length(new2old),length(cellIds)];

        for i=1:size(edges,1)
            for j=1:2 
                newedgename(i,j)= old2new(edges(i,j));
            end
        end
        G=graph(newedgename(:,1),newedgename(:,2));
        bins=conncomp(G);
        % number of connected components 
        nocomp=unique(bins);
        %disp(['# of connected components  ', num2str(length(nocomp))]);
        for i=1:length(nocomp)
            numberOfObjectsInConnectedComponents(i)=sum(nocomp(i)==bins);
        end
        
        
        [sa,sb]=sort(numberOfObjectsInConnectedComponents,'descend');
        
        index=1;
        for i=1:length(sa)
            if sa(i)>1
                LCCIds=find(bins==nocomp(sb(i)));
                LCC{index}=new2old(LCCIds);
                index=index+1;
            end
        end
end




function [neighbor,neighborList,convexVolume] = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points)
%{
Input:
    N - an n*3 matrix with the [x,y,z] coordinates of each of the n points.
    spacing -  1*3 array with the [x,y,z] physical size of the voxels (if N is already in physical units, just use [1,1,1]).
    delta_spacing - a positive real value that specifies the frequency of the sampling of the 3D grid in the output image (the density space / the output argument V).
    exclude_boundary_points - a boolean value that allows the user to choose whether to include or exclude boundary vertices to avoid potential noise in the
    physical boundaries of the sample, default is 'false' (i.e. include boundaries)
Output:
    V - the 3D matrix of estimated densities.
    X, Y, Z - the coordinates of the interpolated space V.
Example:
    N = rand(100, 3);
    spacing = [1 1 1];
    delta_spacing = 0.01;
    exclude_boundary_points = false;
    V = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points);
    figure;
    imagesc(V(:,:,50));
    axis image;
    colormap jet;
%}
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
TRI = delaunay(N(:,1), N(:,2), N(:,3));
[~,convexVolume]=convexHull(delaunayTriangulation(N));

clear neighbor
for i = 1 : size(N,1)
    temp=[];
    for j=1:size(TRI,1)
        for k=1:size(TRI,2)
            if TRI(j,k)==i
                temp=[temp,TRI(j,:)];
            end
        end
    end

   
    neighborList{i}=setdiff(unique(temp),i);
    %neighbor(i,1)=length(neighborList{i});
    ids= neighborList{i};
    dist=pdist(N(ids,1:3));    
    neighbor(i,1)=min(dist);
end




%triplot(tri,x,y); for 2d 
%tetramesh(tri,X); 3d 

[size(N,1), length(unique(TRI(:)))];






end
         
