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




%dlmwrite(inputfile,neighbor,'\t')

%triplot(tri,x,y); for 2d 
%tetramesh(tri,X); 3d 

[size(N,1), length(unique(TRI(:)))];


% fid=fopen(inputfile,'w');
% disp('done')
% 
% for i = 1 : size(N,1)
% 	neighbor(i)=find(any(TRI == i,2));
% 	for j=1:length(neighbor)
% 		fprintf(fid,'%d\t%d\n',i,neighbor(j));
% 	end 
% 	
% 	%sum(tetra_vol(any(TRI == i,2)))
%     	%point_density(i,1) = 4/sum(tetra_vol(any(TRI == i,2)));
% end

%fclose(fid);

%point_density(5,1)
%convexVolume




end
