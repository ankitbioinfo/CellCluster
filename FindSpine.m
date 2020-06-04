% 
% 
% clear all 
% 
% dt_path_wt={ '../data/Nuclei_and_Cells_DT_S18_m6_wt/', '../data/Nuclei_and_Cells_DT_S17_m2_wt/',...  
%              '../data/Nuclei_and_Cells_DT_S84_m3_wt/', '../data/Nuclei_and_Cells_DT_S51_m2_wt/',...
%              '../data/Nuclei_and_Cells_DT_S84_m4_wt/'};
% 
% pt_path_wt = {  '../data/Nuclei_and_Cells_PT_S18_m6_wt/','../data/Nuclei_and_Cells_PT_S17_m2_wt/',...  
%                 '../data/Nuclei_and_Cells_PT_S84_m3_wt/','../data/Nuclei_and_Cells_PT_S51_m2_wt/',...
%                 '../data/Nuclei_and_Cells_PT_S84_m4_wt/'};
% 
% dt_path_mut= {'../data/Nuclei_and_Cells_DT_S17_m1_mut/', '../data/Nuclei_and_Cells_DT_S18_m2_mut/' ,...
%               '../data/Nuclei_and_Cells_DT_S84_m1_mut/', '../data/Nuclei_and_Cells_DT_S84_m5_mut/'};
%     
% pt_path_mut = {'../data/Nuclei_and_Cells_PT_S17_m1_mut/', '../data/Nuclei_and_Cells_PT_S18_m2_mut/',...
%                '../data/Nuclei_and_Cells_PT_S84_m1_mut/', '../data/Nuclei_and_Cells_PT_S84_m5_mut/', };
%                
% du_path_wt={'../data/Nuclei_and_Cells_DU_S51_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m3_wt/'};
% 
%        
% allpath={dt_path_wt; pt_path_wt; dt_path_mut; pt_path_mut};  
% 
% meanOfAllCelltemp=load('meanOfAllCell.dat');
% count=1;
% for gi=1:length(allpath)
% 	for gj=1:length(allpath{gi})
%         GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,1:3);
%         count=count+1;
%     end
% end
% 
% 
% 
% for gi=3%:length(allpath)
%     bonetype=gi;
% 	for gj=1%1:length(allpath{gi})
% 		path=allpath{gi}{gj};
%         s=strsplit(path,'Nuclei_and_Cells_');
%         %ns=strsplit(path,'T_');
%         input2=strcat('Cluster_Structure_Prediction/',s{2});
%         
%         if ~exist([input2],'dir')
%               mkdir([input2]);
%         end
%         
%            alignment=load(['./../../',path,'Alignment_matrix.dat']);
%             vec=alignment(:,1:3);
%              [numbers,txt,raw] = xlsread(['./../../',path,'Tile_coordinates.xlsx']);
%                coordinates = zeros(size(txt,1)-3,5);
%                 for i = 4:size(txt,1),
%                     temp =  char(txt(i,1));
%                     res = strsplit(temp,'_POS');
%                     coordinates(i-3,1) = str2num(char(res(2)));
%                     coordinates(i-3,2:5) = numbers(i-3,:);
%                 end
%                 tile=coordinates(:,2:end);
%         
%         
%         
%         
%         input1=strcat('MakeListColumnarStructurePrediction/',s{2});
%         load([input1,'centroid_and_surface_cells.mat']);
%         
% %         a1=load(['./../../',path,'all_cells_nuclei.mat']);
% %         nuc=a1.all_cells_nuclei; 
% %         nuccent=nuc(:,5:7);
%         a1=load([input1,'NeighboringCell_in_20_micron_cube.dat']);
%         edges=a1(:,[1,2]); [~,ia]=unique(edges,'rows');
%         edges=edges(ia,:);
%         disp('start LCC')
%         LCCall=LargestConnectedComponents(edges);
%         
%         
%         %makeClusterTifFiles(input2,tile,coordinates,path,LCCall,unique_tileid,bonetype,centroid,  GlobalCenter{gi}{gj});
%       
%        
% 
% 
%     end
% end




  plotClusterFigure(input2,LCCall,centroid,nuc);


function plotClusterFigure(olddirectory,LCCall,centroid,csurfVertices)

    directory=strcat(olddirectory,'Figure_matlab/');

      if ~exist([directory],'dir')
              mkdir([directory]);
      end
      

    for k=9%1:length(LCCall)
        LCC=LCCall{k};
        if length(LCC)>=4
            h=figure;

            mu=centroid(LCC,:);
            %t=C.coords(LCC,:);
            plot3(mu(:,1),mu(:,2),mu(:,3),'bo-');
            hold on 
            clear csurf
            for i=1:size(LCC,1)
                %plot_cube(t(i,[3,4,1,2,5,6]).*voxelimit);
                V=csurfVertices{LCC(i)};
                csurf{i}=V;
                if mod(i,2)==0
                    plot3(V(:,1),V(:,2),V(:,3),'b.','markersize',0.5);
                else
                    plot3(V(:,1),V(:,2),V(:,3),'r.','markersize',0.5);
                end
                text(mu(i,1),mu(i,2),mu(i,3),num2str(LCC(i)),'fontsize',16);
                hold on 
            end
            %saveas(h,[directory,'Fig','_',num2str(k),'_',num2str(length(LCC)),'.fig']);
            %saveas(h,[directory,'Fig','_',num2str(k),'_',num2str(length(LCC)),'.png']);
            %close all 
        end
    end
    
    
    
    view(70,31)
    smoothUsingAverage(csurf,'r.-');
    smoothUsingEllipsoid(csurf,'b.-');
          
    
end



function smoothUsingAverage(surface,color)

        csurf=[];
        for i=1:length(surface)
            csurf=[csurf;surface{i}];
        end


        z=min(csurf(:,3)):2:max(csurf(:,3));
        for i=1:length(z)-1
              index=find( (csurf(:,3)>z(i)) &  (csurf(:,3)<=z(i+1)));
              d(i,:)=mean(csurf(index,:));

        end
    
        %plot3(d(:,1),d(:,2),d(:,3),'kx-')
     
         CS = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));
         dd = interp1(CS, d, unique([CS(:)' linspace(0,CS(end),100)]),'pchip');
         plot3(dd(:,1),dd(:,2),dd(:,3),color)
end



function smoothUsingEllipsoid(surface,color)

    warning 'off'
    for i=1:length(surface)
        [ center, radii, evecs, v, chi2 ] = ellipsoid_fit_new( surface{i} ); %radii
		radii=[12,6,3];
        [x,y,z] = ellipsoid(0,0,0,radii(1),radii(2),radii(3),40);
        tt=[x(:),y(:),z(:)]*inv(evecs);
%         nx=reshape(tt(:,1),size(x,1),size(x,2));
%         ny=reshape(tt(:,2),size(x,1),size(x,2));
%         nz=reshape(tt(:,3),size(x,1),size(x,2));
        ellipsoidalCell{i}=tt+center';
        %ankit{i}=evecs; 
        %plot3(ellipsoidalCell{i}(:,1),ellipsoidalCell{i}(:,2),ellipsoidalCell{i}(:,3),'r.');
    end
    
    smoothUsingAverage(ellipsoidalCell,color);
    
end








function LCC=LargestConnectedComponents(edges)
        
        cellIds=unique(edges(:));
        for j=1:length(cellIds)
            old2new(cellIds(j),1)=j;
            new2old(j,1)=cellIds(j);    
        end
        [length(old2new),length(new2old),length(cellIds)];

        for i=1:length(edges)
            for j=1:2 
                newedgename(i,j)= old2new(edges(i,j));
            end
        end
        G=graph(newedgename(:,1),newedgename(:,2));
        bins=conncomp(G);
        % number of connected components 
        nocomp=unique(bins);
        disp(['# of connected components  ', num2str(length(nocomp))]);
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




function saveTiffFiles(tifoutputname,I,characteristics)
    t = Tiff(tifoutputname,'w');
    tagstruct.ImageLength = I.size(1);
    tagstruct.ImageWidth =  I.size(2);
    tagstruct.SubFileType= Tiff.SubFileType.Page;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt; % uint
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.Compression = Tiff.Compression.PackBits;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    for ii=1:I.size(3)
       setTag(t,tagstruct);
       write(t,characteristics(:,:,ii));
       writeDirectory(t);
    end
    close(t)
end








