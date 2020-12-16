

clear all 

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

meanOfAllCelltemp=load('meanOfAllCell.dat');
count=1;
for gi=1:length(allpath)
	for gj=1:length(allpath{gi})
        GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,1:3);
        count=count+1;
    end
end



for gi=1%:length(allpath)
    bonetype=gi;
	for gj=1:length(allpath{gi})
		path=allpath{gi}{gj};
        s=strsplit(path,'Nuclei_and_Cells_');
        %ns=strsplit(path,'T_');
        input2=strcat('Cluster_Structure_Prediction/',s{2});
        
        if ~exist([input2],'dir')
              mkdir([input2]);
        end
        
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
        
        
        
        
        input1=strcat('MakeListColumnarStructurePrediction/',s{2});
        load([input1,'centroid_and_surface_cells.mat']);
        
%         a1=load(['./../../',path,'all_cells_nuclei.mat']);
%         nuc=a1.all_cells_nuclei; 
%         nuccent=nuc(:,5:7);
        a1=load([input1,'NeighboringCell_in_20_micron_cube.dat']);
        edges=a1(:,[1,2]); [~,ia]=unique(edges,'rows');
        edges=edges(ia,:);
        disp('start LCC')
        LCCall=LargestConnectedComponents(edges);
        
        
        makeClusterTifFiles(input2,tile,coordinates,path,LCCall,unique_tileid,bonetype,centroid,  GlobalCenter{gi}{gj});
        plotClusterFigure(input2,LCCall,centroid,nuc);
       


    end
end





function plotClusterFigure(olddirectory,LCCall,centroid,csurfVertices)

    directory=strcat(olddirectory,'Figure_matlab/');

      if ~exist([directory],'dir')
              mkdir([directory]);
      end
      

    for k=1:length(LCCall)
        LCC=LCCall{k};
        if length(LCC)>=4
            h=figure;

            mu=centroid(LCC,:);
            %t=C.coords(LCC,:);
            plot3(mu(:,1),mu(:,2),mu(:,3),'bo-');
            hold on 
            for i=1:size(LCC,1)
                %plot_cube(t(i,[3,4,1,2,5,6]).*voxelimit);
                V=csurfVertices{LCC(i)};
                if mod(i,2)==0
                    plot3(V(:,1),V(:,2),V(:,3),'r.','markersize',0.1);
                else
                    plot3(V(:,1),V(:,2),V(:,3),'g.','markersize',0.1);
                end
                text(mu(i,1),mu(i,2),mu(i,3),num2str(LCC(i)),'fontsize',16);
                hold on 
            end
            saveas(h,[directory,'Fig','_',num2str(k),'_',num2str(length(LCC)),'.fig']);
            saveas(h,[directory,'Fig','_',num2str(k),'_',num2str(length(LCC)),'.png']);
            close all 
        end
    end
    
end


function makeClusterTifFiles(olddirectory,tile,coordinates,path,LCCall,unique_tileid,bonetype,centroid,GlobalCenter)

      
        if (bonetype==3)|(bonetype==1)
                tempCentroid=[centroid(:,1:2),-centroid(:,3)];
        else
                tempCentroid=[centroid(:,1:2),centroid(:,3)];
        end
        newcentroid=tempCentroid-GlobalCenter;


      directory=strcat(olddirectory,'Figure_tif/');
      if ~exist([directory],'dir')
              mkdir([directory]);
      end

     for k=1:length(LCCall)
            LCC=LCCall{k};
            if length(LCC)>=4
                tile_frequency_in_col=unique_tileid(LCC,1);
                aa=unique(tile_frequency_in_col); 
                clear freq 
                for jj=1:length(aa)
                    freq(jj)=sum(tile_frequency_in_col==aa(jj));
                end
                freq=freq/sum(freq);
                [sa,sb]=max(freq); 
                complete=0;
                if freq(sb)==1
                   complete=1;
                end
                   
                fi=aa(sb);
                position=coordinates(fi,1);
                
                load(strcat('./../../',path,'c_n_pos',num2str(position),' (Characteristics).mat'));
                
                zsize1=double(max(C.coords(:,6)));
                zsize2=ceil(abs(tile(fi,4)-tile(fi,3))/0.387);
                zsize=max(zsize1,zsize2);
                [k,zsize1,zsize2];
                
                
                
                I.size= [900, 900, zsize];  I.data_type='uint8';
                characteristics = zeros(I.size, I.data_type);
                
                p1=round(min(newcentroid(LCC,3)));
                p2=round(max(newcentroid(LCC,3)));

		if p1<-50
                
                    for j=1:length(LCC)
                        if fi==unique_tileid(LCC(j),1)
                            i=unique_tileid(LCC(j),2);
                            %[i,position,fi,k]
                            pos=C.coords(i,:);
                            characteristics(pos(1):pos(2),pos(3):pos(4),pos(5):pos(6))= characteristics(pos(1):pos(2),pos(3):pos(4),pos(5):pos(6)) | C.masks{i};
                        end
                    end                      
                    tifname=strcat([directory,'/Output_cluster_', 'w_', num2str(complete), '_cid_', num2str(k),'_com_',num2str(length(LCC)), '_til_',num2str(fi),'_pos_', num2str(p1),'-',num2str(p2)],'.tif');  
		                    
                    saveTiffFiles(tifname,I,255*characteristics)  

		end
                  
            end
        
     end
            
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








