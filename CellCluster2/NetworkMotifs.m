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
%allpath={du_path_wt};
  
% the columns contains the individual nuclei features
% 1 - stack id
% 2 - volume
% 3 - surface area
% 4 - sphericity
% 5-7 - centroid x,y,z coordinates
% 8-10 - PC1 x,y,z orientation
% 11-13 - PC2 x,y,z orientation
% 14-16 - PC3 x,y,z orientation
% 17-19 - PC1,PC2,PC3 latent coefficient
% 20 - Delaunay density

mycolor={'r.','b.','g.','m.','k.'};
nucallcolor={'r--','b--','g--','m--','k--'};

fcelallcolor={'ro-','bo-','go-','mo-','ko-'};
fnucallcolor={'ro--','bo--','go--','mo--','ko--'};




meanOfAllCelltemp=load('meanOfAllCell.dat');
count=1;
for gi=1:length(allpath)
	for gj=1:length(allpath{gi})
        GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,1:3);
        count=count+1;
    end
end


data=[];
 


for gi=1:length(allpath)
    bonetype=gi;
	for gj=1:length(allpath{gi})
        
        [gi,gj]
              
        path=allpath{gi}{gj};
        s=strsplit(path,'Nuclei_and_Cells_');
        input1=strcat('MakeListColumnarStructurePrediction/',s{2});
        input2=strcat('Cluster_Structure_Prediction/',s{2});
        
%         if bonetype==5
%                     ns=strsplit(path,'U_');
%         else
%                     ns=strsplit(path,'T_');
%         end
        %outputpath=strcat('Columnar_Structure_Prediction_8_100/',s{2});
        load(['MakeListColumnarStructurePrediction/',s{2},'centroid_and_surface_cells.mat'],'centroid','nuc');
        
        a1=load([input1,'NeighboringCell_in_20_micron_cube.dat']);
        edges=a1(:,[1,2]); [~,ia]=unique(edges,'rows');
        edges=edges(ia,:);
        disp('start LCC')
        [LCCall,degree,graphlet]=LargestConnectedComponents(edges);
        
         
          directory=strcat('degree_of_the_column/Graphlet/');
          if ~exist([directory],'dir')
                  mkdir([directory]);
          end
        
        save([directory,'graphlet_',s{2}(1:strlength(s{2})-1),'.mat'],'graphlet');
        
      
          directory=strcat('degree_of_the_column/degree_sequence/');
          if ~exist([directory],'dir')
                  mkdir([directory]);
          end
          
        fid=fopen([directory,'degree_',s{2}(1:strlength(s{2})-1),'.dat'],'w');
        for i=1:length(degree)
            for j=1:length(degree{i})
                fprintf(fid,'%d ',degree{i}(j));
            end
            fprintf(fid,'\n');
        end
       
        
        
    end
end

% 
% function degree=findDegrees(edges,LCCall)
% 
%         cellIds=unique(edges(:));
%         for j=1:length(cellIds)
%             old2new(cellIds(j),1)=j;
%             new2old(j,1)=cellIds(j);    
%         end
%         [length(old2new),length(new2old),length(cellIds)]
% 
%         for i=1:length(edges)
%             for j=1:2 
%                 newedgename(i,j)= old2new(edges(i,j));
%             end
%         end
%         G=graph(newedgename(:,1),newedgename(:,2));
% 
%        
%     
%         
%         
% end










function [LCC,degree,Graphlet]=LargestConnectedComponents(edges)
        
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
                degree{index}=G.degree(LCCIds);
                
                
                temp=[];
                for k=1:length(LCCIds)
                    for j=1:length(newedgename)
                        for l=1:2 
                            if newedgename(j,l)==LCCIds(k)
                                temp=[temp;j];
                            end
                        end
                    end
                end
              
                Graphlet{index}=newedgename(temp,:);
               
%                 temp=[temp;find(newedgename(:,1)==LCCIds)];
%                 temp=[temp;find(newedgename(:,2)==LCCIds)];
              
                index=index+1;
                
            end
        end
        
       
        
        
        
        
end
