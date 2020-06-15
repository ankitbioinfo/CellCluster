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
 


for gi=5%:length(allpath)
    bonetype=gi;
	for gj=3%:length(allpath{gi})
        
        [gi,gj]
              
        path=allpath{gi}{gj};
        s=strsplit(path,'Nuclei_and_Cells_');
        input1=strcat('MakeListColumnarStructurePrediction/',s{2});
        input2=strcat('Cluster_Structure_Prediction/',s{2});
        

        load(['degree_of_the_column/graphlet_',s{2}(1:strlength(s{2})-1),'.mat']);
        
        
        
        idlist=[1, 2, 3, 4, 10, 7, 8, 9, 26, 21, 23, 97, 90, 661];
        
        for fi=1:length(idlist)-1
                clear old2new
                clear new2old 
                clear newedgename
            
                edges=unique(graphlet{idlist(fi)},'rows');
                cellIds=unique(edges(:));
                for j=1:length(cellIds)
                    old2new(cellIds(j),1)=j;
                    new2old(j,1)=cellIds(j);    
                end


                for i=1:length(edges)
                    for j=1:2 
                        newedgename(i,j)= old2new(edges(i,j));
                    end
                end



                G=graph(newedgename(:,1),newedgename(:,2));

                h=figure;
                set(gcf, 'PaperSize', [5 5]); %7
                set(gcf, 'PaperPositionMode', 'manual');
                set(gcf, 'PaperPosition', [0 0 5 5]);
                
                
                plot(G,'Layout','force','linewidth',3)
                set(findobj(gcf,'type','axes'),'Visible','off')
                
                
              
                
                
                
                print(gcf, '-dpng', ['Visualize/bad/graphlet_',s{2}(1:strlength(s{2})-1),'_',num2str(idlist(fi))]);
                RemoveWhiteSpace([], 'file', ['Visualize/bad/graphlet_',s{2}(1:strlength(s{2})-1),'_',num2str(idlist(fi)),'.png'] ,'output', ['Visualize/graphlet_',s{2}(1:strlength(s{2})-1),'_',num2str(idlist(fi)),'.png']  );

                %exportgraphics(t,'fourplots.pdf','BackgroundColor','none')
                close all 
        end
    end
end







