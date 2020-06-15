
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
        %GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,1:3);
        GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,[6,9]);
        count=count+1;
    end
end



data=[];
 


for gi=1:length(allpath)
    bonetype=gi;
    sampleROC=[];
	for gj=1:length(allpath{gi})
        
              
        path=allpath{gi}{gj};
        s=strsplit(path,'Nuclei_and_Cells_');
        input1=strcat('MakeListColumnarStructurePrediction/',s{2});
        input2=strcat('Cluster_Structure_Prediction/',s{2});
        

        %outputpath=strcat('Columnar_Structure_Prediction_8_100/',s{2});
        %load(['MakeListColumnarStructurePrediction/',s{2},'centroid_and_surface_cells.mat'],'centroid');
        
%         a1=load([input1,'NeighboringCell_in_20_micron_cube.dat']);
%         edges=a1(:,[1,2]); [~,ia]=unique(edges,'rows');
%         edges=edges(ia,:);
%         disp('start LCC')
%         LCCall=LargestConnectedComponents(edges);
        
        pureColumn=load(['degree_of_the_column/StraightLine/pureColum_',s{2}(1:strlength(s{2})-1),'.dat']);
        
        ROC=load(['RadiusOfCurvature/ROC_',s{2}(1:strlength(s{2})-1),'.dat']);
        
        
        
        limit=GlobalCenter{gi}{gj};
        starting=limit(1)+ 0.2*(limit(2)-limit(1));
        finishing=limit(1)+ 0.6*(limit(2)-limit(1));
        [limit,starting,finishing];
        
        
        
        
        inverseR=[];
        for i=1:length(pureColumn)
            for j=1:length(ROC)
                if pureColumn(i)==ROC(j,1)
                    if ROC(j,2)>=4
                        if ((starting<=ROC(j,5)) & (ROC(j,6)<=finishing))
                                inverseR=[inverseR;ROC(j,3:4)];
                        end
                    end
                end
            end
        end
        
        
       
        
        t=find(~isinf(inverseR(:,2)));
        [gi,gj,length(ROC), length(pureColumn),length(inverseR),length(t)]
          Rsq=corrcoef(inverseR(t,1),inverseR(t,2)); [mean(inverseR(:,1)), Rsq(1,2)]
        
    end
    %gpROC{gi}=sampleROC;
end
        

% 
% for i=[1,3]
%     [y1,x1]=histnorm(gpROC{i},10);
%     plot(x1,y1)
%     hold on 
% end
% 
% 



     
