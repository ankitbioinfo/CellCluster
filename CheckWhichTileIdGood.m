
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



for gi=4%:length(allpath)
    bonetype=gi;
	for gj=4%1:length(allpath{gi})
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
        
        
        for position = coordinates(:,1)'  % intersect(PZ{gi}{gj},coordinates(:,1)')        %setdiff(coordinates(:,1)',RZ{gi}{gj}) 
                load(strcat('./../../',path,'c_n_pos',num2str(position),' (Characteristics).mat'));
                indtemp = find(G.inter.volume_ratio>1);
                fi = find(coordinates(:,1) == position);
                for i=1:length(indtemp)
                    j=indtemp(i);
               
                    cent2=G.cel.centroids(j,:);   
                    cent1=mean(C.surfaces(j).vertices);
                    d=pdist([cent1;cent2]);
                    
                    ankit(i)=d;
                    
                end
                [position,mean(ankit)]
                clear ankit 
        end
                      
        
      
       


    end
end
