
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

     

h1=figure();
set(gcf, 'PaperSize', [13 10]);
set(gcf, 'PaperPosition', [0 0 13 10]);

count=1;
for bonetype=1:5
        path=allpath{bonetype};
        clear col1avg
        clear col2avg
        clear col3avg
        clear col4avg
	       
        subplot(2,3,bonetype)
        for fi=1:length(path)
                b1=load(['../../',path{fi},'all_cells.mat']);
                cel=b1.all_cells;
                clear celcent

                s=strsplit(path{fi},'Nuclei_and_Cells_');
                %input2=strcat('Columnar_Structure_Prediction_8_100/',s{2});
                %doubletIndex=load([input2,'LCC_final.dat']);
                %input1=strcat('MakeListColumnarStructurePrediction/',s{2});    
                %load([input1,'centroid_and_surface_cells.mat']);

                %temp=unique(doubletIndex(:));

                if (bonetype==3)|(bonetype==1)
                    celcent(:,1:3)=[cel(:,5),cel(:,6),-cel(:,7)];
                    %coordinate=[centroid(temp,1),centroid(temp,2),-centroid(temp,3)];

                else
                    celcent(:,1:3)=[cel(:,5),cel(:,6),cel(:,7)];
                    %coordinate=[centroid(temp,1),centroid(temp,2),centroid(temp,3)];
                end


                meanofallcell=mean(celcent);
                celcent=celcent-meanofallcell;
                %doubletCoordinates=coordinate-meanofallcell; 
%                 plot(doubletCoordinates(:,1),doubletCoordinates(:,2),'b.');
%                 hold on 
                data(count,:)=[meanofallcell,min(celcent),max(celcent)];
                count=count+1;
        end


end


% axis image 
% view(90,0);
% xlabel('X');
% ylabel('y');
% zlabel('z');
saveas(h1,['DoubletPositionsIn2D.png']);


function [positioninGrowthPlate,centroidZ,Xinterval,youtput]= makeProfile(centroidZ,data,profilesize,zpositionofcolumn)
                %centroidZ=centroid(:,3);
                minl=min(centroidZ); maxl=max(centroidZ);
                Yinterval=linspace(minl-0.00001,maxl+0.00001,profilesize);
                Xinterval=zeros(1,length(Yinterval));
                binsize=Yinterval(2)-Yinterval(1);
                horizontalCount=cell(1,length(Yinterval));
		positioninGrowthPlate=zeros(size(zpositionofcolumn));
                for k=1:length(Yinterval)
                    horizontalCount{k}=[];
                    for tt=1:length(zpositionofcolumn)
                        if (zpositionofcolumn(tt) > Yinterval(k)) & (zpositionofcolumn(tt) <= Yinterval(k+1))
                            horizontalCount{k}=[horizontalCount{k},tt];
				positioninGrowthPlate(tt,1)=k;
                        end
                    end
                    Xinterval(k)=mean(data(horizontalCount{k}));
                end
                youtput=Yinterval+binsize/2;
end
