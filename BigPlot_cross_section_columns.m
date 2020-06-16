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
% allpath={dt_path_wt; pt_path_wt; dt_path_mut; pt_path_mut; du_path_wt};  
% 
% 
% meanOfAllCelltemp=load('meanOfAllCell.dat');
% count=1;
% for gi=1:length(allpath)
% 	for gj=1:length(allpath{gi})
%         %GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,1:3);
%         GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,[6,9]);
%         count=count+1;
%     end
% end
% 
% 
% 
% for gi=1:length(allpath)
%     	bonetype=gi;
% 	for gj=1:length(allpath{gi})
% 		path=allpath{gi}{gj};
%         s=strsplit(path,'Nuclei_and_Cells_');
%         input2=strcat('Columnar_Structure_Prediction_8_100/',s{2});
%         
%         pureColumn=load(['degree_of_the_column/StraightLine/pureColum_',s{2}(1:strlength(s{2})-1),'.dat']);
%         column_curvature=load(['Curvature/Column_curvature_',s{2}(1:strlength(s{2})-1),'.mat']) ;
%         ROC=load(['RadiusOfCurvature/ROC_',s{2}(1:strlength(s{2})-1),'.dat']);
% 
%         
%         a1=load(['./../../',path,'all_cells_nuclei.mat']);
%         nuc=a1.all_cells_nuclei; 
%         nuccent=nuc(:,5:7);
%         
%         C=mean(nuccent); %zgap=max(nuccent(:,3))-min(nuccent(:,3));
%         
%         limit=GlobalCenter{gi}{gj};
%         starting=limit(1)+ 0.2*(limit(2)-limit(1));
%         finishing=limit(1)+ 0.6*(limit(2)-limit(1));
%         [limit,starting,finishing];
%         
%         
%         
%         nuccent1=nuccent-C;
%         fixzmin=starting; fixzmax=finishing;
%         if (bonetype==3)|(bonetype==1)
%         zmin=-fixzmin;  zmax=-fixzmax;
%         colzmin=fixzmin;colzmax=fixzmax;
%         else
%         zmin=fixzmin;  zmax=fixzmax;  %70 add 
%         colzmin=zmin;  colzmax=zmax;
%         end
%         
%         shp = alphaShape(nuccent1);
% 	   [tetrahedron,V] = alphaTriangulation(shp);
%        trep = triangulation(tetrahedron, V);
%        [tri, V] = freeBoundary(trep);
%        
%         if zmin<zmax
%         PZboundary=find(V(:,3)>zmin & V(:,3)<zmax);
%         nuc_PZ_boundary=find(nuccent1(:,3)>zmin & nuccent1(:,3)<zmax);
%         else
%         PZboundary=find(V(:,3)<zmin & V(:,3)>zmax); 
%         nuc_PZ_boundary=find(nuccent1(:,3)<zmin & nuccent1(:,3)>zmax);
%         end
%         
%         x=V(PZboundary,1); y=V(PZboundary,2);
%         K=boundary(x,y);
%         databoundary{gi}{gj}=[x(K),y(K)];
%         
%         
%         index=1;
%         inverseR=[]; 
%        for i=1:length(pureColumn)
%             for j=1:length(ROC)
%                 if pureColumn(i)==ROC(j,1)
%                     if ROC(j,2)>=3
%                         if ((starting<=ROC(j,5)) & (ROC(j,6)<=finishing))
%                                 inverseR=[inverseR;ROC(j,3:4)];
%                                 twocell{index,1}=column_curvature.CurvatureOfColumn.NormalizedColumnCentroid{pureColumn(i)};
%                                 index=index+1;
%                         end
%                     end
%                 end
%             end
%        end
%         
%         
%        
%       
%         curvatureColormap{gi}{gj}=inverseR(:,1);
%         nucleiDensity{gi}{gj}=nuccent1( nuc_PZ_boundary,1:2);
%         data{gi}{gj}=twocell;
%         %slope{gi}{gj}=sampleSlope;
%         clear twocell
%         %clear sampleSlope 
%         [size(data{gi}{gj},1), size(  curvatureColormap{gi}{gj},1)     ];
%     end
% end





name{2}={ 'S18 m6','S17 m2', 'S84 m3', 'S51 m2', 'S84 m4'  };
name{1}={ 'S18 m6','S17 m2', 'S84 m3', 'S51 m2', 'S84 m4' };  
    
name{3}={ 'S17 m1', 'S18 m2', 'S84 m1', 'S84 m5' }; 
name{4}={ 'S17 m1', 'S18 m2', 'S84 m1', 'S84 m5' }; 
name{5}={ 'S51 m2', 'S84 m2', 'S84 m3'};

tname{5}='DU WT'; tname{4}='PT MT'; tname{3}='DT MT'; tname{2}='PT WT'; tname{1}='DT WT';



mycolor={'r','b','k','c','g'};

h1=figure();
XL=0.07;XR=0.03;XGap=0.02;Row=5;
YT=0.06;YB=0.08;YGap=0.05;Col=5;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [15 12]);
set(gcf, 'PaperPosition', [0 0 15 12]);


U=[1,0,0];      mymap=colormap(redbluecmap(3));

 map3=[0 0 1; 0 1 0; 1 0 0];
 map2=[0 0 1; 0 1 0];
 map4=[0 0 1; 0 1 0; 1 0 0;1 1 0];
 map5=[0 0 1; 0 1 0; 1 0 0;1 1 0;1 0 1];


 curvatureRange=linspace(0,200,50);
 
 
 count=1;


for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        if j<=length(data{i})
                subplot('Position',marray);
                column_coordinate=data{i}{j};
                pzboundary=databoundary{i}{j};
                nuclei=nucleiDensity{i}{j};
                curvatureValue=curvatureColormap{i}{j};
               
                
                
                mudata=[];
                %[i,j,size(yr)]
                for k=1:size(column_coordinate,1)
                       PC=pca(column_coordinate{k});
                       GC=mean(column_coordinate{k});
                       V=column_coordinate{k};
                       mudata(k,:)=GC(1:2);
                       colLen(k,1)=size(column_coordinate{k},1);
                end
                
                C=mean(nuclei);
                
                
                %pzboundary=pzboundary-C;
                goodmin=min(pzboundary);
                goodmax=max(pzboundary);
                
                
                if size(mudata,1)>1
                    %mudata=mudata-C;
                    goodindex=find((mudata(:,1)>goodmin(1))&(mudata(:,1)<goodmax(1))&(mudata(:,2)>goodmin(2))&(mudata(:,2)<goodmax(2)));
                    curvatureValue=curvatureValue(goodindex,:);
                    colLen= colLen(goodindex,:);
                    [th,r]=cart2pol(mudata(goodindex,1),mudata(goodindex,2));
                    
                    ankit(count,:)=[min(curvatureValue), max(curvatureValue), length(unique(colLen))]; count=count+1;
                    
                    
                    sizeIndex=zeros(size(curvatureValue));
                    for k=2:length(curvatureRange)
                        allInd=find((curvatureValue>curvatureRange(k-1)) & (curvatureValue<=curvatureRange(k) ) );
                        sizeIndex(allInd,1)=curvatureRange(k);
                    end
                    
                    
                    
                    polarscatter(th,r,10*sizeIndex,  colLen ,'marker','.');
                    number= length(unique(colLen));
                    cbh=colorbar; 
                    if number==2
                        cbr=colormap(gca,map2); 
                    elseif number==3
                        cbr=colormap(gca,map3); 
                    elseif number==4
                        cbr=colormap(gca,map4); 
                    else 
                        cbr=colormap(gca,map5);
                    end
                    
                    cbh.Ticks=unique(colLen');
                %slope=atan2(mudata(:,2),mudata(:,1));
                %polarhistogram(slope,12,'FaceColor','b','normalization','probability','FaceAlpha',.3); 
                hold on 
                end
                
            
                [th,r]=cart2pol(pzboundary(:,1),pzboundary(:,2));
                polarplot(th,r,'color','m','marker','.','linestyle','-','linewidth',1);
                
               
                
                
                
                %xr=[yr(:,1:3);yr(:,4:6)];
                set(gca,'fontsize',7);
                %[i,j,max(ankit(:,1)),min(ankit(:,1))];
                %dlmwrite(['GlobalXangle/cosine_ratio_',num2str(i),'_',num2str(j),'.dat'],ankit,'\t');
                %[i,j,100*sum(ankit>75)/length(ankit),100*sum(ankit<15)/length(ankit)];
                %clear ankit 
                
            
            %axis([min(xr(:,1)),max(xr(:,1)),min(xr(:,2)),max(xr(:,2))])
            %axis([-400,400,-400,400])
            title(strcat(tname{i}, '-',name{i}{j}),'fontweight','normal','fontsize',7)
            %legend(p,name{chro},'location','northwest');
            

                         
        end
        
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end


saveas(h1,['column_cs_spatial_.png']);

%print(gcf,'Cross_section_50-250_mu','-dpdf','-r300');



function [xrange,yrange]=make_hist(x)
    interval=[0:10:100];%       [0:2:20,25:5:95];
    y=zeros(1,length(interval)-1);
    for i=2:length(interval)
        y(i-1)=sum(find((x<interval(i)) & (x>=interval(i-1))));
    end
    yrange=(y/sum(y))';
    xrange=(interval(1:end-1))';
  
end


function chooseColor=pickcolor(value,map)
        colormap(map); 
        %map=(colormap(jet(6)));    
        theta=[0,15,30,45,60,75,90];
        %theta=linspace(0,90,7);
        %theta=[0,5,10,15,20,25,30,45,60,75,90];            
        for i=2:length(theta)
            if value<=theta(i)
                chooseColor=map(i-1,:);
                break 
            end
        end
end
    
