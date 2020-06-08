
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
 


for gi=3%:length(allpath)
    bonetype=gi;
	for gj=1:length(allpath{gi})
        
        h=figure;
              
        path=allpath{gi}{gj}
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
        LCCall=LargestConnectedComponents(edges);
        
        
    
        a1=load(['./../../',path,'all_cells_nuclei.mat']);
        nucleiFeatures=a1.all_cells_nuclei; 
        nuccent=nucleiFeatures(:,5:7);
             
        C=mean(nuccent); zgap=max(nuccent(:,3))-min(nuccent(:,3));
        nuccent1=nuccent-C;
        
        s=strsplit(path,'Nuclei_and_Cells_');
        outputpath=strcat('Columnar_Structure_Prediction_8_100/',s{2});
       
      
        if (bonetype==3)|(bonetype==1)
        %zmin=-zmin;  zmax=-zmax;
        bonetypeMatrix=[1,0,0;0 1 0;0, 0,-1];
        nuccent1=nuccent1*bonetypeMatrix;
        else
        %zmin=zmin+70;  zmax=zmax+70; 
        bonetypeMatrix=[1,0,0;0 1 0;0, 0,1];
        nuccent1=nuccent1*bonetypeMatrix;
        end

        [CurvatureAxisLine,~]=movingMeanAverage(nuccent1);
        dlmwrite(['Curvature/CurvLine_',s{2}(1:strlength(s{2})-1),'.dat'],CurvatureAxisLine,'\t');
        
        
        %plot3(nuccent(:,1),nuccent(:,2),-nuccent(:,3),'b.')
%         if mod(bonetype,2)==0
%                  plot3(nuccent1(:,1),nuccent1(:,2),nuccent1(:,3),'b.')
%         else
%                  plot3(nuccent1(:,1),nuccent1(:,2),nuccent1(:,3),'b.')
%         end
%         axis image 
%         xlabel('X')
%         ylabel('Y')
%         zlabel('P-D axis')
        
        
        
       plotClusterFigure(input2,LCCall,centroid,nuc, C, bonetypeMatrix,'r.');

     
        
       
       
       shp = alphaShape(nuccent1);
	   [tetrahedron,V] = alphaTriangulation(shp);
       trep = triangulation(tetrahedron, V);
       [tri, V] = freeBoundary(trep);
        
        plot3( V(:,1), V(:,2), V(:,3),'b.','markersize',0.01);
        hold on 
        
        plot3( CurvatureAxisLine(:,1), CurvatureAxisLine(:,2), CurvatureAxisLine(:,3),'r.-');
        %plot3( CurvatureAxisLineSmooth(:,1), CurvatureAxisLineSmooth(:,2), CurvatureAxisLineSmooth(:,3),'g.-');

        % hold on 
        % plot3(Xdata,Ydata,Zdata,'r-')
        title('Curvature axis line','fontweight','normal')
        axis image
        %view(36,14)
        xlabel('X')
        ylabel('Y')
        zlabel('P-D axis')
        %zlabel('Proximal-Distal Axis')
        view(-38,30)
        saveas(h,['Curvature/curvatureAxisLine_',s{2}(1:strlength(s{2})-1),'.png'])
        saveas(h,['Curvature/curvatureAxisLine_',s{2}(1:strlength(s{2})-1),'.fig'])
        close all 
        
        end
        
        
        
        
     
end






function [CurvatureAxisLine,dd]=movingMeanAverage(nuccent)
    zmin=min(nuccent(:,3)); zmax=max(nuccent(:,3));
    interval=zmin-0.001:1:zmax+0.001;
    count=1;
     for i=1:length(interval)
            left=interval(i)-100;
            right=interval(i)+100;
            index=find((nuccent(:,3)>=left) & (nuccent(:,3)<right));
            if length(index)>=100
                CurvatureAxisLine(count,:)=mean(nuccent(index,:));
                count=count+1;
            end
     end
     
%        d=CurvatureAxisLine;
%        CS = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));
%        dd = interp1(CS, d, unique([CS(:)' linspace(0,CS(end),100)]),'pchip');
         dd=[];
     
     
end



function plotClusterFigure(olddirectory,LCCall,centroid,csurfVertices,globalCenter, bonetypeMatrix,mycolor)

      directory=strcat(olddirectory,'Figure_matlab/');

      if ~exist([directory],'dir')
              mkdir([directory]);
      end
      

    for k=1:length(LCCall)
        LCC=LCCall{k};
        if length(LCC)>=4
            mu=(centroid(LCC,:)-globalCenter)* bonetypeMatrix;
            %plot3(mu(:,1),mu(:,2),mu(:,3),'bo-');
            hold on 
            clear csurf
            for i=1:size(LCC,1)
                %plot_cube(t(i,[3,4,1,2,5,6]).*voxelimit);
                V=(csurfVertices{LCC(i)}-globalCenter)* bonetypeMatrix;
                csurf{i}=V;
               
                %plot3(V(:,1),V(:,2),V(:,3),mycolor,'markersize',0.5);
              
                %text(mu(i,1),mu(i,2),mu(i,3),num2str(LCC(i)),'fontsize',16);
                hold on 
            end
            
             smoothUsingAverage(csurf,'g.-');
             %smoothUsingEllipsoid(csurf,'b.-');
            
            
            %saveas(h,[directory,'Fig','_',num2str(k),'_',num2str(length(LCC)),'.fig']);
            %saveas(h,[directory,'Fig','_',num2str(k),'_',num2str(length(LCC)),'.png']);
            %close all 
        end
    end
    
    
    
   
          
    
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
        
         d = d( ~any( isnan( d ) | isinf( d ), 2 ),: );
     
         CS = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));
         dd = interp1(CS, d, unique([CS(:)' linspace(0,CS(end),100)]),'pchip');
         plot3(dd(:,1),dd(:,2),dd(:,3),color,'linewidth',5)
end






function [data,centroidZ,Xinterval,youtput]= makeProfile(centroidZ,data,profilesize)
                %centroidZ=centroid(:,3);
                minl=min(centroidZ); maxl=max(centroidZ);
                Yinterval=linspace(minl-0.00001,maxl+0.00001,profilesize);
                Xinterval=zeros(1,length(Yinterval));
                binsize=Yinterval(2)-Yinterval(1);
                horizontalCount=cell(1,length(Yinterval));
                for k=1:length(Yinterval)
                    horizontalCount{k}=[];
                    for tt=1:length(centroidZ)
                        if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) <= Yinterval(k+1))
                            horizontalCount{k}=[horizontalCount{k},tt];
                        end
                    end
                    Xinterval(k)=sum(data(horizontalCount{k}));
                end
                youtput=Yinterval+binsize/2;
end


function b=binarySum(a)
          for i=1:size(a,1)
              if ((a(i,1)==1) & (a(i,2)==1))
                  b(i)=1;
              elseif ((a(i,1)==0) & (a(i,2)==1))
                  b(i)=2;
              elseif ((a(i,1)==0) & (a(i,2)==0))    
                  b(i)=3;
              elseif ((a(i,1)==1) & (a(i,2)==0))
                  b(i)=4;
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

