
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
        
        pureColumn=load(['degree_of_the_column/StraightLine/pureColum_',s{2}(1:strlength(s{2})-1),'.dat']);
        
        
    
        a1=load(['./../../',path,'all_cells_nuclei.mat']);
        nucleiFeatures=a1.all_cells_nuclei; 
        nuccent=nucleiFeatures(:,5:7);
             
        C=mean(nuccent); zgap=max(nuccent(:,3))-min(nuccent(:,3));
        
        
        s=strsplit(path,'Nuclei_and_Cells_');
        outputpath=strcat('Columnar_Structure_Prediction_8_100/',s{2});
       
      
        if (bonetype==3)|(bonetype==1)
        %zmin=-zmin;  zmax=-zmax;
        bonetypeMatrix=[1,0,0;0 1 0;0, 0,-1];
        nuccent1=(nuccent-C)*bonetypeMatrix;
        GlobCent=GlobalCenter{gi}{gj};
        nuccent2=(nuccent-GlobCent*bonetypeMatrix)*bonetypeMatrix;
        else
        %zmin=zmin+70;  zmax=zmax+70; 
        bonetypeMatrix=[1,0,0;0 1 0;0, 0,1];
        nuccent1=(nuccent-C)*bonetypeMatrix;
        GlobCent=GlobalCenter{gi}{gj};
        nuccent2=(nuccent-GlobCent*bonetypeMatrix)*bonetypeMatrix;
        end
        
        %[C, GlobalCenter{gi}{gj}]
        %[min(nuccent1(:,3)),  min(nuccent2(:,3)), max(nuccent1(:,3)),  max(nuccent2(:,3))]
        
        [CurvatureAxisLine,~]=movingMeanAverage(nuccent2);
        %dlmwrite(['Curvature/CurvLine_',s{2}(1:strlength(s{2})-1),'.dat'],CurvatureAxisLine,'\t');
        
         bone_curvature=CurvatureAxisLine;
        
        
        
       CC=plotClusterFigure(bone_curvature, pureColumn,LCCall,centroid,nuc, C, bonetypeMatrix,'r.');
       %save(['Curvature/Column_curvature_',s{2}(1:strlength(s{2})-1),'.mat'], 'CurvatureOfColumn');
       dlmwrite(['RadiusOfCurvature/ROC_',s{2}(1:strlength(s{2})-1),'.dat'], [CC.RadiusOfCurvature,CC.sphFit(:,[2]),CC.ColumnPos_PDaxis],'\t');

       
       
       RadOfCurv=CC.sphFit; t=find(~isnan(RadOfCurv(:,4)));       
       method1=corrcoef(RadOfCurv(t,1),RadOfCurv(t,2));
       method2=corrcoef(RadOfCurv(t,3),RadOfCurv(t,4));
       
       [gi,gj,method1(1,2),method2(1,2)]
      
       
       shp = alphaShape(nuccent1);
	   [tetrahedron,V] = alphaTriangulation(shp);
       trep = triangulation(tetrahedron, V);
       [tri, V] = freeBoundary(trep);
        
        plot3( V(:,1), V(:,2), V(:,3),'b.','markersize',0.01);
        hold on 
        
        plot3( CurvatureAxisLine(:,1), CurvatureAxisLine(:,2), CurvatureAxisLine(:,3),'r.-');
      
        title('Curvature axis line','fontweight','normal')
        axis image
        xlabel('X')
        ylabel('Y')
        zlabel('P-D axis')
        view(-38,30)
        
        savefolder='Curvature/figure_Curvature_column/';
        if ~exist([savefolder],'dir')
              mkdir([savefolder]);
        end
        
        saveas(h,[savefolder,'curvatureAxisLine_',s{2}(1:strlength(s{2})-1),'.png'])
        saveas(h,[savefolder,'curvatureAxisLine_',s{2}(1:strlength(s{2})-1),'.fig'])


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



function output=plotClusterFigure(bone_curv,pureColumn,LCCall,centroid,csurfVertices,globalCenter, bonetypeMatrix,mycolor)
  
    compareAlgorithm=[];
    RadiusOfCurvature=[];
    for k=1:length(LCCall)
        LCC=LCCall{k};
        if length(LCC)>=3
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

             
            curv_column1=findCurvatureUsingCentroid_method1(mu);
            curv_column2=findCurvatureUsingCentroid_method2(mu);  
            RadiusOfCurvature=[ RadiusOfCurvature;[k,length(LCC), curv_column1]];

             for i=1:length(pureColumn)
                 if k==pureColumn(i)
                        if curv_column1>5
                         plot3(mu(:,1),mu(:,2),mu(:,3),'k.-','linewidth',1);
                        else
                         plot3(mu(:,1),mu(:,2),mu(:,3),'m.-','linewidth',1);
                        end
                 end
             end

             
                Z=mu(:,3);      
                ColumnPos_PDaxis(k,:)=[min(Z),max(Z)];
                boneCurvIndex=[];
                for j=1:length(bone_curv)-1
                    if (bone_curv(j,3)<=min(Z) &  min(Z)<=bone_curv(j+1,3))
                        boneCurvIndex=[boneCurvIndex;(j:j+1)'];
                    end
                    if (bone_curv(j,3)<=max(Z) &  max(Z)<=bone_curv(j+1,3))
                        boneCurvIndex=[boneCurvIndex;(j:j+1)'];
                    end
                end
                   bone_curv_along_col=bone_curv(unique(boneCurvIndex),:);
               
             if  size(bone_curv_along_col,1)>2
                curv_bone1=findCurvatureUsingCentroid_method1(bone_curv_along_col);
                curv_bone2=findCurvatureUsingCentroid_method2(bone_curv_along_col);
             else
                 curv_bone1=0;   curv_bone2=0;
             end
             
             
             compareAlgorithm=[ compareAlgorithm;[ curv_column1,  curv_bone1,curv_column2,  curv_bone2]];       
             smoothline{k}=smoothUsingAverage(csurf,'g.-',0);
             columnCentroid{k}=mu;
             
             %smoothUsingEllipsoid(csurf,'b.-');
            
            
          
        end
    end
    
    
    
   output.CurvatureOfColumn=smoothline;
   output.NormalizedColumnCentroid=columnCentroid;       
   output.sphFit=compareAlgorithm;
   output.RadiusOfCurvature=RadiusOfCurvature;
   output.ColumnPos_PDaxis=ColumnPos_PDaxis;
end



function  Bone_curvature_def=findCurvatureUsingCentroid_method1(mu)
                       [rcol1,a,b,c]=sphereFit_Prem(mu); 

                       clear theta
                       for i=1:size(mu,1)
                           vec1=mu(i,:)-[a,b,c];
                           for j=i+1:size(mu,1)
                               vec2=mu(j,:)-[a,b,c];
                               theta(i,j)=angleCompute(vec1', vec2');
                           end
                       end

                       [sa,sb]=sort(theta(:),'descend');
                       %[size(mu,1), size(theta)]
                       [x,y] = ind2sub(size(theta),sb(1:2));
                        PlanePoints=mu(unique([x;y]),:);
                        P1=PlanePoints(1,:); P2=PlanePoints(2,:); P3=PlanePoints(3,:); 
                        points=mu(unique([x;y]),:)';  

    % %                    figure
    % %                    plot3(mu(:,1),mu(:,2),mu(:,3),'k.','linewidth',1);
    % %                    hold on 
    % %                    plot3(a,b,c,'bo','linewidth',1);
    % %                    plot3(mu([x;y],1),mu([x;y],2),mu([x;y],3),'ro','linewidth',1);
    % %                    
    % %                    [X,Y,Z] = sphere; radius=rcol1;
    % %                     X = X * radius; Y = Y * radius; Z = Z * radius;
    % %                     surf(X+a,Y+b,Z+c, 'edgecolor','none','facealpha',0.1)              
    % %                     fill3(points(1,:),points(2,:),points(3,:),'r'); alpha(0.2)

                        normal = cross(P1-P2, P1-P3);                     
                        d=point_to_plane_distance( [a,b,c], normal, P1 );

                        Bone_curvature_def= rcol1-d; %1/rcol1;

end



function  Bone_curvature_def=findCurvatureUsingCentroid_method2(mu)
                      warning 'off';
                      [center,rcol1]=sphereFit_Alan(mu); a=center(1);b=center(2);c=center(3); 
                       clear theta
                       for i=1:size(mu,1)
                           vec1=mu(i,:)-[a,b,c];
                           for j=i+1:size(mu,1)
                               vec2=mu(j,:)-[a,b,c];
                               theta(i,j)=angleCompute(vec1', vec2');
                           end
                       end
                       [sa,sb]=sort(theta(:),'descend');
                       [x,y] = ind2sub(size(theta),sb(1:2));
                        PlanePoints=mu(unique([x;y]),:);
                        P1=PlanePoints(1,:); P2=PlanePoints(2,:); P3=PlanePoints(3,:); 
                        points=mu(unique([x;y]),:)';  
                        normal = cross(P1-P2, P1-P3);                     
                        d=point_to_plane_distance( [a,b,c], normal, P1 );
                        Bone_curvature_def= rcol1-d; %1/rcol1;

end



 


function dd=smoothUsingAverage(surface,color,plotfig)

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
         
         if plotfig==1
         plot3(dd(:,1),dd(:,2),dd(:,3),color,'linewidth',5)
         end
end





function value=angleCompute(u,v) 
         value=atan2(norm(cross(u,v)),dot(u,v));
         value=180/pi*value;
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

