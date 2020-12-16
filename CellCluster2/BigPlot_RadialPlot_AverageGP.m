
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
        %GlobalCenter{gi}{gj}=meanOfAllCelltemp(count,[6,9]);
        count=count+1;
    end
end

GPname={'DT_WT','PT_WT','DT_MT','PT_MT','DU_WT'};
tname={'90-180','0-90','180-270','270-360'};


for gi=1:length(allpath)
    	bonetype=gi;
        clear boneData
        clear columnData
        clear Yall
        samples=1:length(allpath{gi});
	for gj=samples
		path=allpath{gi}{gj};
        s=strsplit(path,'Nuclei_and_Cells_');
        input2=strcat('Columnar_Structure_Prediction_8_100/',s{2});
          
        a1=load(['./../../',path,'all_cells_nuclei.mat']);
        nuc=a1.all_cells_nuclei; 
        nuccent=nuc(:,5:7);
        
        C=mean(nuccent); %zgap=max(nuccent(:,3))-min(nuccent(:,3));
        
        limit=GlobalCenter{gi}{gj};
        starting=limit(1)+ 0.2*(limit(2)-limit(1));
        finishing=limit(1)+ 0.6*(limit(2)-limit(1));
        [limit,starting,finishing];
        
        
        
        GlobCent=GlobalCenter{gi}{gj};
        
        
        fixzmin=starting; fixzmax=finishing;
        if (bonetype==3)|(bonetype==1)
        zmin=-fixzmin;  zmax=-fixzmax;
        colzmin=fixzmin;colzmax=fixzmax;
        bonetypeMatrix=[1,0,0;0 1 0;0, 0,-1];
        nuccent1=(nuccent-C)*bonetypeMatrix;
        %nuccent2=(nuccent-GlobCent*bonetypeMatrix)*bonetypeMatrix;
        else
        zmin=fixzmin;  zmax=fixzmax;  %70 add 
        colzmin=zmin;  colzmax=zmax;
        bonetypeMatrix=[1,0,0;0 1 0;0, 0,1];
        nuccent1=(nuccent-C)*bonetypeMatrix;
        %nuccent2=(nuccent-GlobCent*bonetypeMatrix)*bonetypeMatrix;
        end
        
       shp = alphaShape(nuccent1);
	   [tetrahedron,V] = alphaTriangulation(shp);
       trep = triangulation(tetrahedron, V);
       [tri, V] = freeBoundary(trep);
       

       %shift the centroids of the bones 
       shift=1;
       
       
        close all 
       
        [boneData,columnData]=plot2dColumns(V,s,shift);
      
         boneplot{gi}{gj}=  boneData;
         colplot{gi}{gj}=columnData;
         lengthofgraphletmotif=length(columnData)
         
         for fi=1:4 
              normalize=boneData.max_radialDistance;
              boneX=boneData.boneaxis;
              Yall{fi}(:,gj)=boneData.boundary_radialDistance{fi}/normalize;  
             
         end      
  
        
    end
    


     
     for fi=1:4
        temp(fi)=max(mean(Yall{fi},2));
     end

	%factor used to make avg growth plate radius 1. 
        factor=1/max(temp);


    
      if true
        h2=figure;
        
        set(gcf, 'PaperSize', [10 7]); %7
        set(gcf, 'PaperPosition', [0 0 10 7]);
        
       
       for fi=1:4
            subplot(2,2,fi)
            
            X=[];Y=[];
            for gj=samples
                boneData=boneplot{gi}{gj};
                columnData=colplot{gi}{gj};
               
               for j=1:length(columnData)
                    normalize=boneData.max_radialDistance;
                    index=find(columnData(j).boundary_radialDistance{fi});
                    
                    X=[X;columnData(j).boneaxis(index)];
                    Y=[Y;columnData(j).boundary_radialDistance{fi}(index)/normalize];
                    
%                     if columnData(j).length==3
%                         plot(columnData(j).boneaxis(index),columnData(j).boundary_radialDistance{fi}(index)/normalize,'mo-','markersize',1,'markerfacecolor','m')
%                     elseif columnData(j).length==4
%                         plot(columnData(j).boneaxis(index),columnData(j).boundary_radialDistance{fi}(index)/normalize,'ko-','markersize',3,'markerfacecolor','k')
%                     elseif columnData(j).length==5
%                         plot(columnData(j).boneaxis(index),columnData(j).boundary_radialDistance{fi}(index)/normalize,'go-','markersize',5,'markerfacecolor','g')
%                     else 
%                         plot(columnData(j).boneaxis(index),columnData(j).boundary_radialDistance{fi}(index)/normalize,'ro-','markersize',7,'markerfacecolor','r')
%                     end
%                     
%                     hold on 
%                     if j==1
%                         %d=[boneData.boneaxis, boneData.boundary_radialDistance{fi}/normalize];                    
%                         %plot(d(:,1),d(:,2),'b.-')
%                     end
               end
                
            end
            
            
              boneY=mean(Yall{fi},2);
              xy=gridAveragePlot(X,Y,boneX,factor*boneY);
              %plot(X,Y,'r.') hold on 
              scatter(xy(:,1),xy(:,2),5*xy(:,3),'markerfacecolor','b'); hold on 
           
              plot(boneX,factor*boneY,'k-'); % factor multiply to make one. 
        
              if fi>2
              xlabel('Bone P-D axis (0:RZ, 1:HZ)')
              end
              if ((fi==1)|(fi==3))
              ylabel('Bone radius (0:center, 1:Periphery)')
              end
              axis([0,1,0,1])
              title(tname{fi},'fontweight','normal')
              box on 
              hold off
        end
    
            directory='AvgGP';
            if ~exist([directory],'dir')
              mkdir([directory]);
            end
            saveas(h2,[directory,'/normalizedColumnsPos_',  GPname{gi},'.png'])
            
      end
    
end








function mydata=gridAveragePlot(x,y,bx,by)
    
        binsize=51;
        gx=linspace(0,1,binsize+1);
        gy=linspace(0,1,binsize+1);
        
        
        for i=1:length(gx)-1
            ngx(i)=mean(gx(i:i+1));
            ngy(i)=mean(gy(i:i+1));
        end
        
        
        
        index=0;
        for i=1:binsize
            for j=1:binsize
                index=index+1;
                mat(index,:)=0;
                newdata(index,:)=[ngx(i),ngy(j)];
            end
        end
                
         
        
        for i=1:binsize
           for j=1:binsize
               index=index+1;
               bonedata(i,j)=0;
                for k=1:length(bx)
                    if ((gx(i)<= bx(k))&(bx(k)< gx(i+1)) & (gy(j)<= by(k)) &(by(k)< gy(j+1)))
                        bonedata(i,j)=1;
                    end
                end
            end
        end
      
         
        for i=1:binsize
            if sum(bonedata(i,:))>0
                 temp=max(find(bonedata(i,:)));
                 boneboundary(i,1)=temp;
            else
                 boneboundary(i,1)=0;
            end
        end
              
        
        index=0;
        for i=1:binsize
           for j=1:binsize
               index=index+1;
                for k=1:length(x)
                    if ((gx(i)<= x(k))&(x(k)< gx(i+1)) & (gy(j)<= y(k)) &(y(k)< gy(j+1)))
                        %  [x(k),y(k), gx(i:i+1),gy(j:j+1)]
                        % Remove the centroid to outside of bone boundary 
                         if j<boneboundary(i)
                           mat(index,:)=mat(index,:)+1;
                         end
                    end
                end
            end
        end
        
        %mat=mat/ max(mat(:));
      
        index=1;
        for i=1:length(newdata)
            if mat(i)>0
                mydata(index,:)=[newdata(i,:),mat(i)];
                index=index+1;
            end
        end
        
        
end








function [boneData, columnData]=plot2dColumns(V,s,shift)
     

        no_of_bins=50;
        PDaxis_bone=linspace(0,1,no_of_bins);
        PDaxis1_bone=linspace(min(V(:,3)),max(V(:,3)),no_of_bins+1);
        
        no_of_bins=101;
        PDaxis_col=linspace(0,1,no_of_bins);
        PDaxis1_col=linspace(min(V(:,3)),max(V(:,3)),no_of_bins+1);
       
       % pureColumn=load(['degree_of_the_column/StraightLine/pureColum_',s{2}(1:strlength(s{2})-1),'.dat']);
        column_curvature=load(['Curvature/Column_curvature_',s{2}(1:strlength(s{2})-1),'.mat']) ;
        %ROC=load(['RadiusOfCurvature/ROC_',s{2}(1:strlength(s{2})-1),'.dat']);
        bone_curvature=load(['Curvature/CurvLine_',s{2}(1:strlength(s{2})-1),'.dat']);
        
        
        
        %read graphlet motifs 
        fid=fopen(['degree_of_the_column/StraightLine/Colum_GraphletID_',s{2}(1:strlength(s{2})-1),'.dat']);
        tline=fgetl(fid);
        count=1;
        while ischar(tline)          
                for j=1:length(tline)
                    gcolid{count}=tline;
                end
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
        
        pureColumn=[];
        for i=1:length(gcolid)
            if (strcmp(gcolid{i},'G1') | strcmp(gcolid{i},'G3') | strcmp(gcolid{i},'G9') | strcmp(gcolid{i},'N1') | strcmp(gcolid{i},'M1'))
            %if (strcmp(gcolid{i},'N1') | strcmp(gcolid{i},'M1'))
                pureColumn=[pureColumn;i];
            end
        end

       % [length(tpureColumn), length(pureColumn)];
      

        %USE THIS COMMAND FOR SHIFTING THE CENTROIDS 
        if shift==1
            nV=centroid_shift_of_bone(bone_curvature,V);
        else
            nV=V;
        end
      
        [boneaxis,boundary_radialDistance, max_radialDistance,mainX,mainY]=find_external_boundaries(PDaxis_bone,PDaxis1_bone,nV); 
        %radialDistance=  sqrt(V(:,1).^2 + V(:,2).^2);
        
        boneData.boneaxis=boneaxis;
        boneData.boundary_radialDistance=boundary_radialDistance;
        boneData.X=mainX;
        boneData.Y=mainY;
        boneData.max_radialDistance=max_radialDistance;      

        
%        figure    
%         for i=1:4
%             subplot(2,2,i)
%             %k=boundary(boneaxis(ind{i}),radialDistance(ind{i}));
%             X=boneaxis;
%             Y=boundary_radialDistance{i}/max_radialDistance;
%             plot(X,Y,'b.-')
%         end
        
%             plot3(nV(:,1),nV(:,2),nV(:,3),'b.')
%             hold on 
%             plot3(bone_curvature(:,1),bone_curvature(:,2),bone_curvature(:,3),'r.')

      
            for j=1:length(pureColumn)
                cent=column_curvature.CurvatureOfColumn.NormalizedColumnCentroid{pureColumn(j)};         
                %plot3(cent(:,1),cent(:,2),cent(:,3),'o')
                %USE THIS COMMAND FOR SHIFTING THE CENTROIDS 
                if shift==1
                        ncent=centroid_shift_of_bone(bone_curvature,cent);
                else
                        ncent=cent;
                end
                [boneaxis,boundary_radialDistance, max_radialDistance,X,Y]=find_external_boundaries(PDaxis_col,PDaxis1_col,ncent); 
                columnData(j).boneaxis=boneaxis;
                columnData(j).boundary_radialDistance=boundary_radialDistance;
                columnData(j).X=X;
                columnData(j).Y=Y;
                columnData(j).length=size(cent,1);
            end
        
           
            
%           h1=figure;
%           for fi=1:4
%                subplot(2,2,fi)
%                for j=1:length(columnData)
%                     plot(columnData(j).X{fi},columnData(j).Y{fi},'ro-','markersize',5)
%                     hold on 
%                     if j==1
%                     plot( boneData.X{fi}, boneData.Y{fi},'b.')
%                     end
%                end
%                 title(tname{fi},'fontweight','normal')
%           end
            
 
          
        
        
         
         
end




function nv=centroid_shift_of_bone(bone_curvature,V)
	x=bone_curvature(:,3);
	nV=zeros(size(V)); 
	for i=1:size(V,1)
        for j=1:length(x)
             dist(j)=abs(V(i,3)-x(j));
        end
        [sa,sb]=min(dist);
        localcenter=bone_curvature(sb,[1:2]);    
        nv(i,[1,2])=V(i,[1,2])-localcenter;
        nv(i,3)=V(i,3);
    end
end 


function [Xaxis,Yaxis,rmax,X,Y]=find_external_boundaries(PDaxis,PDaxis1,V)
    	%PD axis  =  # of bins 
        %PD axis 1 = # of bins +1 
        %PD axis 2 = # of bins for bone centroid shift cutoff
       
	 
        [TH,R] = cart2pol(V(:,1),V(:,2)); TH=TH*180/pi;
        ind{2}= find((TH>=0) & (TH<90));
        ind{1}= find((TH>=90) & (TH<180));
        ind{4}= find((TH<0) & (TH>=-90));
        ind{3}= find((TH<-90) & (TH>=-180));

        radialDistance=  sqrt(V(:,1).^2 + V(:,2).^2);
        rmax=max(radialDistance);

         for j=1:length(V)
            for i=1:length(PDaxis1)-1
                    if ((PDaxis1(i)<V(j,3)) & (PDaxis1(i+1)>=V(j,3)))
                    boneaxis(j,1)=PDaxis(i);
                    end
            end             
         end
     
        for i=1:4
            %subplot(2,2,i)
            %k=boundary(boneaxis(ind{i}),radialDistance(ind{i}));
            X{i,1}=boneaxis(ind{i});
            Y{i,1}=radialDistance(ind{i});
            %plot(X,Y,'b.')
        end
    

	for j=1:4
		data=V(ind{j},3);
		RD=radialDistance(ind{j});
		Yaxis{j}=zeros(length(PDaxis),1);
		for i=1:length(PDaxis1)-1
		        index=find((PDaxis1(i)<data) & (PDaxis1(i+1)>=data));
		        if length(index)>0
		        Yaxis{j}(i,1)=max(RD(index));        
		        end
	    	end

    end
	Xaxis=PDaxis';
end


% 
% 
% function nv=centroid_shift_of_bone2(PDaxis2,bone_curvature,V)
% 	x=bone_curvature(:,3);
% 	nV=zeros(size(V)); 
% 	for i=1:length(PDaxis2)-1
% 		index=find((PDaxis2(i)<x) & (PDaxis2(i+1) >=x));
% 		localcenter=mean(bone_curvature(index,[1:2]));    
% 		index=find((PDaxis2(i)<=V(:,3)) & (PDaxis2(i+1)>=V(:,3)));
%         	nv(index,[1,2])=V(index,[1,2])-localcenter;
% 		nv(index,3)=V(index,3);
%     end
% end 




   
        



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
    
