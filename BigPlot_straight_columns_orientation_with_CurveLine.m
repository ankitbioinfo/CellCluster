
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

       
allpath={dt_path_wt; pt_path_wt; dt_path_mut; pt_path_mut};  



for gi=1:length(allpath)
    bonetype=gi;
	for gj=1:length(allpath{gi})
		path=allpath{gi}{gj};
        s=strsplit(path,'Nuclei_and_Cells_');
        %ns=strsplit(path,'T_');
        %input2=strcat('Columnar_Structure_Prediction_8_100/',s{2});
        %input1=strcat('MakeListColumnarStructurePrediction/',s{2});   

      
    
        column_curvature=load(['Curvature/Column_curvature_',s{2}(1:strlength(s{2})-1),'.mat']) ;
        bone_curvature=load(['Curvature/CurvLine_',s{2}(1:strlength(s{2})-1),'.dat']);
        pureColumn=load(['degree_of_the_column/StraightLine/pureColum_',s{2}(1:strlength(s{2})-1),'.dat']);
      
        
        a1=load(['./../../',path,'all_cells_nuclei.mat']);
        nuc=a1.all_cells_nuclei; 
        nuccent=nuc(:,5:7);
        
        C=mean(nuccent); zgap=max(nuccent(:,3))-min(nuccent(:,3));
        nuccent1=nuccent-C;
        fixzmin=50; fixzmax=250;
        if (bonetype==3)|(bonetype==1)
        zmin=-fixzmin;  zmax=-fixzmax;
        colzmin=fixzmin;colzmax=fixzmax;
        else
        zmin=fixzmin+70;  zmax=fixzmax+70;
        colzmin=zmin;  colzmax=zmax;
        end
        
        shp = alphaShape(nuccent1);
	   [tetrahedron,V] = alphaTriangulation(shp);
       trep = triangulation(tetrahedron, V);
       [tri, V] = freeBoundary(trep);
       
        if zmin<zmax
        PZboundary=find(V(:,3)>zmin & V(:,3)<zmax);
        nuc_PZ_boundary=find(nuccent1(:,3)>zmin & nuccent1(:,3)<zmax);
        else
        PZboundary=find(V(:,3)<zmin & V(:,3)>zmax); 
        nuc_PZ_boundary=find(nuccent1(:,3)<zmin & nuccent1(:,3)>zmax);
        end
        
        x=V(PZboundary,1); y=V(PZboundary,2);
        K=boundary(x,y);
        databoundary{gi}{gj}=[x(K),y(K)];
        
        
                   
        for cutoff=[3,4,5]
            output=findSlopeOfColumn(bone_curvature, column_curvature.CurvatureOfColumn,pureColumn,cutoff);
            dataslope{gi}{gj}{cutoff}=output.sampleSlope;
            data{gi}{gj}{cutoff}=output.twocell;
        end
            

         nucleiDensity{gi}{gj}=nuccent1( nuc_PZ_boundary,1:2);
     
        
    end
end





name{2}={ 'S18 m6','S17 m2', 'S84 m3', 'S51 m2', 'S84 m4'  };
name{1}={ 'S18 m6','S17 m2', 'S84 m3', 'S51 m2', 'S84 m4' };  
    
name{3}={ 'S17 m1', 'S18 m2', 'S84 m1', 'S84 m5' }; 
name{4}={ 'S17 m1', 'S18 m2', 'S84 m1', 'S84 m5' }; 
name{5}={ 'S51 m2', 'S84 m2', 'S84 m3'};

tname{5}='DU WT'; tname{4}='PT MT'; tname{3}='DT MT'; tname{2}='PT WT'; tname{1}='DT WT';



mycolor={'','','m','k','g','r','y'};

h1=figure();
XL=0.02;XR=0.03;XGap=0.02;Row=4;
YT=0.06;YB=0.03;YGap=0.07;Col=5;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [19.2 14.4]);
set(gcf, 'PaperPosition', [0 0 19.2 14.4]);


U=[1,0,0];      mymap=colormap(redbluecmap(6));
 map=(colormap(jet(6)));
 
 


for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        if j<=length(data{i})
                subplot('Position',marray);
             
                pzboundary=databoundary{i}{j};
                nuclei=nucleiDensity{i}{j};
                
                for cutoff=3
                
                slope=dataslope{i}{j}{cutoff};
                column_coordinate=data{i}{j}{cutoff};
                
                mudata=[];
                for k=1:size(column_coordinate,1)
                       PC=pca(column_coordinate{k});
                       GC=mean(column_coordinate{k});
                       V=column_coordinate{k};
                       mudata(k,:)=GC(1:2);
                end
                
                C=mean(nuclei);
                
                pzboundary=pzboundary-C;
                bmin=min(pzboundary);
                bmax=max(pzboundary);
                
                [thb,rb]=cart2pol(pzboundary(:,1),pzboundary(:,2));
                
                
                if size(mudata,1)>1
                    mudata=mudata-C;
                    correctIndex=find( (mudata(:,1)>bmin(1))&(mudata(:,1)<bmax(1))&(mudata(:,2)>bmin(2))&(mudata(:,2)<bmax(2)));
                    [th,r]=cart2pol(mudata(correctIndex,1),mudata(correctIndex,2));
                    slope=slope(correctIndex,1);
                    clear myslope 
                    myslope=[];
                    for p=1:length(slope)
                        myslope(p)=slope(p);
                        if slope(p)>90
                            myslope(p)=180-slope(p);
                        end
                    end
                    if length(myslope)>0
                    [int,orient]=avgOrientation(th,myslope);
                    statistical_test{i}{j}=[int(1:end-1)', orient(1:end-1)'];
                    polarplot(int,orient,'color',mycolor{cutoff},'marker','.','linestyle','-');
                    end
             
                hold on 
                end
                
          
           
                set(gca,'fontsize',7);
            
                
            
                
                end
                
          
            title(strcat(tname{i}, '-',name{i}{j}),'fontweight','normal','fontsize',7)    
                         
        end
        
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end


saveas(h1,['column_cs_mean_orient_curv_3.png']);

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
    

function [int,y]=avgOrientation(th,myslope)
    interval=linspace(-pi-0.001,pi+0.001,13);
    %interval=[interval,interval(end)+interval(2)-interval(1)];
    for j=1:length(interval)-1
        x=[];
        for i=1:length(th)
            if ((th(i)>interval(j)) & (th(i) <interval(j+1)))
                x=[x,myslope(i)];
            end
        end
        if length(x)==0
            y(j)=0;
        else
             y(j)=mean(x);
        end
        int(j)=mean([interval(j),interval(j+1)]);
    end
    int(j+1)=int(1);
    y(j+1)=y(1);
end
    



