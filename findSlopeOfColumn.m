function  output=findSlopeOfColumn(bone_curv,col_curv_structure,truecolumn,cutoff)
%function  output=findSlopeOfColumn(gcolid,LCCIndex,doubletIndex,curvatureLine,colzmin,colzmax,cutoff)

      col_curv=col_curv_structure.CurvatureOfColumn;
      col_cent=col_curv_structure.NormalizedColumnCentroid;

      %[size(col_curv), size(truecolumn)]
      
      twocell=[];
      sampleSlope=[];
 	 count=1;
	 %for k=1:length(col_curv)
      for k1=1:length(truecolumn)
           k=truecolumn(k1);
           %colid=col_curv{k};  % use this for column curvature 
           colid=col_cent{k}; % use this for centroid of column 
           if size(colid,1)>=cutoff
                coordinate=colid;   
              
                
                Z=coordinate(:,3);   
                
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
                
                
               if size( bone_curv_along_col,1)>=3
                     
%               [min(curvatureLine(:,3)),max(curvatureLine(:,3)), min(boneCurv(:,3)),max(boneCurv(:,3))]
%               [min(Z),max(Z)]
                %[size(boneCurv),size(coordinate)]
                PCbon=pca(bone_curv_along_col);
                PCcol=pca(coordinate);
                        
                angle_diff=angleCompute(PCbon(:,1),PCcol(:,1));
                
                
                %[colzmin,colzmax];
                colzmin=-1000;
                colzmax=1000;
              
                    Zindex=find(Z>colzmin & Z<colzmax);
                    if length(Zindex)~=0
                          twocell{count,1}=coordinate;
                          %sampleSlope(count,1)=allSlope(k,1);
                          sampleSlope(count,1)=angle_diff;
                          count=count+1;
                    end

               end
         
           end
           
          
     end 
     
     output.twocell=twocell;
     output.sampleSlope=sampleSlope;
     
end


function value=angleCompute(u,v) 
         value=atan2(norm(cross(u,v)),dot(u,v));
         value=180/pi*value;
         if value>90
             value=180-value;
         end
end