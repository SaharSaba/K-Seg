







function BW2 = step()



clear all; close all; clc;


imgc = niftiread('Circle3.nii');
   
MiddImg = permute( imgc, [2 3 1 ] );
[a11 b11 r] = size(MiddImg);


% MiddImg1= MiddImg(:,:,len);
BW2=zeros(size(MiddImg));
midd=zeros(size(MiddImg));

for dim= 1:r
   
    BWW(:,:,dim) = fibermetric(MiddImg(:,:,dim) ,1,'ObjectPolarity', 'dark', 'StructureSensitivity',10);
   
end
%% 
                   
%        figure, imshow( BWW,[]);


   
        for dim=1:r
            
            BW = BWW(:,:,dim);
            L = bwlabel(BW);
            CC = bwconncomp(BW);
            bw_label= labelmatrix(CC);
            rp = regionprops(bw_label,'FilledArea','BoundingBox','Centroid','PixelIdxList');
            [maxx , inde]= max([rp.FilledArea]);
            [m,n] = size(rp);
            array= zeros([1,m]);
            
            for i=1 : m
                if ((rp(i).FilledArea)< maxx )
                    
                    array(1,i)=0;
                else
                    array(1,i)= rp(i).FilledArea;
                    
                end
            end
            
            arr = find(array==0);
            indexesOfMinregion1 = [rp.FilledArea]~= array;
            pixelsNotToShow1 = vertcat(rp(indexesOfMinregion1).PixelIdxList);
            BW(pixelsNotToShow1) = 0;
                       
            mid(:,:,dim)= BW;
            
           
     
        end 
        
%         figure, imshow( mid,[]);
%       niftiwrite(mid,'Midd.nii');

      
      
for dim= 1:r    
            
 midd(:,:,dim) = fibermetric(mid(:,:,dim) ,5,'ObjectPolarity', 'dark', 'StructureSensitivity',10);

 
end



%         figure, imshow( midd,[]);
        
        
        
        

for  dim=1 :r 

    for i=1 : a11
         for j=1 : b11
        
                if (midd(i,j,dim)>0)
                    
                    BW2(i,j,dim) = MiddImg(i,j,dim);
                    
                end
            
         end
         
    end
    
end 


%        figure, imshow( BW2(:,:,8),[]);


BW2 = imfill(BW2,18,'holes') ;


BW2 = permute( BW2, [3 1 2] );

    
   
niftiwrite(BW2,'Midd.nii');






end

