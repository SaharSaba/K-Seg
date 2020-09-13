
% Before everything, I want to say there is no need to crop the image and
% divide it into 2 parts(cervical and torasic). The only thing to do is
% just press the Run buttom and give an input image to the algoritm and then get an 
% output image as" Spinalcord.nii". 
% But in some cases the algorithm couldnot reach to the ideal output
%(espeially in large size of images) , so here you can crop the image using 
% "FSL toolbox" and then give it to the algoritm.  


clear all; close all; clc;

%% Extracting Region of interest in the left–right direction

tic

% Getting nifti images from local computer
[MainImg,PathName,FilterIndex] = uigetfile({'*.nii';'*.nii.gz'});
imgM=niftiread(strcat(PathName,MainImg));


% Keep all the metadata information( size, dimentions)
info = niftiinfo(strcat(PathName,MainImg));

[a , b, r] = size(imgM);

if (a<b)
    
    % Rotate each image to the saggital view to find the symmetry slice using MI
    
    img = permute( imgM, [3 2 1 ] );
    img = imrotate(img ,90 );
    
else
     img = imgM;
    
end


MatSim= zeros(r,r);
[a , b, r] = size(img);
k=1;
len=10;

% Extract the middle slice of image by mutual information(MI) method
for dim =1 :r
    
    if (dim ~=r)
        
        ref = img(:,:,dim);
        im = img(:,:,r);
        r= r-1;
        CR  = MI(ref,im);
        MatSim(dim,k)=CR;
        k=k+1;
        
    end
    
end

% The maximum value of MI shows the middle sagittal view.
% The symmetry of the body helps to perform next phases better.
maxx = max(max(MatSim));
[simetry ,y ] = find(MatSim==maxx);


% After we find the symmetry body slice, we need to add and minus middle
% slice number with "len=10" inorder to expand the ROI. 

roi = img(:,: ,simetry-len:simetry+len );
niftiwrite(roi,'ROI.nii' );


%%


img = niftiread('ROI.nii');
[a11 b11 r] = size(img);
LineMat= zeros(size(img));


% applying a canny filter to extract the 3D edges and Hough line transform
% to remove the extra parts from the ROI in the Anterior-Posterior direction
BW1 =  edge3 (img,'approxcanny',0.1 );


for dim= 1:r
    
    
    re = (img(:,:,dim));
    BW= BW1(:,:,dim);
    [H,T,R] = hough(BW);
    
    % Finding all the lines with length=50
    P  = houghpeaks(H,50);
    
    lines = houghlines(BW,T,R,P );
    [x ,y ] = size(lines);
    L= bwlabel(BW);    
    [m ,n] = size(L);    
    h=1;
    for i = 1 :y      
        
    % Connecting starting and ending points to form vertical lines alongside
    %of the spinal cord in saggital view.
        
        xy = [lines(i).point1; lines(i).point2];
        x11= xy(1,1);
        y11= xy(1,2);
        x22= xy(2,1);
        y22= xy(2,2);
        
        for k=1 :m
            for j=1 :n
                
                lab= L(k,j);
                if(  (j==x11   && k==y11  ) || (j==x22   &&  k==y22) )
                    
                    ind(1,h) = lab;
                    h=h+1;
                    
                end
            end
        end
    end
    
    [R T] = size(ind);
    [x,y] = size(L);
    L2 = zeros(size(L));
    
    
    % Removing the vertical lines which are out of the ROI and maintaining
    % those regions that are surrounding of the spinal cord.
    for i = 1 :x
        for j = 1 :y
            for k = 1 :T
                
                if (L(i,j)== ind (1,k))
                    
                    L2(i,j)=1;
                    
                end
            end
        end
    end
    
    
    % Combining each slice with the result of the previous step to find maximum
    % average that shows the brightest regions of CSf .
    C = imfuse(re , L2 ,'blend');
    
    avg = mean(mean(C));
    [ro co]=size(C);
    
    for i=1 : ro
        
        for j=1: co
            
            if (C(i,j)>=avg)
                
                li(i,j)=1;
                
            else
                
                li(i,j)=0;
                               
            end
                        
        end
    end
    
    LineMat (:,:,dim)= li;
    
    
end

niftiwrite(LineMat,'LineMat.nii');

%% Spinal Cord Detection by Circular Hough Transform

img = niftiread('LineMat.nii');
imgT2 = niftiread('ROI.nii');

% Rotate each image to the axial view in order to find the circular shape
img = permute( img, [ 3 1 2 ] );
imgT2 = permute( imgT2, [ 3 1 2 ] );
[a11 b11 r] = size(img);

% Extracting the center of the axial images
Xcenter = round(a11/2);
Ycenter = round(b11/2);
ImgCenMain = [Ycenter Xcenter];
mask2=zeros(size(img));
distance =zeros([r 1]);
mask3=zeros(size(img));

for dim=1:r
    
    
    % Empty all the variables in each repetition of the loop.
    dis=[];
    
    
    % Finding all circles in each axial slice with radious= 7-9 because those
    % circles are similar to the cord shape.
    [centers,radii] = imfindcircles(img(:,:,dim),[7 9],'Sensitivity',0.97);
    
    num = numel(centers);    
    num = num/2;
        
    if(num >1 && num ~=0 )
        
        
        %         figure, imshow( imgT2(:,:,dim),[]);
        %
        %         viscircles(centers,radii,'EdgeColor','r');
        
        
        for j=1: num
        
            % Calculating the distance of all the circles in each slice to
            % the center and saving thier radious.
            dis(j,1) = norm ( ImgCenMain - centers(j,:) );
            
            dis(j,2) = radii(j,:);
            
            
        end
        
      % Discovering a suitable and nearset circle to the center of the
      % image as start point

        MinDis = min(dis(:,1));
        
        [xx yy]= find( dis(:,:) ==MinDis);
        
        distance(dim ,1) = MinDis;
        
        AllCenters (dim,1) = centers(xx,1);
        
        AllCenters (dim,2) = centers(xx,2);
        
        AllRadii(dim,1) = radii(xx,1);
        
        Slice(dim,1)= dim;
        
        
    elseif (num==1)
        
        %         figure, imshow( imgT2(:,:,dim),[]);
        %
        %         viscircles(centers,radii,'EdgeColor','r');
        
        distance(dim ,1) = norm ( ImgCenMain - centers(1,:) );
        
        AllCenters (dim,1) = centers(1,1);
        
        AllCenters (dim,2) = centers(1,2);
        
        AllRadii(dim,1) = radii;
        
        Slice(dim,1)= dim;
        
        
    end
    
end

distance(distance==0)= NaN;

X = distance;

res = ~any(~isnan(X(:)));


if ( res ~=1 )
    
     % Discovering a suitable and nearset circle to the center of the
     % image as start point
    MinDis = min(distance(:,1));
    
    [xx yy]= find( distance(:,1) ==MinDis);
    
    ImgCen =  AllCenters(xx, :);
    
    radi = AllRadii(xx,:);
    
    StartPoint =  Slice(xx,1);
    
    ImgCen1 = ImgCen;
    
%     
%     
%     figure, imshow( imgT2(:,:,StartPoint),[]);
%     
%     viscircles(ImgCen,radi,'EdgeColor','r');
%     
    
    % Here we can start from the slice that include start point circle to the
    % end of the slices.
    
    % Usualy thoracic regions of the spine (e.g., the middle slices of the
    % selected volumes) yielded better results because of the higher contrast
    % and shape of the spinal cord. So the startpoint set to the middle part of
    % the image. Here we need to divide image into 2 parts to compare detected
    % circles with startpoint.
    
    
    
    for dim=StartPoint:r
        
        
%       figure, imshow( imgT2(:,:,dim),[]);
        
        centers = [];
        radii = [];      
        dis1= [];
        
        
        % In each slice, it is better to find several circles with radious=6-12
        % due to the diversity of cord shape alongside of the spinal cord region.(
        % cervical, thorasic, lumbar)
        
        [centers,radii] = imfindcircles(img(:,:,dim),[6 9],'Sensitivity',0.99);
        num = numel(centers);
        
        if(num >0)            
            
            % In each slice we have access to several circles, so we need to find the
            % closest circle to startpoint circle coordinate.
            
            for i= 1: num/2
                dis1(i,1) = norm ( ImgCen - centers(i,:) );
            end
                        
            MinDis = min(dis1(:,1));
            
            [xMinDis yy]= find( dis1(:,1) ==MinDis);            
            
            cen =  centers(xMinDis, :); 
                
            % Sometimes finding circle do not show the surrounding of the area and
            % it causes to miss some parts of the cord, so we can increase the radious
            % of the discovered circle.
            radi = radii(xMinDis,:)+2;
            
            ImgCen = cen;
            
         % viscircles(cen,radi,'EdgeColor','r');
            
        % Creat a mask around the finding circle and keep it in the new image .
            mask = createCirclesMask(img(:,:,dim) ,cen,radi);            
                 
            
        % After finding the optimal circle in each slice, we need to draw a mask
       % around it with the intensities of main image. 
            
            for k=1 :a11
                for l=1 :b11
                    if(   mask(k,l)==1 )
                        
                        mask3(k,l,dim)= imgT2(k,l,dim);
                    end
                end
            end
            
            
        end
    end
    
    
    % We repeat the pervious step code here for downside of the image.
    for dim=StartPoint:-1:1       
        
                
%       figure, imshow( imgT2(:,:,dim),[]);
        
        centers = [];
        radii = [];
        dis1= [];
        [centers,radii,metric] = imfindcircles(img(:,:,dim),[3 12],'Sensitivity',0.99);
        num = numel(centers);        
        
        if(num >0)
            
            for i= 1: num/2
                dis1(i,1) = norm ( ImgCen1 - centers(i,:) );
            end      
            
            
            MinDis = min(dis1(:,1));
            [xMinDis yy]= find( dis1(:,1) == MinDis);
            
            cen =  centers(xMinDis, :);
                
            radi = radii(xMinDis,:)+2;
                        
            ImgCen1 = cen;
            
%           viscircles(cen,radi,'EdgeColor','r');
            
            mask = createCirclesMask(img(:,:,dim) ,cen,radi);
                                  
            for k=1 :a11
                for l=1 :b11
                    if(  mask(k,l)==1 )
                        
                        mask3(k,l,dim)= imgT2(k,l,dim);
                    end
                end
            end                                  
            
        end       
    end  
    
    % Finaly, we achive to  all detected circles around the spinal cord
    % with original intensities.
    niftiwrite(mask3,'Circle.nii');
    
    
    %% Anisotropic Diffusion Filter  AND  k-means Algorithm
    
    % The Anisotropic Diffusion Filter (AD) allows us to
    % combine the two most important attributes of the denoising algorithm: edge
    % preservation and noise removal.
    
 img = niftiread('Circle.nii');
[a b r] = size(img);
SpinalCord= zeros(size(img));


num_iter = 4;

% integration constant
delta_t = 1/44; 

% gradient modulus threshold that controls the conduction
kappa = 20;

% conduction coefficient functions
option = 1; 

% 3x1 vector column with the x, y and z dimensions
voxel_spacing = ones(3,1); 
img = anisodiff3D(img, num_iter, delta_t, kappa, option, voxel_spacing);

 % Median filter to smooth image better. 
img = medfilt3(img);

% K-means algorithm is sensitive to noise, so after the filtering image,
% we can cluster the image properly.

for dim =1 : r
    int =( img(:,:,dim));
    disk = int(:);
    disk(disk==0)= NaN;
    X = disk;
    res = ~any(~isnan(X(:)));
    if ( res ~=1 )
        
    % The number of clusters that we consider to segment spinal cord and canal
    % areas is (k=3). The first cluster is related to the spinal cord region,
    % the second class includes the spinal canal area and the third cluster is
    % referred to keep the same false areas which added on both SC and CSF
    % throughout the detecting circles. By switching among these clusters the
    % automatic SC and CSF segmentation can be achieved as well.
        
        
        
        [idx,C] = kmeans(X,3,'Replicates',10);
        idxx = reshape (idx , [a b]);
        minn = min(C(:,1));
        maxx = max(C(:,1));
        
        % This part can determine the cord and canal and also is realy
        % important to specify in T1 and T2 image by swiching between max and min value.
        % In T1 images, the cord has highest intentensity so we need to choose "maxx" value.
        %but in T2 it is important to change this line, so this line is the most important 
        %difference between T1 & T2. 
      
        [Xcsf yy]= find( C(:,1)== maxx );
        for i=1 :a
            for j=1 :b
                if (idxx(i,j) == Xcsf )
                    SC(i,j)= img(i,j,dim);
                else
                    SC(i,j )=0;
                end
            end
        end
        SpinalCord (:,:,dim) = SC ;
    end
end

niftiwrite(SpinalCord,'K-SpinalCord.nii');
    
       
%% 
    
% After discovering the 2 regions( spinal cord and canal using k-means),
% we can still see some unwanted parts in images, to remove these parts we
% use “regionprops” method in sagittal and axial views.

    
    SC = niftiread('K-SpinalCord.nii');
    
    SC = permute( SC, [ 2 3 1] );
    [a11 b11 r] = size(SC);
    SCSag  = zeros(size(SC));
    
    for dim=1:r
        
        BW = SC(:,:,dim);
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
        SCSag(:,:,dim)= BW;             
        
    end
    
    SCSag = permute( SCSag, [3 1 2] );
    niftiwrite(SCSag,'SCSagittal.nii');
    
    %%
    %Remove extra points AND cleaning in Axial view
    
    SC = niftiread('SCSagittal.nii');
    [a11 b11 r] = size(SC);
    SCax  = zeros(size(SC));
    
    for dim=1:r
        BW = SC(:,:,dim);
        L = bwlabel(BW);
        CC = bwconncomp(BW);
        bw_label= labelmatrix(CC);
        rp = regionprops(bw_label,'FilledArea','BoundingBox','Centroid','PixelIdxList');
            
  % Extracting the maximum part of the image that related to length of cord and canal.
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
        
        % Removing extra parts from the array.
        arr = find(array==0);
        indexesOfMinregion1 = [rp.FilledArea]~= array;
        pixelsNotToShow1 = vertcat(rp(indexesOfMinregion1).PixelIdxList);
        BW(pixelsNotToShow1) = 0;
        SCax(:,:,dim)= BW;
        
    end
    
    % Calculate of the area of this mask to estimate of cross section.
    % here we have access to the extracted pure spinal cord and in this
    % step we can calculate spinal cord cross sectional area.
    for i=1 : r
        
        im = (SCax(:,:,i));
        total = bwarea(im);
        Cross (i,1)= total;
        
        
    end
    
    niftiwrite(SCax,'SCAxial.nii');
    
    %% 
    
%Convert detected image as segmentation spinal cord and canal to
% original size of image 


    [aa , bb, rr] = size(imgM);
    ImgSeg1 = zeros (size(imgM ));    
    [a b r] = size(SCax);
    
    LT = simetry -len;
    RT = simetry +len;
    
    if (rr<r)
        
        imgT2 = permute( SCax, [ 2 3 1] );
        [i j k] = size(imgT2);
        ImgSeg1 ( : , : , LT:RT ) = imgT2(: , :, 1:k );
        
        
    else
        
        [a b r] = size(SCax);
        ImgSeg1 (LT:RT , : , :) = SCax(1:a , :, :);
                        
    end
    
    type = info.Datatype;
    comp = strcmp(type,'single');
    
    if (comp ==0 )
        
        ImgSeg = int16(ImgSeg1);
        
    else
        
        ImgSeg = single(ImgSeg1);
        
    end
        
    % It is the Main output as Sinal cord segmentation. 
    niftiwrite(ImgSeg,'SpinalCord.nii',info);
    
    toc
    
else
    % If the algorithm could not find an appropriate circle as startpoint 
    % we need to expand the ROI space with increasing "len=12" and then 
    % repeat all the above steps again. 
   
    wideT1(imgM, info);
    
end




