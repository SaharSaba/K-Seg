



clear all; close all; clc;

%% Extracting Region of interest in the left–right direction

tic
[MainImg,PathName,FilterIndex] = uigetfile({'*.nii';'*.nii.gz'});
imgM=niftiread(strcat(PathName,MainImg));

% info = niftiinfo(strcat(PathName,MainImg));


% files = gunzip('*.gz');

[a , b, r] = size(imgM);



if (a<b)
    
    
    img = permute( imgM, [3 2 1 ] );
    
    img = imrotate(img ,90 );
    
else
    
    img = imgM;
    
end


MatSim= zeros(r,r);
[a , b, r] = size(img);
k=1;
len=8;

for dim =1 :r
    
    if (dim ~=r)
        
        ref = img(:,:,dim);
        im = img(:,:,r);
        r= r-1;
        CR  = mi(ref,im);
        MatSim(dim,k)=CR;
        k=k+1;
        
    end
    
end

maxx = max(max(MatSim));
[simetry ,y ] = find(MatSim==maxx);
roi = img(:,: ,simetry-len:simetry+len );
niftiwrite(roi,'ROI.nii' );


%%
%Removing unwanted regions inside of the ROI in the Anterior–Posterior direction


img = niftiread('ROI.nii');
[a11 b11 r] = size(img);

LineMat= zeros(size(img));


for dim= 1:r
    
    
    BW(:,:,dim) = fibermetric(img(:,:,dim) ,15, 'StructureSensitivity',100);
    
    
end


LineMat= im2double(BW);
niftiwrite(LineMat,'LineMat.nii');




%% Spinal Cord Detection by Circular Hough Transform

img = niftiread('LineMat.nii');
imgT2 = niftiread('ROI.nii');
img = permute( img, [ 3 1 2 ] );
imgT2 = permute( imgT2, [ 3 1 2 ] );
[a11 b11 r] = size(img);
Xcenter = round(a11/2);
Ycenter = round(b11/2);
ImgCenMain = [Ycenter Xcenter];
mask2=zeros(size(img));
% ImgCen=zeros(size(img));

mask3=zeros(size(img));


for dim=1:r
    
    dis=[];
    [centers,radii] = imfindcircles(img(:,:,dim),[7 9],'Sensitivity',0.97);
    num = numel(centers);
    
    for i=1 : num/2
        
        %
        %         figure, imshow( imgT2(:,:,dim),[]);
        %
        %         viscircles(centers(i,:),radii(i,:),'EdgeColor','r');
        %
        %
        mask = createCirclesMask(img(:,:,dim) ,centers(i,:),radii(i,:));
        
        
        for k=1 :a11
            
            for l=1 :b11
                
                if(mask(k,l)==1)
                    
                    mask2(k,l,dim) = imgT2(k,l,dim);
                    
                end
                
            end
        end
        
        
    end
    
    
end


for i=1 : r
    
    im = (mask2(:,:,i));
    total = bwarea(im);
    CrossSec(i,1)= total;
    
    
end




for i=1 : r
    
    if (CrossSec(i,1)>50 )
        
        
        CrossSection(i,1)= CrossSec(i,1);
    else
        
        CrossSection(i,1)=0;
    end
    
end


CrossSection(CrossSection==0)= NaN ;

AvgCross = nanmean (CrossSection);

maxx = max( CrossSection );

StartPoint = find (CrossSec == maxx );

StartPoint = StartPoint(1,1);

[ImgCen,radiii] = imfindcircles(img(:,:,StartPoint),[7 9],'Sensitivity',0.97);

ImgCen1 = ImgCen;


for dim=StartPoint:r
    
    
    %              figure, imshow( imgT2(:,:,dim),[]);
    
    centers = [];
    radii = [];
    
    
    dis1= [];
    [centers,radii,metric] = imfindcircles(img(:,:,dim),[6 12],'Sensitivity',0.99);
    num = numel(centers);
    
    if(num >0)
        
        
        for i= 1: num/2
            dis1(i,1) = norm ( ImgCen - centers(i,:) );
        end
        
        MinDis = min(dis1(:,1));
        [xMinDis yy]= find( dis1(:,1) ==MinDis);
        
        
        cen =  centers(xMinDis, :);
        radi = radii(xMinDis,:)+30;
        
        
        ImgCen = cen;
        
        %                  viscircles(cen,radi,'EdgeColor','r');
        
        mask = createCirclesMask(img(:,:,dim) ,cen,radi);
        
        cc= round(ImgCen1(1));
        
        
        for k=1 :a11
            for l=1 :b11
                if(   mask(k,l)==1 )
                    
                    mask3(k,l,dim)= imgT2(k,l,dim);
                end
            end
        end
        
        
        
    end
end

for dim=StartPoint:-1:1
    
    
    %                 figure, imshow( imgT2(:,:,dim),[]);
    
    centers = [];
    radii = [];
    dis1= [];
    [centers,radii,metric] = imfindcircles(img(:,:,dim),[6 12],'Sensitivity',0.99);
    num = numel(centers);
    
    
    if(num >0)
        
        for i= 1: num/2
            dis1(i,1) = norm ( ImgCen1 - centers(i,:) );
        end
        
        MinDis = min(dis1(:,1));
        [xMinDis yy]= find( dis1(:,1) == MinDis);
        
        cen =  centers(xMinDis, :);
        radi = radii(xMinDis,:)+30;
        
        
        ImgCen1 = cen;
        
        %                    viscircles(cen,radi,'EdgeColor','r');
        
        mask = createCirclesMask(img(:,:,dim) ,cen,radi);
        cc= round(cen(1));
        
        
        
        
        for k=1 :a11
            for l=1 :b11
                if(  mask(k,l)==1 )
                    
                    mask3(k,l,dim)= imgT2(k,l,dim);
                end
            end
        end
        
        
        
        
        
    end
    
    
end



niftiwrite(mask3,'Circle.nii');






%% Anisotropic Diffusion Filter  AND  k-means Algorithm


img = niftiread('Circle.nii');
[a b r] = size(img);
SpinalCord= zeros(size(img));

num_iter = 4;
delta_t = 1/44;
kappa = 20;
option = 1;
voxel_spacing = ones(3,1);
img = anisodiff3D(img, num_iter, delta_t, kappa, option, voxel_spacing);


for dim =1 : r
    
    int =( img(:,:,dim));
    disk = int(:);
    disk(disk==0)= NaN;
    X = disk;
    res = ~any(~isnan(X(:)));
    
    if ( res ~=1 )
        
        [idx,C]  = kmeans(X,3);
        idxx = reshape (idx , [a b]);
        minn  = min(C(:,1));
        maxx  = max(C(:,1));
        
        [Xcsf yy]= find( C(:,1)== maxx  );
        
        for i=1 :a
            for j=1 :b
                if (idxx(i,j) == Xcsf  )
                    
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


img = niftiread('K-SpinalCord.nii');

[a11 b11 r] = size(img);
Xcenter = round(a11/2);
Ycenter = round(b11/2);
ImgCen = [Ycenter Xcenter];

mask3=zeros(size(img));


for dim=1:r
    
    
    %              figure, imshow( imgT2(:,:,dim),[]);
    
    centers = [];
    radii = [];
    dis1= [];
    
    
    
    
    [centers,radii] = imfindcircles(img(:,:,dim),[7 9],'Sensitivity',0.99);
    num = numel(centers);
    
    if(num >0)
        
        for i= 1: num/2
            dis1(i,1) = norm ( ImgCen - centers(i,:) );
        end
        
        MinDis = min(dis1(:,1));
        [xMinDis yy]= find( dis1(:,1) ==MinDis);
        
        cen =  centers(xMinDis, :);
        radi = radii(xMinDis,:)+5;
        
        ImgCen = cen;
        
        %                  viscircles(cen,radi,'EdgeColor','r');
        
        mask = createCirclesMask(img(:,:,dim) ,cen,radi);
        
        cc= round(cen(1));
        
        cenn(dim,1) =cc;
        for k=1 :a11
            for l=1 :b11
                if(  cc>l &&  mask(k,l)==1 )
                    
                    mask3(k,l,dim)= img(k,l,dim);
                    
                    
                end
            end
        end
        
    end
end

%
% for dim= 1:r
%
%     mask33(:,:,dim) = fibermetric(mask3(:,:,dim) ,5,'ObjectPolarity', 'dark', 'StructureSensitivity',100);
%
% end


niftiwrite(mask3,'Circle2.nii');



%%


xxx = abs(diff(cenn));

xxx2 = abs(diff(xxx));

[ bb, rr] = size(xxx2);


bbb = round(bb/2);
k=1;

for i=bbb : bb
    
    if( xxx2(i,1)>4  )
        
        xxx3(k,1) = i;
        k=k+1;
        
    end
    
end

minDis = ( xxx3(1,1));


mask4 = zeros (size(mask3 ));

for i =1 :r
    
    if(i<minDis)
        
        
        mask4(:,:,i)= mask3(:,:,i);
        
    end
    
end


niftiwrite(mask4,'Circle3.nii');

%%



mid = step();


niftiwrite(mid,'Circle4.nii');


[aa , bb, r] = size(mid);


for i=1 : r
    
    
    im = (mid(:,:,i));
    
    total = bwarea(im);
    CrossSecDisk(i,1)= total;
    
    
end


for i=1 : r
    
    
    if ( i> minDis )
        
        CrossSecDisk(i,1) = 0;
    end
    
end




CrossSecDisk(CrossSecDisk==0)= NaN;

TF= islocalmin(CrossSecDisk,'MinProminence',15);

DiskVer = find(TF==1);

imgT2 = niftiread('ROI.nii');
imgT2 = permute( imgT2, [ 3 1 2 ] );



[aaa , rr] = size(DiskVer);

for i=1 : aaa
    
    
    SlicNum = DiskVer(i,1);    
    
    PxlRed =  imgT2 (:,:,SlicNum) >0 ;
     
          
    imgT2 (:,:,SlicNum)= PxlRed; 
    
       
    
end

niftiwrite(imgT2,'Circle5.nii');


%%


%         %Convert to original size of image and apply Medial filter
%
imgM=niftiread(strcat(PathName,MainImg));
info = niftiinfo(strcat(PathName,MainImg));
[aa , bb, rr] = size(imgM);
ImgSeg1 = zeros (size(imgM ));

LT = simetry -len;
RT = simetry +len;

if (aa>rr)
    
    imgT2 = permute( imgT2, [ 2 3 1] );
    [a b r] = size(imgT2);
    ImgSeg1 (1:a , 1:b , LT:RT) = imgT2(1:a , 1:b , 1:r);
    
else
    
    
    %     img = permute( imgM, [3 2 1 ] );
    
    imgT2 = imrotate(imgT2 ,90 );
    imgT2 = permute( imgT2, [2 1 3 ] );
    ImgSeg1 (LT:RT , 1:b , 1:r) = imgT2(1:a , 1:b , 1:r);
    
    
end

%         ImgSeg = medfilt3(ImgSeg);
ImgSeg = int16(ImgSeg1);
niftiwrite(ImgSeg,'Canal.nii',info);
% ImgSeg1 = logical(ImgSeg);
% ImgSeg = int16(ImgSeg1);
% niftiwrite(ImgSeg,'SpinalCordLogical.nii',info);





toc



%%




