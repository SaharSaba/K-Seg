function KC_seg(imgM, info)

tic
[a , b, r] = size(imgM);
img = imgM;

if (a<b)
    img = permute( img, [ 2 3 1] );
end

MatSim= zeros(r,r);
[a , b, r] = size(img);
k=1;
len=12;

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

BW1 = edge3 (img,'approxcanny',0.1 );

for dim= 1:r
    
    
    re = (img(:,:,dim));
    BW= BW1(:,:,dim);
    [H,T,R] = hough(BW);
    P  = houghpeaks(H,50);
    lines = houghlines(BW,T,R,P );
    [x ,y ] = size(lines);
    L= bwlabel(BW);
    
    [m ,n] = size(L);
    
    h=1;
    for i = 1 :y
        
        
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
    
    for i = 1 :x
        for j = 1 :y
            for k = 1 :T
                
                if (L(i,j)== ind (1,k))
                    
                    L2(i,j)=1;
                    
                end
            end
        end
    end
    
    
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
img = permute( img, [ 3 1 2 ] );
imgT2 = permute( imgT2, [ 3 1 2 ] );
[a11 b11 r] = size(img);
Xcenter = round(a11/2);
Ycenter = round(b11/2);
% ImgCenMain = [Ycenter Xcenter];
mask2=zeros(size(img));


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
    
    if (CrossSec(i,1)>100 )
        
        CrossSection(i,1)= CrossSec(i,1);
    else
        
        CrossSection(i,1)=0;
    end
    
end

CrossSection(CrossSection==0)= NaN ;

AvgCross = nanmean (CrossSection);

minn = min( CrossSection );

StartPoint = find (CrossSec == minn );

[ImgCen,radiii] = imfindcircles(img(:,:,StartPoint),[7 9],'Sensitivity',0.97);

ImgCen1 = ImgCen;

for dim=StartPoint:r
    
    
    %  figure, imshow( imgT2(:,:,dim),[]);
    
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
        radi = radii(xMinDis,:);
        
        if (AvgCross >200 )
            radi = radi+5;
            
        end
        
        
        ImgCen = cen;
        
        %  viscircles(cen,radi,'EdgeColor','r');
        
        mask = createCirclesMask(img(:,:,dim) ,cen,radi);
        
        for k=1 :a11
            for l=1 :b11
                if(mask(k,l)==1)
                    
                    mask2(k,l,dim) = imgT2(k,l,dim);
                    
                end
            end
        end
    end
end

for dim=StartPoint:-1:1
    
    
    %     figure, imshow( imgT2(:,:,dim),[]);
    
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
        radi = radii(xMinDis,:);
        
        
        if (AvgCross >200 )
            radi = radi+5;
            
        end
        
        
        
        ImgCen1 = cen;
        
        %    viscircles(cen,radi,'EdgeColor','r');
        
        mask = createCirclesMask(img(:,:,dim) ,cen,radi);
        for k=1 :a11
            for l=1 :b11
                if(mask(k,l)==1)
                    
                    mask2(k,l,dim) = imgT2(k,l,dim);
                end
            end
        end
    end
end


niftiwrite(mask2,'Circle.nii');


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
        
        [idx,C]  = kmeans(X,3,'Replicates',10);
        idxx = reshape (idx , [a b]);
        minn  = min(C(:,1));
        maxx  = max(C(:,1));
        
        [Xcsf yy]= find( C(:,1)~= minn &  C(:,1)~= maxx  );
        
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


%% Remove extra points AND cleaning in sagittal view


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
    SCax(:,:,dim)= BW;
    
end


for i=1 : r
    
    im = (SCax(:,:,i));
    total = bwarea(im);
    Cross (i,1)= total;
    
    
end


niftiwrite(SCax,'SCAxial.nii');


%%
%Convert to original size of image and apply Medial filter

% imgM=niftiread(strcat(PathName,MainImg));
% info = niftiinfo(strcat(PathName,MainImg));

[aa , bb, rr] = size(imgM);
ImgSeg = zeros (size(imgM ));

LT = simetry -len;
RT = simetry +len;

if (aa>rr)
    SCax = permute( SCax, [ 2 3 1] );
    [a b r] = size(SCax);
    ImgSeg (1:a , 1:b , LT:RT) = SCax(1:a , 1:b , 1:r);
    
else
    ImgSeg (LT:RT , 1:b , 1:r) = SCax(1:a , 1:b , 1:r);
end

ImgSeg = medfilt3(ImgSeg);
ImgSeg = int16(ImgSeg);
niftiwrite(ImgSeg,'SpinalCord.nii',info);
ImgSeg1 = logical(ImgSeg);
out = int16(ImgSeg1);
niftiwrite(out,'SpinalCordLogical.nii',info);


toc


end


