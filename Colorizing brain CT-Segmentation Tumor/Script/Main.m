clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;
fontSize = 20;

grayImage = dicomread('2.DCM');

imshow(grayImage);


imgHdr = dicominfo('2.dcm');

[m, n, x] = size(grayImage);

set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'name' ,'Demo by ImageAnalyst','numbertitle','off') 


%///////////////////////////////////////////////////////


%  »œÌ· ⁄‰«’— HU »Â ÅÌò”· 


huMat = int16(zeros (512));

for i=1 :m
    for j=1 : n 
        
        pixle_value = int16(grayImage(i,j));
        
        huPixel =  pixle_value -  1024 ; %  hu = pixel_value * slope + intercept
        
        huMat(i,j) =huPixel;
        
    end
    
end 



subplot(1, 2, 1);

%imshow(grayImage, []);
imshow(huMat, []);
title('original gray scale image', 'FontSize', fontSize);
uiwait(msgbox('Draw the ROI'));


roi = imfreehand();  
binaryImage= roi.createMask();  


roival = huMat(binaryImage);

maxroi =  max (roival);
minroi =  min (roival);

d= maxroi - minroi ;

tumorCell = zeros(size(huMat));


for  i=1 :m
 
    for j= 1 :n
        
        pixel=huMat(i,j);
       
        
         if ((pixel<=maxroi) & ( pixel>=minroi)) 

              new_pixel=pixel;
        
         else
               new_pixel = 0;
         
        end
             
         tumorCell(i,j)=new_pixel;

      end
end



binaryImage1 = binaryImage | tumorCell  ;

subplot(1, 2, 2);
imshow(tumorCell, []);
title('GrayImage TumorCell ', 'FontSize', fontSize);


%subplot(2, 2, 3);
%imshow(binaryImage, []);
%title('Your ROI', 'FontSize', fontSize);



%//////////////////////////////////

% «ÌÃ«œ  ’ÊÌ— —‰êÌ 

IMG = zeros(m,n,3);

for  ii=1 :m
 
    for jj= 1 :n

         
       if tumorCell(ii,jj)~=0 
        
           a= tumorCell(ii,jj);
           a= double(a);
           d= double(d);
                     
           v = a./d ;
                    
         if( v < 0.5)
                     
              
            p1=2.*v;
            p2=1-p1;
            p3=255.*p2;% red
            p4=255.*p1;% green
            
            if (p3>p4)              

                
                IMG(ii,jj,1) = p3;
            
                IMG(ii,jj,2) = 0;

                IMG(ii,jj,3) = 0;
            
            elseif (p3<p4)
                
                
                
                IMG(ii,jj,1) = 0;
            
                IMG(ii,jj,2) = p4;

                IMG(ii,jj,3) = 0;
            end
            
         elseif(v>=0.5)
             
             
            p1=v-0.5;
            p2= p1.*2;
            p3= 1- p2;
            p4 = 255.*p3;
            p5 = v-0.5;
            p6=p5.*2;
            p7=255.*p6;
            
            
            if (p4>p7)
                
                IMG(ii,jj,1)=0;
            
                IMG(ii,jj,2)=p4;
            
                IMG(ii,jj,3)=0;
                
                
            elseif(p4<p7)
                
                IMG(ii,jj,1)=0;
            
                IMG(ii,jj,2)=0;
            
                IMG(ii,jj,3)=p7;
            
            end
            
        end
        end
    end
end


rgbImage = cat(3, IMG(:,:,1),IMG(:,:,2), IMG(:,:,3));
figure,subplot(1,2,1) , imshow(rgbImage,[]);
title('RGB Image ', 'FontSize', fontSize);




%/////////////////////////////////////////////////////

%  ⁄ÌÌ‰ Ìò ‰ﬁÿÂ œ—  ’ÊÌ— —‰êÌ »Â ⁄‰Ê«‰   seed 


% Label the binary image and computer the centroid and center of mass.
%labeledImage = bwlabel(binaryImage);
%measurements = regionprops(binaryImage, grayImage, ...
  %  'area', 'Centroid', 'WeightedCentroid');
%area = measurements.Area;%  ⁄œ«œ ò· ÅÌò”· Â« œ —  roi
%centroid = measurements.Centroid;% „—ò“ Àﬁ·  roi 


%centerOfMass = measurements.WeightedCentroid %  „—ò“ œ—  roi 

%r= centerOfMass(: , 1)

%c= centerOfMass(: , 2)
 
 [c , r , p  ] = impixel(rgbImage);
 
 seedImage = rgbImage ;
 
for i=1 : m
    for j=1 :n
        
        if(  i== r && j==c )
            
            
            seedImage(i,j,: ) = 0;
                               
        end
        
    end
end
 


seedImage = cat(3, seedImage(:,:,1),seedImage(:,:,2), seedImage(:,:,3));
subplot(1,2,1) ,  imshow(seedImage,[]);
title('Seed Point Image ', 'FontSize', fontSize);

 
%////////////////////////////////////////////////////////////

 
 
%% Segmentation  part 


% „Õ«”»Â ›«’·Â seed  « Â— ÅÌò”·

D= zeros(m,n);

for i= 1 : m 
    for j=1 : n
        
        
        Px =i ;Py=j ; 
        xSeed = r; ySeed = c;
        D(i,j) = sqrt((Px-xSeed)^2 + (Py-ySeed)^2); % pythagoras theorem
               
    end
end



%//////////////////////////////////////////////

redChannel = rgbImage(:, :, 1);
greenChannel = rgbImage(:, :, 2);
blueChannel = rgbImage(:, :, 3);

% Create a new color channel images for the output.
outputImageR = tumorCell;
outputImageG = tumorCell;
outputImageB = tumorCell;


% Transfer the colored parts.
outputImageR(binaryImage) = redChannel(binaryImage);
roiRed = outputImageR(binaryImage);
sumR = sum (roiRed);
Ar = sumR ./ numel (roiRed);


outputImageG(binaryImage) = greenChannel(binaryImage);
roiGreen = outputImageG(binaryImage);
sumG = sum (roiGreen);
Ag = sumG ./ numel (roiGreen);



outputImageB(binaryImage) = blueChannel(binaryImage);
roiBlue = outputImageB(binaryImage);

sumB = sum (roiBlue);
Ab = sumB ./ numel (roiBlue);

% Convert into an RGB image
% Ãœ« ò—œ‰  ‰Â« „‰ÿﬁÂ  roi  »’Ê—  —‰êÌ 




outputRGBImage = cat(3, outputImageR, outputImageG, outputImageB);
 %subplot(1,2,2), imshow(outputRGBImage,[]);
 %title('RGB ROI ', 'FontSize', fontSize);

 
%///////////////////////////////////////////////////// 

% caculate color factor

colorMatris = zeros (m,n);

Pr = zeros(m,n);
Pg = zeros(m,n);
Pb = zeros(m,n);


for i=1 : m 
    for j = 1 : n 
        
        Pr(i,j) = rgbImage(i, j, 1);
       
        Pg(i,j) = rgbImage(i, j, 2);

        Pb(i,j) = rgbImage(i, j, 3);
        
        colorPixel= ((Pr(i,j)-Ar).^2 + (Pg(i,j)-Ag).^2 +(Pb(i,j) - Ab).^2 ).^0.5;
        
        colorMatris(i,j) = colorPixel./255;
        
       
        
    end
end


tumorPixels = zeros(m,n);
Tp =zeros(m,n);

for i=1 : m
    for j=1 : n

        reg_mean = (r+c)./2;
        Tp(i,j) = (colorMatris(i,j).* D(i,j))./reg_mean ;
        
        if ( Tp(i,j)<= 0.1)
            
            
             tumorPixel =Tp(i,j);
              
                   
         else
               tumorPixel = 0;
         
        end
         tumorPixels(i,j)=tumorPixel;

        end
        
end 
 
segmentImage = tumorCell & tumorPixels ;
%figure , subplot(1,2,1), imshow(segmentImage,[]);
%title('Segmented Tumor ', 'FontSize', fontSize);


%////////////////////////////////////////////////////

% colorize segmented roi

seg_RGB = zeros(m,n,3); 

for i=1: m 
    for j= 1: n 
        
        if (segmentImage(i,j) ~= 0)
            
            seg_RGB(i,j,1)=255;
      
            seg_RGB(i,j,2)=0;
      
            seg_RGB(i,j,3)=0;
        
        end
    end 
end  

seg_RGB = cat(3, seg_RGB(:,:,1),seg_RGB(:,:,2), seg_RGB(:,:,3));

%subplot(1,2,2) , imshow(seg_RGB,[]);
%title('segmentRGB ', 'FontSize', fontSize);

redChannel = seg_RGB(:, :, 1);
greenChannel = seg_RGB(:, :, 2);
blueChannel = seg_RGB(:, :, 3);

% Create a new color channel images for the output.

grayImage1 = im2double(grayImage);
outputImageR = grayImage1;
outputImageG = grayImage1;
outputImageB = grayImage1;

% Transfer the colored parts.
outputImageR(segmentImage) = redChannel(segmentImage);

outputImageG(segmentImage) = greenChannel(segmentImage);

outputImageB(segmentImage) = blueChannel(segmentImage);

outputRGBseg = cat(3, outputImageR, outputImageG, outputImageB);
subplot(1,2,2), imshow(outputRGBseg,[]);
title('RGBsegment ', 'FontSize', fontSize);

 
 
 