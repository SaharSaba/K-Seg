function [score,  y, x ] = Generalized_hough_transform(Is,Itm) 
%Find template/shape Itm in greyscale image Is using generalize hough trasform

%Use generalized hough transform to find Template/shape binary image given in binary image Itm in grayscale image Is (greyscale image)

%Is is greyscale  picture were the template Itm should be found 
%Itm is bool edge image of the template with edges markedd ones
%Return the x,y location  cordniates  which gave the best match 
%Also return the score of each this point (number of point matching)

%The x,y are the cordinates in image Is in which the  the top left edge of image Itm (1,1) should be positioned in order to give the best match

%Is=imread('');
%Itm=imread('');
%if nargin<3 thresh=1;end;

%--------------------------create edge and system edge images------------------------------------------------------------------------------------------------------------------------

%Is=rgb2gray(Is);

Iedg=edge(Is,'canny'); % Take canny edge images of Is with automatic threshold
%}
%--------------------------------------------------------------------------------------------------------------------------------------
[y x]=find(Itm>0); % find all y,x cordinates of all points equal 1 inbinary template image Itm
nvs=size(x);% number of points in the  template image
%-------------------Define Yc and Xc ----------------------------------------------
Cy=1;%round(mean(y));% find object y center, note that any reference point will do so the origin of axis hence 1 could be used just as well
Cx=1;%round(mean(x));% find object z center, note that any reference point will do so the origin of axis hence 1 could be used just as well

%------------------------------create gradient map of Itm, distrobotion between zero to pi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GradientMap = gradient_direction( Itm );

%%%%%%%%%%%%%%%%%%%%%%%Create an R-Table of Itm gradients to  parameter space in parameter space.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------create template descriptor array------------------------------------
MaxAngelsBins=30;% devide the angel space to MaxAngelsBins uniformed space bins
MaxPointsPerangel=nvs(1);% maximal amount of points corresponding to specific angel

PointCounter=zeros(MaxAngelsBins);% counter for the amount of edge points associate with each angel gradient
Rtable=zeros(MaxAngelsBins,MaxPointsPerangel,2); % assume maximum of 100 points per angle with MaxAngelsBins angles bins between zero and pi and x,y for the vector to the center of each point
% the third adimension are vector from the point to the center of the vessel

%------------------fill the  angel bins with points in the Rtable---------------------------------------------------------
for f=1:1:nvs(1)
    bin=round((GradientMap(y(f), x(f))/pi)*(MaxAngelsBins-1))+1; % transform from continues gradient angles to discrete angle bins and 
    PointCounter(bin)=PointCounter(bin)+1;% add one to the number of points in the bin
    if (PointCounter(bin)>MaxPointsPerangel)
        disp('exceed max bin in hugh transform');
    end;
    Rtable(bin, PointCounter(bin),1)= Cy-y(f);% add the vector from the point to the object center to the bin
    Rtable(bin, PointCounter(bin),2)= Cx-x(f);% add the vector from the point to the object center to the bin
end;
%plot(pc);
%pause;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%create and populate hough space%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------use the array in previous image to identify the template in the main image Is----------------------------------------
[y x]=find(Iedg>0); % find all edg point in the in edge image Iedg of the main image Is
np=size(x);% find number of edge points Is edge image
GradientMap=gradient_direction(Is); % create gradient direction  map of the Is
Ss=size(Is); % Size of the main image Is
houghspace=zeros(size(Is));% the hough space assume to be in size of the image but it should probably be smaller
    for f=1:1:np(1)
          bin=round((GradientMap(y(f), x(f))/pi)*(MaxAngelsBins-1))+1; % transform from continues gradient angles to discrete angle bins and 

          for fb=1:1:PointCounter(bin)
              ty=Rtable(bin, fb,1)+ y(f);
              tx=Rtable(bin, fb,2)+ x(f);
               if (ty>0) && (ty<Ss(1)) && (tx>0) && (tx<Ss(2))  
                   houghspace(Rtable(bin, fb,1)+ y(f), Rtable(bin, fb,2)+ x(f))=  houghspace(Rtable(bin, fb,1)+ y(f), Rtable(bin, fb,2)+ x(f))+1; % add point in were the center of the image should be according to the pixel gradient
               end;        
          end;
    end;

%{
%====================================show The Hough Space in color==================================================================================================
imtool(houghspace);
imshow(houghspace,[]);
colormap jet
colorbar
pause
%}

%============================================Find best match in hough space=========================================================================================

%---------------------------------------------------------------------------normalized according to template size (fraction of the template points that was found)------------------------------------------------------------------------------------------------
Itr=houghspace;%./(sum(sum(Itm))); % Itr become the new score matrix
%---------------------------------------------------------------------------find  the location best score all scores which are close enough to the best score
%imtool(Itr,[]);
mx=max(max(Itr));% find the max score location
[y,x]=find(Itr==mx);% 


%[y,x]=find(Itr>=thresh*mx,  10, 'first'); % find the location first 10 best matches which their score is at least thresh percents of the maximal score and pot them in the x,y array

score=Itr(y,x); % find max score in the huogh space 

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ Is ] = gradient_direction( i3 )
%return map of the absolute direction from -pi/2 to pi/2  of gradient in every point of the gradient  i3(only half circle does not have negative directionss
%-------------------------------------------------------------------
Dy=imfilter(double(i3),[1; -1],'same');%x first derivative  sobel mask
Dx=imfilter(double(i3),[1  -1],'same');% y sobel first derivative
%Is=atan2(Dy,Dx)+pi();
Is=mod(atan2(Dy,Dx)+pi(), pi());%atan(Dy/Dx);%note that this expression can reach infinity if dx is zero mathlab aparently get over it but you can use the folowing expression instead slower but safer: 
%mod(atan2(Dy,Dx)+pi(), pi());%gradient direction map going from 0-180
%--------------------show the image-----------------------------------------------
%{
imshow(Is,[]);% the ,[]  make sure the display will be feeted to doube image
colormap jet
colorbar
pause;
%}
end

