
function im=imresize3D(im, newsiz)


siz=size(im);    % current size

if sum(siz==newsiz)~=3 % check if the new size is not the same

coef=isnan(im);  % if we've used masking through nans, save them
im(coef)=0;      % set nans to zeros

M=siz./newsiz;   % the fraction between current and old sizes 
 

% if resizing to smaller, remove the higher frequency content (low pass filtering through fft)
if sum(M<1)==0,    
    fim=fftshift(fftn(im));
    fim2=zeros(siz);
    range1=round(0.5*(1-1./M).*(siz-1)+1);
    range2=round(0.5*(1+1./M).*(siz-1)+1);
    fim2(range1(1):range2(1),range1(2):range2(2),range1(3):range2(3))=fim(range1(1):range2(1),range1(2):range2(2),range1(3):range2(3));
    im=abs(ifftn(fim2));

    % put nans back
    im(coef)=nan;
end

% prepair the coordinates of the smaller image
x=linspace(1,siz(2),newsiz(2));
y=linspace(1,siz(1),newsiz(1));
z=linspace(1,siz(3),newsiz(3));
[x,y,z]=meshgrid(x,y,z);

% linear interpolate

im= interp3(im,x,y,z,'linear');


end