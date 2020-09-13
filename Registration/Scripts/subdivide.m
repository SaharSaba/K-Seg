
function Y=subdivide(X, M)

[mg,ng, kg, tmp, tmp]=size(X);

x=squeeze(X(:,:,:,1,:));
y=squeeze(X(:,:,:,2,:));
z=squeeze(X(:,:,:,3,:));

%% append one row and one column to X to get 2*size 

col= x(:,end,:,:)+1 ;
x=horzcat(x,col);
row= x(end,:,:,:) ;
x=vertcat(x,row);
x1 = x(:,:,1,:);
x= cat (3,x,x1);


%% append one row and one column to Y to get 2*size 

col= y(:,end,:,:) ;
y=horzcat(y,col);
row= y(end,:,:,:)+1 ;
y=vertcat(y,row);
y1 = y(:,:,1,:);
y= cat (3,y,y1);

%% z

col= z(:,end,:,:) ;
z=horzcat(z,col);
row= z(end,:,:,:) ;
z=vertcat(z,row);
h1= kg+1;
for i=1:mg+1
    for j=1:ng+1           
         z(i,j,h1,:) =h1;
       
    end 
end 


%%
[mg,ng, kg, tmp, tmp]=size(x);

M=1;
xnew=zeros(mg, 2*ng-2, kg, M);
ynew=zeros(mg, 2*ng-2, kg, M);
znew=zeros(mg, 2*ng-2, kg, M);

xfill=(x(:,1:end-1,:,:)+x(:,2:end,:,:))/2;
yfill=(y(:,1:end-1,:,:)+y(:,2:end,:,:))/2;
zfill=(z(:,1:end-1,:,:)+z(:,2:end,:,:))/2;


for i=1:ng-1
   xnew(:,2*i-1:2*i,:,:)=cat(2,x(:,i,:,:), xfill(:,i,:,:)); 
   ynew(:,2*i-1:2*i,:,:)=cat(2,y(:,i,:,:), yfill(:,i,:,:)); 
   znew(:,2*i-1:2*i,:,:)=cat(2,z(:,i,:,:), zfill(:,i,:,:));
end

x=xnew(:,1:end,:,:);
y=ynew(:,1:end,:,:); 
z=znew(:,1:end,:,:); 
%

xnew=zeros(2*mg-2, 2*ng-2, kg, M);
ynew=zeros(2*mg-2, 2*ng-2, kg, M);
znew=zeros(2*mg-2, 2*ng-2, kg, M);


xfill=(x(1:end-1,:,:,:)+x(2:end,:,:,:))/2;
yfill=(y(1:end-1,:,:,:)+y(2:end,:,:,:))/2;
zfill=(z(1:end-1,:,:,:)+z(2:end,:,:,:))/2;

for i=1:mg-1
   ynew(2*i-1:2*i,:,:,:)=cat(1,y(i,:,:,:), yfill(i,:,:,:)); 
   xnew(2*i-1:2*i,:,:,:)=cat(1,x(i,:,:,:), xfill(i,:,:,:)); 
   znew(2*i-1:2*i,:,:,:)=cat(1,z(i,:,:,:), zfill(i,:,:,:));
end

x=xnew(1:end,:,:,:);
y=ynew(1:end,:,:,:);
z=znew(1:end,:,:,:);
%
xnew=zeros(2*mg-2, 2*ng-2, 2*kg-2, M);
ynew=zeros(2*mg-2, 2*ng-2, 2*kg-2, M);
znew=zeros(2*mg-2, 2*ng-2, 2*kg-2, M);


xfill=(x(:,:,1:end-1,:)+x(:,:,2:end,:))/2;
yfill=(y(:,:,1:end-1,:)+y(:,:,2:end,:))/2;
zfill=(z(:,:,1:end-1,:)+z(:,:,2:end,:))/2;

for i=1:kg-1
   ynew(:,:,2*i-1:2*i,:)=cat(3,y(:,:,i,:), yfill(:,:,i,:)); 
   xnew(:,:,2*i-1:2*i,:)=cat(3,x(:,:,i,:), xfill(:,:,i,:)); 
   znew(:,:,2*i-1:2*i,:)=cat(3,z(:,:,i,:), zfill(:,:,i,:));
end

x=xnew(:,:,1:end,:);
y=ynew(:,:,1:end,:);
z=znew(:,:,1:end,:);

Y=cat(5, x, y, z);
Y=permute(Y,[1 2 3 5 4]);
Y=2*Y-1;