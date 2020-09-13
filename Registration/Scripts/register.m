

function [result , res]=register(refim, im, main, optim)

tic

dimen=size(refim);
M=ceil(dimen/2^(main.subdivide-1));

[x,y,z]=meshgrid(1:M(1),1:M(2),1:M(3));

[main.mg, main.ng , main.kg]=size(x);

main.siz=[main.mg main.ng main.kg];

main.X=cat(4,x,y,z); 
main.Xgrid=main.X; 

for level=1:main.subdivide
    
    if level==2
        main.mg=2*main.mg; 
        main.ng=2*main.ng;
        main.kg=2*main.kg;
        main.alfa =300;       

        
    elseif   level==3 
        main.mg=2*main.mg; 
        main.ng=2*main.ng;
        main.kg=2*main.kg;
        main.alfa =200;       

        
    elseif   level==4
        main.mg=2*main.mg; 
        main.ng=2*main.ng;
        main.kg=2*main.kg; 
        main.alfa =100;       
        
    end
        
    main.level=level;
    main.siz=[main.mg main.ng main.kg];
    main.K=initK(main.siz); % Init Laplacian eigenvalues (used for regularization)
       
    main.refimsmall=imresize3D(refim,main.siz); % resize images
%     main.refimsmall(main.refimsmall<0)=0;  
%     main.refimsmall(main.refimsmall>1)=1;
     
    imsmall=imresize3D(im,main.siz); 
%     imsmall(imsmall<0)=0;
%     imsmall(imsmall>1)=1;
    
    [gradx, grady ,gradz]=gradient(imsmall);
    main.imsmall=cat(4,imsmall, gradx,grady,gradz); 
    
    
    [main.X, result]=registration(main.X, main, optim);

    if level<main.subdivide,
        main.X =subdivide(main.X, 1);
        main.Xgrid=subdivide(main.Xgrid, 1);
    end;
%     
end
res=main.X;

disptime(toc);
