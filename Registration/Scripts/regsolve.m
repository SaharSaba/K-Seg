


function [X,f]=regsolve(X,T,main,optim, mode)


if ~exist('mode','var'), mode=0; end;

X=X-main.Xgrid;

if mode,
    dux=dctn(X(:,:,:,1));
    duy=dctn(X(:,:,:,2));
    duz=dctn(X(:,:,:,3));
    
    f=0.5*main.alfa*sum(main.K(:).*(dux(:).^2+duy(:).^2+duz(:).^2));
else
   
    X=X-T*optim.tau;

    diax=main.alfa*optim.tau*main.K+1;
    dux=dctn(X(:,:,:,1))./diax;
    duy=dctn(X(:,:,:,2))./diax;
    duz=dctn(X(:,:,:,3))./diax;
    f=0.5*main.alfa*sum(main.K(:).*(dux(:).^2+duy(:).^2+duz(:).^2)); % compute the value of the regularization term
    dux=idctn(dux);
    duy=idctn(duy);
    duz=idctn(duz);
   
    X=cat(4,dux,duy,duz)+main.Xgrid;
    
    
    
end 
