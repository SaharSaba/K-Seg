

clear all; close all; clc;
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

% Main settings
main.similarity='CR';  % similarity measure,
main.subdivide=3;      % use 3 hierarchical levels
main.alfa =20;       % transformation regularization weight,
main.single=1;



% Optimization settings
optim.maxsteps =50 ;   % maximum number of iterations at each hierarchical level
optim.stopping_criterion =1e-5;  % tolerance (stopping criterion)
optim.tau =1;        % initial optimization step size
optim.anneal= 0.98;   % annealing rate on the optimization step


load Float.mat
im = squeeze(D) ;


load Refrence.mat
ref = squeeze(D) ;


eta_bef= check_similarity(ref , im)


[newim ,res]=register(ref,im, main, optim);


%     x= res(:,:,:,1);
%     y = res(:,:,:,2);
% %     z = res(:,:,:,3);
%  
% [m,n,z]= size(newim);
% 
% z= z-2;
% i=1;
% for k= 7: z
%     
%     imm (:,:,i,:) = newim(:,:,k,:);
%     i=i+1;
% 
% end 
% 
% o=1;
% for k=7: z
%     
%     ref1 (:,:,o,:) = ref(:,:,k,:);
%     o=o+1;
%     
% end 

eta_bef1= check_similarity(ref, newim)


%% imfuse 2 registerd volume


[x1,y1,z1] = size(ref);
[x2,y2,z2] = size(newim);

% h1 = figure ;
for k = 1:z1
    
    ref1 =ref(:,:,k);
    im1  = im(:,:,k);
    newim1 = newim(:,:,k);
    subplot(2,2,1); imshow(ref1,[]);
    subplot(2,2,2); imshow(im1 ,[]);
    subplot(2,2,3); imshowpair(ref1 ,im1);
    subplot(2,2,4); imshowpair(ref1 ,newim1);
    C = imfuse(ref1,newim1,'blend');
    M2(:,:,k)=C;
    pause(0.3);
             
end


% close(h1);
% 
% 
% %% saving registerd  image to a .mat  file and fusion 
% 
% D=zeros(x1,y1,1,z1);
%     
%     for n=1:z1
%         D(:,:,:,n)=M2(:,:,n);
%     end
%     
% save fusion D
% load fusion.mat
% registerd = squeeze(D);
% 
% registerd  = imresize3D(registerd , [384 384 136]);
% 
% 
%  figure;
% 
% 
% [x,y,z1] = size(registerd);
% filename = 'registerd';
% 
% 
% for k = 5:z1-2
% 
%     imshow(registerd(:,:,k,1),[])
%     drawnow
% 
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1;
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
% 
%     pause(.3)
% end
% 
% 
