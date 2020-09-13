

clear all; close all; clc;


% load fisheriris

meas = [randn(100,3)*0.75-ones(100,3);
    randn(100,3)*0.8+ones(100,3)];
    



[cidx2,cmeans2] = kmeans(meas,3,'dist','sqeuclidean');
% [silh2,h] = silhouette(meas,cidx2,'sqeuclidean');

ptsymb = {'g*','r+','b^'};
for i = 1:3
    clust = find(cidx2==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i},'MarkerSize',8, 'LineWidth',1.5 );
    hold on
end
plot3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'ko','MarkerSize',12,'LineWidth',2);
% plot3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'ko','MarkerSize',10,'LineWidth',2);
% plot3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'ko','MarkerSize',10,'LineWidth',2);
legend({'Cluster 1','Cluster 2','Cluster 3','Centroids'},'FontSize',14);

xlabel('Dim 1','FontSize',14,'FontWeight','bold');
ylabel('Dim 3','FontSize',14,'FontWeight','bold');
zlabel('Dim 2', 'FontSize',14,'FontWeight','bold');

title('K-Means Clusters','FontSize',18);
set(gca,'FontSize',15);



view(-137,10);
grid on