%Project: Point Set Registration
clear all;
addpath(['.' filesep 'data']);
addpath(['.' filesep 'data' filesep 'bunny' filesep 'data']);
load('ridges_1.mat');
load('ridges_2.mat');

% % AFFINE TEST
X = [100, 0 ; 0, 100; 50, 50];
Y = [200, 100; 150, 150; 100, 200];
[ynew, R, s, t] = PointSetRegistration(X, Y, 1)
% 
% figure;
% hold on;
% plot(X(:,1),X(:,2),'.r');
% hold on;
% plot(Y(:,1),Y(:,2),'.g');
% 
% figure;
% hold on;
% plot(X(:,1),X(:,2),'.r');
% hold on;
% plot(ynew(:,1),ynew(:,2),'.g');



% AFFINE TEST Image
xground=xground(3:10:2996,:);
yground=yground(3:10:2996,:);

yground = bsxfun(@plus, xground, [10, 30]);

figure;
hold on;
plot(xground(:,1),xground(:,2),'.r'); title('X (r), Y(g)');
plot(yground(:,1),yground(:,2),'.g');


[ynew, R, s, t]=PointSetRegistration(xground,yground,1);
figure;
hold on;
plot(xground(:,1),xground(:,2),'.r'); title('X (r), Ynew(g)');
plot(ynew(:,1),ynew(:,2),'.g');









%coord=readplyfile('bun000.ply');