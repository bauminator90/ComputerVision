%Project: Point Set Registration

% % %Exercise 1: Standford Bunny
% % clear all;
% % addpath(['.' filesep 'data']);
% % addpath(['.' filesep 'data' filesep 'bunny' filesep 'data']);
% % 
% % %Read-in Bunnypoints
% % coord1=readplyfile('bun000.ply');
% % coord2=readplyfile('bun045.ply');
% % 
% % coord1=coord1(1:100:end);
% % coord2=coord2(1:100:end);
% % coord1=coord1(1:1203);
% % 
% % %Rigid:
% % figure;
% % hold on;
% % scatter(coord1(:,1),coord1(:,2),coord1(:,3),'.r'); title('Bunnydata (r), Other Bunnypoints(g)');
% % scatter(coord2(:,1),coord2(:,2),coord2(:,3),'.g');
% % 
% % [coord3, R, s, t] = PointSetRegistration(coord1, coord2, 1);
% % figure;
% % hold on;
% % scatter(coord1(:,1),coord1(:,2),coord1(:,3),'.r'); title('Bunnydata (r),  Aligned Bunnypoints(g)');
% % scatter(coord3(:,1),coord3(:,2),coord3(:,3),'.g');
% % 



% % %Exercise 3: Given Ridges-Matrices
% % 
% % %%RIGID:
% % 
%First image:
clear all;
addpath(['.' filesep 'data']);
load('ridges_1.mat');

figure;
hold on;
plot(xground(:,1),xground(:,2),'.r'); title('X (r), Y(g)');
plot(yground(:,1),yground(:,2),'.g');
 
figure;
[ynew, R, s,t]=PointSetRegistration(xground,yground,1);
figure;
hold on;
plot(xground(:,1),xground(:,2),'.r'); title('X (r), Ynew(g)');
plot(ynew(:,1),ynew(:,2),'.g');
% % 
% % 
% % %Second image:
% % load('ridges_2.mat');
% % 
% % figure;
% % hold on;
% % plot(xground(:,1),xground(:,2),'.r'); title('X (r), Y(g)');
% % plot(yground(:,1),yground(:,2),'.g');
% %  
% % [ynew, R, s]=PointSetRegistration(xground(1:10:end),yground(1:10:end),3);
% % figure;
% % hold on;
% % plot(xground(1:10:end,1),xground(1:10:end,2),'.r'); title('X (r), Ynew(g)');
% % plot(ynew(:,1),ynew(:,2),'.g');



