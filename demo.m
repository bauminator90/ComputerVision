%Project: Point Set Registration
clear all;
addpath(['.' filesep 'data']);
addpath(['.' filesep 'data' filesep 'bunny' filesep 'data']);
load('ridges_1.mat');
load('ridges_2.mat');

%rigid transformation
ynew=PointSetRegistration(xground,yground,1);

X = [1, 0; 0, 1];
Y = [1, 0; 0, 0; 0, 1];
Y = X;
ynew = PointSetRegistration(X, Y, 1);


% alpha=25;
% Rot=[1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
% x=[1,1,1;1,2,3;4,3,6];
% y=Rot*x+[5,19,1;5,19,1;5,19,1];
% 
% [ynew,B,s,t]=PointSetRegistration(x,y,3);

coord=readplyfile('bun000.ply');