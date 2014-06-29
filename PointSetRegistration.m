function [NewPoints,par1,par2,par3]=PointSetRegistration(X,Y,opt)
% Input:
% X - NxD is the first point set of N points with the dimensionality D.
%     This is the dataset.
% Y - MxD is the second point set of M points with the dimensionality D.
%     These are the GMM centroids
% opt - scalar to choose the option of the method.
%     1 is for the rigid transformation, 2 for the non-rigid.


[X, scale_x, t_x] = normalise(X);
[Y, scale_y, t_y] = normalise(Y);

%Initialization

w=0.01;
D=size(X,2);
N=size(X,1);
M=size(Y,1);

A = sum(X .* X, 2);
B = -2*X*Y';
C = sum(Y .* Y, 2);
sig2 = bsxfun(@plus, A, B);
sig2 = bsxfun(@plus, C', sig2);
sig2 = sum(sig2(:)) / (D*N*M);

R=eye(D);
t=zeros(D,1);
s=1;



switch opt
    case 1
        [NewPoints,par1,par2,par3]=RigidPointSet(X,Y,D,M,N,R,s,t,sig2,w);
        
        %Update the Rotation and Translation:
        R=par1;
        s=par2;
        t=par3;
              
        R=bsxfun(@rdivide, s*R*scale_y', scale_x);
        %t=bsxfun(@rdivide, s*R*t_y', scale_x)+bsxfun(@rdivide, t, scale_x)+bsxfun(@rdivide, t_x', scale_x);
        t=bsxfun(@rdivide, t, scale_x);
        par1=R;
        par3=t;
        
    case 2
        [NewPoints,par1,par2,par3]=AffinePointSet(X,Y,D,M,N,R,t,sig2,w);
    case 3
        %Initialization only used in this case:
        W=zeros(M,D);
        beta=0.5;
        alpha=0.5;
        G=zeros(M,M);
        for i=1:M
            for j=1:M
                G(i,j)=exp(-1/(2*beta^2)*norm(Y(i,:)-Y(j,:))^2);
            end
        end
                
        [NewPoints,par1,par2]=NonRigidPointSet(X,Y,D,M,N,sig2,w,W,G,alpha);
end



NewPoints = bsxfun(@plus, bsxfun(@rdivide, NewPoints, scale_x), t_x);
