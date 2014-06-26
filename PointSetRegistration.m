function [NewPoints,par1,par2,par3]=PointSetRegistration(X,Y,opt)
% Input:
% X - NxD is the first point set of N points with the dimensionality D.
%     This is the dataset.
% Y - MxD is the second point set of M points with the dimensionality D.
%     These are the GMM centroids
% opt - scalar to choose the option of the method.
%     1 is for the rigid transformation, 2 for the non-rigid.


%Initialization

D=size(X,2);
N=size(X,1);
M=size(Y,1);
R=eye(D);
t=zeros(1,D);
s=1;

w=0.5;

A = sum(X .* X, 2);
B = -2*X*Y';
C = sum(Y .* Y, 2);
sig = bsxfun(@plus, A, B);
sig = bsxfun(@plus, C', sig);
sig = sum(sig(:)) / (D*N*M);

switch opt
    case 1
        [NewPoints,R,s,t]=RigidPointSet(X,Y,D,M,N,R,s,t,sig,w);
    case 2
        [NewPoints,B,s,t]=AffinePointSet(X,Y,D,M,N,R,t,sig,w);
    case 3
        %Initialization only used in this case:
        W=zeros(M,D);
        beta=0.5;
        alpha=0.5;
        G=zeros(M,M);
        for i=1:M
            for j=1:M
                G(i,j)=exp(-1/(2*beta^2)*norm(Y(i,:)-Y(j,:))^2)
            end
        end
                
        [NewPoints,G,W]=NonRigidPointSet(X,Y,D,M,N,sig,w,W,G,alpha);
end

