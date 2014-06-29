function [NewPoints,G,W]=NonRigidPointSet(X,Y,D,M,N,sig2,w,W,G,alpha)
    NewPoints=ones(N,D);
    
    iter=0;
       
    while (iter<1000) & (abs(sig2)>10*eps)
        iter=iter+1;
              
        %E-Step:

        % Fill P with all cominations of (x_n - sRy_m+t)^2.
        T = bsxfun(@plus, Y,G*W);
        A = sum(X .* X, 2);
        B_m = -2*X*T';
        C = sum(T .* T, 2);
        P = bsxfun(@plus, A, B_m);
        P = bsxfun(@plus, C', P);
        
     
        % Transform every element p=exp(-1/(2*sig2) * p);
        P = exp(-P/(2*sig2));
        % The denominator is specific to each column. Sum over the rows.
        % Add constant term
        denom = sum(P,2) + (2*pi*sig2)^(D/2)*w/(1-w)*M/N;
        assert(length(denom) == N);
        % Divide each column by the denominator for that column.
        P = bsxfun(@rdivide, P, denom);
        P=P';
        
        %M-Step:
        diagonal=diag(P*ones(N,1));
        diagonal1=diag(P'*ones(M,1));
        
        A=(diagonal*G+alpha*sig2);%*((diagonal)^-1));
        B=P*X-diagonal*Y;
        %B=((diagonal)^-1)*P*X-Y;
        W=A\B;
         
        
        N_P=ones(M,1)'*P*ones(N,1);
        T=Y+G*W;

        sig2=1/(N_P*D)*(trace(X'*diagonal1*X)-2*trace((P*X)'*T)+trace(T'*diagonal*T));

        NewPoints=T;
        
    end
end