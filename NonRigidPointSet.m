function [NewPoints,G,W]=NonRigidPointSet(X,Y,D,M,N,sig,w,W,G,alpha)
    NewPoints=ones(N,D);
    
    iter=0;
       
    while (iter<1000) & (abs(sig)>10*eps)
        iter=iter+1;
              
        %E-Step:

        % Fill P with all cominations of (x_n - sRy_m+t)^2.
        P = pdist2(X, bsxfun(@plus, Y, G*W), 'euclidean') .^2;
        % Transform every element p=exp(-1/(2*sig) * p);
        P = exp(-P/(2*sig));
        % The denominator is specific to each column. Sum over the rows.
        % Add constant term
        denom = sum(P,2) + (2*pi*sig)^(D/2)*w/(1-w)*M/N;
        assert(length(denom) == N);
        % Divide each column by the denominator for that column.
        P = bsxfun(@rdivide, P, denom);
        P=P';
        
        %M-Step:
        diagonal=diag(P*ones(N,1));
        diagonal1=diag(P'*ones(M,1));
        
        A=(G+alpha*sig*((diagonal)^-1));
        B=((diagonal)^-1)*P*X-Y;
        W=A\B;
        
        
        
        N_P=ones(M,1)'*P*ones(N,1);
        T=Y+G*W;

        sig=1/(N_P*D)*(trace(X'*diagonal1*X)-2*trace((P*X)'*T)+trace(T'*diagonal*T));

        NewPoints=T;
        
    end
end