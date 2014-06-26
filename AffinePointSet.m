function [NewPoints,B,s,t]=AffinePointSet(X,Y,D,M,N,B,t,sig,w)
    NewPoints=ones(N,D);
    s=1;
    
    iter=0;
       
    while (iter<1000) & (abs(sig)>10*eps)
        iter=iter+1;
        %E-Step:

        % Fill P with all cominations of (x_n - sRy_m+t)^2.
        P = pdist2(X, bsxfun(@plus, s*Y*B', t), 'euclidean') .^2;
        % Transform every element p=exp(-1/(2*sig) * p);
        P = exp(-P/(2*sig));
        % The denominator is specific to each column. Sum over the rows.
        % Add constant term
        denom = sum(P,2) + (s*pi*sig)^(D/2)*w/(1-w)*M/N;
        assert(length(denom) == N);
        % Divide each column by the denominator for that column.
        P = bsxfun(@rdivide, P, denom);
        P=P';
        
        %M-Step:
        
        %calculate the means and subtract them.
        N_P=ones(M,1)'*P*ones(N,1);
        m_x=1/(N_P)*X'*P'*ones(M,1);
        m_y=1/(N_P)*Y'*P*ones(N,1);
        X1=X-ones(N,1)*m_x';
        Y1=Y-ones(M,1)*m_y';
        
        %compute the new rotation matrix.
        diagonal=diag(P*ones(N,1));
        
        C=(Y1'*diagonal*Y1);
        B=(X1'*P'*Y1)/C;
                
        %compute t and sig
        t=(m_x-B*m_y)';
        diagonal1=diag(P'*ones(M,1));
        sig=1/(N_P*D)*(trace(X1'*diagonal1*X1)-trace(X1'*P'*Y1*B'));
                        
        NewPoints=Y*B'+ones(M,1)*t;
    end
end