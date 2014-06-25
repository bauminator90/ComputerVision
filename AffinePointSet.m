function [NewPoints,B,s,t]=AffinePointSet(X,Y,D,M,N,B,t,sig,w)
    NewPoints=ones(N,D);
    s=1;
    
    %ist das gut?
    while norm(NewPoints-X)>0.2
        
        %E-Step:
        P=zeros(M,N);
        for m=1:M
            for n=1:N
                den=0;
                for k=1:M
                   den=den+exp(-1/(2*sig)*(norm(X(n,:)-(((B*Y(k,:)')')+t))^2));
                   den=den+(2*pi*sig)^(D/2)*w/(1-w)*M/N;
                end
                
                                
                P(m,n)= exp(-1/(2*sig)*(norm(X(n,:)-(((B*Y(m,:)')')+t))^2)) / den;
            end
        end
        
        %M-Step:
        
        %calculate the means and subtract them.
        N_P=ones(M,1)'*P*ones(N,1);
        m_x=1/(N_P)*X'*P'*ones(M,1);
        m_y=1/(N_P)*Y'*P*ones(N,1);
        X1=X-ones(N,1)*m_x';
        Y1=Y-ones(M,1)*m_y';
        
        %compute the new rotation matrix.
        diagonal=diag(P*ones(N,1));
        B=(X1'*P'*Y1)*(Y1'*diagonal*Y1)^-1;
        
        %compute t and sig
        t=(m_x-B*m_y)';
        diagonal1=diag(P'*ones(M,1));
        sig=1/(N_P*D)*(trace(X1'*diagonal1*X1)-trace(X1'*P'*Y1*B'));
        
        NewPoints=Y*B'+ones(M,1)*t;
    end
end