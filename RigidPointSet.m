function [NewPoints,R,s,t]=RigidPointSet(X,Y,D,M,N,R,s,t,sig,w)
    NewPoints=ones(N,D);
    
    %ist diese Bedingung gut?
    while norm(NewPoints-X)>0.2
        norm(NewPoints-X)
        
        %E-Step:
        P=zeros(M,N);
        for m=1:M
            for n=1:N
                den=0;
                for k=1:M
                   den=den+exp(-1/(2*sig)*(norm(X(n,:)-(s.*((R*Y(k,:)')')+t))^2));
                   den=den+(2*pi*sig)^(D/2)*w/(1-w)*M/N;
                end
                
                                
                P(m,n)= exp(-1/(2*sig)*(norm(X(n,:)-(s.*((R*Y(m,:)')')+t))^2)) / den;
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
        A=X1'*P'*Y1;
        [U,S,V]=svd(A);
        C=eye(D);
        C(D,D)=det(U*V');
        R=U*C*V';
        
        %compute s,t and sig
        diagonal=diag(P*ones(N,1));
        s=trace(A'*R)/trace(Y1'*diagonal*Y1)';
        t=(m_x-s.*R*m_y)';
        sig=1/(N_P*D)*(trace(X1'*diagonal*X1-s*trace(A'*R)));
        
        NewPoints=s.*Y*R'+ones(M,1)*t;
    end
end