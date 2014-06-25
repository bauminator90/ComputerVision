function [NewPoints,G,W]=NonRigidPointSet(X,Y,D,M,N,sig,w,W,G,alpha)
    NewPoints=ones(N,D);
    
    while norm(NewPoints-X)>0.2
       
        %E-Step:
        P=zeros(M,N);
        for m=1:M
            for n=1:N
                den=0;
                for k=1:M
                   den=den+exp(-1/(2*sig)*(norm(X(n,:)-Y(k,:)+G(k,:)*W)^2));
                   den=den+(2*pi*sig)^(D/2)*w/(1-w)*M/N;
                end
                
                                
                P(m,n)= exp(-1/(2*sig)*(norm(X(n,:)-Y(m,:)+G(m,:)*W)^2))/ den;
            end
        end 
     
        
        
        %M-Step:
        
        %Lösen der Gleichung. Wie geht das am besten??
        
        N_P=ones(M,1)'*P*ones(N,1);
        T=Y+G*W;
        diagonal=diag(P*ones(N,1))
        diagonal1=diag(P'*ones(M,1));
        sig=1/(N_P*D)*(trace(X'*diagonal1*X)-2*trace((P*X)'*T)+trace(T'*diagonal*T);
        
    end
end