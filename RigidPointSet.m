function [NewPoints,R,s,t]=RigidPointSet(X,Y,D,M,N,R,s,t,sig2,w)
    NewPoints=ones(N,D);
    iter=0;
       
        gif = 'animatedgif.gif';
    while (iter<1000) && (abs(sig2)>1e-5)
        iter=iter+1
              
        %E-Step:

        % Fill P with all cominations of (x_n - sRy_m+t)^2.
        T = bsxfun(@plus, s*R*Y', t)';
        A = sum(X .* X, 2);
        B = -2*X*T';
        C = sum(T .* T, 2);
        P = bsxfun(@plus, A, B);
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
        
        %calculate the means and subtract them.
        N_P=ones(M,1)'*P*ones(N,1);
        m_x=1/(N_P)*X'*P'*ones(M,1);
        m_y=1/(N_P)*Y'*P*ones(N,1);
        X1=X-ones(N,1)*m_x';
        Y1=Y-ones(M,1)*m_y';
        
        %compute the new rotation matrix.
        A=X1'*P'*Y1;
        [U,~,V]=svd(A);
        C=eye(D);
        C(D,D)=det(U*V');
        R=U*C*V';
        
        %compute s,t and sig2
        diagonal=diag(P*ones(N,1));
        diagonal1=diag(P'*ones(M,1));
        s=trace(A'*R)/trace(Y1'*diagonal*Y1);
        t=(m_x-s.*R*m_y);
        sig2=1/(N_P*D)*(trace(X1'*diagonal1*X1)-s*trace(A'*R));
        
        if (sig2 < -1e-10)
            fprintf('Oops, shouldn''t have happened. sig2: %f\n', sig2);
            assert(sig2 > 0);
        end
        
        NewPoints=s.*Y*R'+ones(M,1)*t';
%         if (mod(iter,100)==1)
%             Nrm=normalise(X);
%             plot(Nrm(:,1),Nrm(:,2),'.r');
%             plot(NewPoints(:,1),NewPoints(:,2),'.g');
%             title( sprintf('k = %d', iter) );
%     
%             Copy the new image in the gif.
%             frame = getframe(1);
%             new_image = frame2im(frame);
%             [ind_image,colormap] = rgb2ind(new_image,256);
%             if iter == 1
%                 imwrite(ind_image,colormap,gif,'gif','LoopCount',Inf,'DelayTime',1);
%             else
%                 imwrite(ind_image,colormap,gif,'gif','WriteMode','append','DelayTime',1);
%             end
%         end
    end
 
end