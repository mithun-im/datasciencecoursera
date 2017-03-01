function [lb, ub] = function123(N,M,Qi,Qo,P,count,Pi,Po)
    c = randn(N,1);
    
    cvx_begin sdp 
        variable x(1,N) complex
        maximize real(x*c)
        subject to
            norm(x*Qi') <= sqrt(Pi);
            norm(x*Qo') <= sqrt(Po);
    cvx_end
    y = x;

    for k = 1:count
        cvx_begin sdp 
            variable x(1,N) complex
            maximize log_det(0.5*(kron(eye(M),x)*P*(kron(eye(M),y))' + (kron(eye(M),x)*P*(kron(eye(M),y))')'))
            subject to
                norm(x*Qi') <= sqrt(Pi);
                norm(x*Qo') <= sqrt(Po);
        cvx_end
        y = x;
    end
    
    cvx_begin
        variables X(N,N) Z(M,M)
        maximize log_det(Z)
        subject to 
            X == hermitian_semidefinite(N)
            Z == hermitian_semidefinite(M)
            for i = 1:M
                for j = i:M
                    Z(i,j) == trace(P((i-1)*N+1:i*N,(j-1)*N+1:j*N)*X)
                end
            end
            real(trace(Qi'*Qi*X)) <= Pi
            real(trace(Qo'*Qo*X)) <= Po
    cvx_end        
end