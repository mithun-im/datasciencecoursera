    ## This is a markdown file
clc
clear all

N = 2;
M = 2;

for i = 1:N
    for j = 1:M
        od = randi([1,1],1,1);
        on = od - 1;
        rd = 0.25*randi([-5,-1],1,od) + 1i*randi([-1,0],1,od);
        rd = unique([rd conj(rd)]);
        rn = 0.5*randi([-5,-1],1,on);
        pn(i,j).mat = poly(rn); pd(i,j).mat = poly(rd);
        G(i,j) = (pd(i,j).mat(end)/pn(i,j).mat(end))*tf(pn(i,j).mat,pd(i,j).mat);
    end
end


for i = 1:N
    for j = 1:M
        BW(i,j) = bandwidth(G(i,j));
    end
end

f = ceil(max(BW(:)));

C = [];
Cval = [];
syms s
syms w real

for  i = 1:N
    for j = 1:M
        n = sym(strcat('cn_',num2str(i),num2str(j),'_'), [1 length(pn(i,j).mat)],'real');
        d = sym(strcat('cd_',num2str(i),num2str(j),'_'), [1 length(pd(i,j).mat)-1],'real');
        C = [C n d];
        Cval = [Cval (pd(i,j).mat(end)/pn(i,j).mat(end))*pn(i,j).mat pd(i,j).mat(2:end)];
        GS(i,j) = poly2sym(n,s)/(poly2sym([1 d],s));
    end
end

for i = 1:length(C)
    GD(i).res = subs(diff(GS,C(i)),s,1i*w);
    GD(i).res = subs(GD(i).res,C,Cval);
end

GS_sub = subs(GS,s,1i*w);
GS_sub = subs(GS_sub,C,Cval);


L = 4;
B = chebyshevT([0:L], w);
B = subs(B, w, w/(2*pi*f));
H = [];
for i = 1:L
    H = [H sym(strcat('H',num2str(i),'_'), [N,N], 'real')];
end

W = [];
for i = 1:L
    W = [W B(i)*eye(N)];
end


Z = simplify(W'*W);
temp = partfrac(Z, 'FactorMode', 'complex');
temp1 = int(temp);
Integrand = eval((subs(temp1,w,2*pi*f) - subs(temp1,w,-2*pi*f)));
InputPowerMat = [];
for l = 1:N
    InputPowerMat = blkdiag(InputPowerMat,Integrand);
end
InputPowerMat = 0.5*(InputPowerMat+InputPowerMat');
InputPowerMat(abs(InputPowerMat)<=1e-6) = 0;
eig(InputPowerMat)'

Z = simplify(W'*(GS_sub'*GS_sub)*W);
temp = partfrac(Z, 'FactorMode', 'complex');
temp1 = int(temp);
Integrand = eval((subs(temp1,w,2*pi*f) - subs(temp1,w,-2*pi*f)));
OutputPowerMat = [];
for l = 1:N
    OutputPowerMat = blkdiag(OutputPowerMat,Integrand);
end
OutputPowerMat = 0.5*(OutputPowerMat+OutputPowerMat');
OutputPowerMat(abs(OutputPowerMat)<=1e-6) = 0;
eig(OutputPowerMat)'


for i = 1:length(C)
    for j = 1:length(C)
        Z = simplify(W'*GD(j).res'*GD(i).res*W);
        temp = partfrac(Z, 'FactorMode', 'complex');
        temp1 = int(temp);
        Integrand = eval((subs(temp1,w,2*pi*f) - subs(temp1,w,-2*pi*f)));
        Msym(i,j).mat = [];
        for l = 1:N
            Msym(i,j).mat = blkdiag(Msym(i,j).mat,Integrand);
        end
        [i j]
    end            
end

M = zeros(length(C)*length(H(:)),length(C)*length(H(:)));
l = length(C);
h = length(H(:));
for i = 1:l
    for j = 1:l
        M((i-1)*h+1:i*h,(j-1)*h+1:j*h) = Msym(i,j).mat;
    end
end
M = 0.5*(M+M');


[lb, ub, x] = function123(length(H(:)),length(C),chol(InputPowerMat),chol(OutputPowerMat),P,20,100,100)
H = reshape(x,[M L*M]);

omega = linspace(-2*pi*f,2*pi*f,1024);
phi =  1i*2*pi*rand(M,1024);

Filter = H*W';
u = zeros(M,1024);
for i = 1:M
    temp = 0;
    for j = 1:M
        temp = temp + subs(Filter(i,j),omega).*exp(phi(j,:));
    end
    u(i,:) = ifft(double(temp));
end

y = lsim(G,u,linspace(0,20,1024));
