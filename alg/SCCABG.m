function [ S,y ] = SCCABG( S,A_ini,c,gamma,mu )
% S: n*k Initial bipartite graph
% A_ini: k*k Initial similarity matrix of clusters
% c: Number of clusters
% gamma, mu: hyperparameters
maxiter=20;
[n,~]=size(S);
SS=S;
lambda=1;
rho=0.5;
[~,v]=size(S);
options = optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off');
W=zeros(n,v);
Stmp=zeros(n,v);
A=A_ini;
for iter=1:maxiter
    if min(min(W))<0.99
    parfor i=1:n
        C=-2.*S(i,:)'*S(i,:);
        xx=S(i,:).^2;
        C=C+repmat(xx,v,1);
        C=C+repmat(xx',1,v);
        C=C.*A;
        C=diag((S(i,:)-SS(i,:)).^2)+gamma.*C;
        C=C.*2;
        f=-rho.*ones(v,1);
        C=(C+C')./2;
        x=quadprog(C,f,[],[],[],[],zeros(v,1),ones(v,1),W(i,:),options);
        W(i,:)=x';
    end
    end
    
    d1 = sum(S,2);
    D1 = spdiags(1./max(sqrt(d1),eps),0,n,n);
    d2 = sum(S,1);
    D2 = spdiags(1./max(sqrt(d2'),eps),0,v,v);
    SS1 = D1*S*D2;
    SS2 = SS1'*SS1; 
    SS2 = full(SS2); 
    [V, ev0, ev]=eig1(SS2,c);
    U=(SS1*V)./(ones(n,1)*sqrt(ev0'));
    U = sqrt(2)/2*U; V = sqrt(2)/2*V;
    U_old = U;
    V_old = V;
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 < c-0.0000001
        lambda = 2*lambda;
    elseif fn2 > c+1-0.0000001
        lambda = lambda/2;   U = U_old; V = V_old;
    else
        break;
    end
    
    U1 = D1*U;
    V1 = D2*V;
    dist = L2_distance_1(U1',V1');

    FF=-2.*W.^2.*SS+lambda.*dist;
    parfor i=1:n
        WW=W(i,:)'*W(i,:).*A;
        H=-2.*gamma.*WW;
        H=H+diag(W(i,:).^2);
        H=H+2.*gamma.*diag(sum(WW,2));
        H=H.*2;
        ff=FF(i,:);
        H=(H+H')./2;
        x=quadprog(H,ff',[],[],[],[],zeros(v,1),ones(v,1),S(i,:),options);
        Stmp(i,:)=x';
        

    end
    S=Stmp;
    tmp=S.*W;
    tmp1=S.*tmp;
    tmp2=tmp1'*W;
    tmp3=W'*tmp1;
    tmp4=tmp'*tmp;
    tmp5=tmp3+tmp2-2*tmp4;
    A=A_ini-mu.*tmp5;
    A(A<eps)=0;
    
    rho=rho*2;
    
end
S(S<0.0001)=0;
SS0=sparse(n+v,n+v); SS0(1:n,n+1:end)=S; SS0(n+1:end,1:n)=S';
[clusternum, y]=graphconncomp(SS0);


end

