clear;
rand('state',0);
cc=parcluster('local');
cc.NumWorkers=12;
parpool(cc,cc.NumWorkers);
name='ALLAML_uni';
m=20;
iter2=1;
load([name '.mat']);
c=length(unique(y));
n=length(y);
for iter=1:10
    SS=[];
    for j=1:m
        idxi=(iter-1)*20+j;
        YY{j}=sparse(Yi{idxi});
        SS=[SS,YY{j}];
    end
    A=SS'*SS;
    for gam=1:6
        for mu=1:7
            gamma=10^(gam-6);
            mu1=10^(mu-4);
            [iter,gam,mu]
            [~,ypred]  = SCCABG( SS,A,c,gamma,mu1 );
            ypred = ypred(1:n)';
            res=ClusteringMeasure(y,ypred)
            result_acc(iter,gam,mu)=res(1);
            result_nmi(iter,gam,mu)=res(2);
        end
    end
end