function [ YY,n_sub,k_sub] = combine( Y,m )
%COMBINE 此处显示有关此函数的摘要
%   此处显示详细说明
S=Y*Y';
[ii,jj]=find(S==m);
[n,k]=size(Y);
n_sub=[];
j=1;
processed=[];
anchor=[];
for i=1:length(jj)
    if ~ismember(jj(i),processed)
        idx=ii(find(jj==jj(i)));
        n_sub{j}=idx';
        anchor(j)=idx(1);
        processed=[processed, idx'];
        j=j+1;
    end
end
YY=Y(anchor,:);

SS=Y'*Y;
S2=bsxfun(@minus,SS,diag(SS));
S2=bsxfun(@minus,S2,diag(SS)');
S2=S2+SS;
[ii2,jj2]=find(S2==0);
k_sub=[];
j=1;
processed=[];
anchor2=[];
for i=1:length(jj2)
    if ~ismember(jj2(i),processed)
        idx=ii2(find(jj2==jj2(i)));
        k_sub{j}=idx';
        anchor2(j)=idx(1);
        processed=[processed, idx'];
        j=j+1;
    end
end
YY=YY(:,anchor2);
end

