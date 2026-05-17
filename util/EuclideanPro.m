function x = EuclideanPro( w,v)
%EUCLIDEANPRO 此处显示有关此函数的摘要
%   此处显示详细说明
%   min sum(w.*(x-v).^2)
%   s.t. sum(x)=1

w2=sqrt(w);
c=w2.*v;
a=1./max(w2,eps);
s=w.*v;
%s=c./max(a,eps);
maxs=max(s);
theta=maxs;
%eps=1e-9;
%idx=find(s>=theta);
idx=s>=theta;
s1=sum(idx);
ac=a(idx).*c(idx);
sumac=sum(ac);
suma2=max(sum(a(idx).^2),eps);
iter=0;
while true
    iter=iter+1;
    
        
    theta=(sumac-1)/(suma2);
    theta2=theta-(1e-9);
%    idx1=find(s>=theta-(1e-9));
    idx1=s>theta2;
    ac=a.*c;
%    ac=a(idx1).*c(idx1);
    sumac=sum(ac(idx1));
    aa=a.*a;
    sumaa=sum(aa(idx1));
    suma2=max(sumaa,eps);
    fihh=-theta*suma2+sumac-1;
    s2=sum(idx1);
    if abs(fihh)<=eps || s1==s2
        break;
    end
    if iter>100
        break;
    end
    s1=s2;
%    idx=idx1;
end
x=max(0,c-theta*a);
x=x./w2;
    

end

