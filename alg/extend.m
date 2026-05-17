function [ y_extend ] = extend( y,n_sub,n )
%EXTEND 此处显示有关此函数的摘要
%   此处显示详细说明

y_extend=zeros(n,1);
for i=1:length(y)
    y_extend(n_sub{i})=y(i);
end

end

