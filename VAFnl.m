function [VAF] = VAFnl(y,z)

y=y(:);
z=z(:);
y=y-mean(y);
z=z-mean(z);

VAF=(1-(sum((y-z).^2))/sum(y.^2))*100;

if VAF<0
    VAF=0;
end
end