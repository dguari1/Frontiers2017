function B = generate_B_splines(q,centers,sd)

q=q(:);
number=length(centers);
B=zeros(length(q),number);
for i=1:number
    B(:,i)=(1/(2*pi*sd)^(1/2))*exp(-(0.5/(sd^2))*((q-centers(i)).*(q-centers(i))));
end

return