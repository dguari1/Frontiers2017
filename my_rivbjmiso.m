function [params_est,Ps,x,error,sigmasq,k] = my_rivbjmiso(y,u,na,nb,init_params,numiter,tol)
%Implemented by Diego Guarin
%v1.0 november 25/2012
%v1.1 november 27/2012
% log - Modified to handle multiple inputs with the same denominator
% y(k)=(B1/A)*u1(k)+...(Bm/A)*um(k) + (C/D)e(k)
% where e(k) is i.i.d noise. 

if nargin<4
    disp('not enough input arguments');
    return
end

y=y(:);
%check if there are several inputs
[N,n]=size(u);
if n>N
    u=u';
    [N,n]=size(u);
end
%N -> number of elements
%n -> number of inputs
if length(nb)~=n
    disp('The size of nb does not match with the number of inputs');
    return
end


if nargin<7
    tol=0.000001;
end
if nargin<6
    numiter=20;
end
if nargin<5
    X=regress_matrix (y,u,na,nb);
    params_est=X\y;
    params_est=params_est(:);
else
    if isempty(init_params)
            X=regress_matrix (y,u,na,nb);
            params_est=X\y;
            params_est=params_est(:);
    else
        params_est=init_params(:);
    end
end


A=zeros(na+1,1);
B=cell(n,1);
A=[1; params_est(1:na)];
aux=0;
for p=1:n
    B{p}=params_est(na+1+aux:na+aux+nb(p));
    aux=aux+nb(p);
end

for k=1:numiter

    params_est_old=params_est;
    A=polystab(A);
%    Aold=A';
%    Bold=B;
    temp=0;
    for i=1:n
        temp=temp+filter(B{i},A,u(:,i));
        u1(:,i)=filter(1,A,u(:,i));
    end
    x=temp; %instrumental variable
    x1=filter(1,A,x);
    y1=filter(1,A,y);

    
    Regress_X=regress_matrix(x1,u1,na,nb);
    Regress_Y=regress_matrix(y1,u1,na,nb);
    
    params_est=(Regress_X'*Regress_Y)\(Regress_X'*y1);
    
    A=[1; params_est(1:na)];
    aux=0;
    for p=1:n
        B{p}=params_est(na+1+aux:na+aux+nb(p));
        aux=aux+nb(p);
    end
    
    if k>1
       if max(abs(params_est-params_est_old))<tol
           break
       end
    end
end

temp=0;
for i=1:n
    temp=temp+filter(B{i},A,u(:,i));
    u1(:,i)=filter(1,A,u(:,i));
end
x=temp; %instrumental variable 
x1=filter(1,A,x);
y1=filter(1,A,y);

Regress_Y=regress_matrix(y1,u1,na,nb);

Regress_X=regress_matrix(x1,u1,na,nb);


error=y-x;
sigmasq=cov(error(na+1:end));
Ps=sigmasq*(eye(na+sum(nb))/(Regress_X'*Regress_X));
%Pr=sigmasq*(eye(na+sum(nb))/(Regress_X'*Regress_Y));
end



function X = regress_matrix (y,u,na,nb)
%model 
%y(k)=-a(1)y(k-1)-...-a(na)y(k-na)+b(0)u(k)+...b(nb)u(k-nb)+e(k)

y=y(:); %turnign y into a column vector
n=length(nb);
Y=[];
U=[];
if na>0
    Y=y;
    for j=2:na+1
        Y(:,j)=[0;Y(1:end-1,j-1)];
    end
    Y=-Y; %I want to preserve the sign of the coefficients
    Y(:,1)=[];
end
for p=1:n
    Utemp=[];
    if nb(p)>0
        Utemp=u(:,p);
        for j=2:nb(p)
            Utemp(:,j)=[0;Utemp(1:end-1,j-1)];
        end
    end
    U=[U Utemp];
end
X=[Y U]; %regressor matrix
end