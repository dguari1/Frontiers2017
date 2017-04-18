function [x_est,u_bar_old,a_est,b_est,TV_Nonlinearity,params_lin,params_nl,Basis]=Hammer_TV_rivbj_2ndorder_ens (y,u,M,Basis,init_cond)
%Identification of Hammerstein systems with Time Varying linear element and
%Time Varying static nonlinearity
%
%Method was introduced in: An Instrumental Variable Approach for the 
%Identification of Time-Varying, Hammerstein Systems 
%Authors: Diego L. Guarin  Robert E. Kearney 
%Preprints of the 17th IFAC Symposium on System Identification
%
%Continuous-time model of the form 
% u(t) --> g(u(t),t) -->  u_bar(t)
% d^2x(t)/dt^2 + 2z(t)w(t)dx(t)/dt + w^2x(t)=G(t)w(t)^2u_bar(t)
% y(t)=x(t)+e(t)


[N,~]=size(y);  %Number of data points 

na=M(1);
nb=M(2);
nn=M(3);

%checking the basis
if  length(Basis)~=3
        error('Not enouth basis functions')
end
if length(Basis{1})~=(na)
        error('Each parameter requieres a set of basis functions')
end
if length(Basis{2})~=(nb)
        error('Each parameter requieres a set of basis functions')
end
if length(Basis{3})~=(nn+1)
        error('Each parameter requieres a set of basis functions')
end
    
    %some of the basis must be vectors of 1. 
for i=1:na
    if sum(Basis{1}{i}(:,1))~=N
            error('The first basis function used to represent the elements of the denominator have to be 1');
    end
end    
%How many elements are used to represente each element?
for i=1:na
        [~,N_basis{1}(i)]=size(Basis{1}{i});
end
for i=1:nb
        [~,N_basis{2}(i)]=size(Basis{2}{i});
end    
for i=1:nn+1
        [~,N_basis{3}(i)]=size(Basis{3}{i});
end     


if isempty(init_cond)
    disp('Initial conditions not provided, will be estimated by considering the system as linear, TV (this might produce errors).')
end

[x_est,u_bar_old,a_est,b_est,TV_Nonlinearity,params_lin,params_nl,iter_alg,iter_riv,mse] = Main_Iteration(y,u,M,Basis,N_basis,init_cond);

end

function [x_est,u_bar_old,a_est,b_est,TV_Nonlinearity,params_lin,params_nl,iter_alg,iter_riv,mse] = Main_Iteration(y,u,M,Basis,N_basis,init_cond)

    [N, trials]=size(y);
    
    na=M(1);
    nb=M(2);
    nn=M(3);
    Max_iter=50;
    Max_iter_riv=50;

    %Computing initial conditions of these are not provided
    if isempty(init_cond)
        
            %new input formed as n_u(k)=u(k)+2*u(k-1)+u(k-2)  to reduce complexity
            %of discrete model. This only works with second order models,
            %for higher order it should be modified
            clear INN n_u
            INN=regress_matrix(zeros(N,1),vec(u),0,3);
            n_u=INN(:,1)+2*INN(:,2)+INN(:,3);
            [params_est,~,~] = my_rivbjmiso(vec(y),n_u,na,nb);
            b_est=params_est(na+1:end)';
            b_est=repmat(b_est,N,1);
            F=polystab([1 params_est(1:na)']);
            a_est=F(2:end);
            a_est=repmat(a_est,N,1);
    else %if initial conditions are provided we can use them
        F=[1 init_cond(1:2)'];
        aux=na+1;
        for i=1:na
            index_den=[i aux:(aux-1)+(N_basis{1}(i)-1)];
            aux=aux+N_basis{1}(i)-1;
                a_est(:,i)=Basis{1}{i}*init_cond(index_den);
                clear index_den
        end
        for i=1:nb
            index_num=aux:(aux-1)+N_basis{2}(i);
            aux=aux+N_basis{2}(i);
            b_est(:,i)=(Basis{2}{i}*init_cond(index_num));
            clear index_num
        end
    end
        
    %Approximating the shape of the static-nonlinearity with an expansion
    %of Chebyshev polynomials 
    
    if nn>0  %if th
    
        %The nonlinearity is represented by a expanding the input using
        %Chebyshev polynomials 
        avg = (max(vec(u))+min(vec(u)))/2;
        rng = max(vec(u)) - min(vec(u));
        un = (vec(u) - avg)*2/rng;
        U = multi_tcheb(un,nn);
        
        %is possible to use other expansions but these are not yet
        %implemented. For example, radial basis functions ...
        % U = radial_basis_func(vec(u),nn+1);
        
        
    
        %U contains the basis function expansion of the input
        %Applying the basis function to U to create a new set of inputs,
        %see equation 21 in the paper
        NLU_all=[];
        aux=1;
        for p=1:trials
            clear U_temp 
            NLU=[];
            U_temp=U(aux:N+aux-1,:);
            for i=1:nn+1
                NLU=[NLU repmat(U_temp(:,i),1,N_basis{3}(i)).*Basis{3}{i}];
            end
            NLU_all=[NLU_all;NLU];
            aux=aux+N;
        end
    end
        
    clear mse
    for iter_alg=1:Max_iter   
        %%%%%%%%%%%%--------------------%%%%%%%%%%%%%
        %here we estimate the static-nonlinearity, it is assumed that we
        %have some information about the linear component (either estimated
        %or provided by the initial conditions)
        if nn>0
            %now we filter the new inputs with the current estimate of the
            %linear component. See equation 21 in the paper
            %tic
            [total_size,elem_nol_size]=size(NLU_all);
            NLU_fil=zeros(total_size,elem_nol_size);
            for k=1:elem_nol_size
                clear TEMP_matrix TEMP_matrix2 INN new_NLU_all
                TEMP_matrix=reshape(NLU_all(:,k),[],trials);
                TEMP_matrix2=evaluate_filter_TV(TEMP_matrix,b_est,a_est);
                NLU_fil(:,k)=vec(TEMP_matrix2);
            end
            %toc

            %Now we estimate the weights associates to the
            %static-nonlinearity, see equation 22 in the paper
            [~,params_nl,noise_var,lk,SIGMA]=TV_Bayes(vec(y),[],NLU_fil);

            sign_carrier=sign(params_nl(1));
            norm_carrier=norm(params_nl);
            params_nl=(params_nl.*sign_carrier)./norm_carrier;   
    
            %an estimate the intermediate signal, see equation 16 in the
            %paper
            u_bar_old=NLU_all*params_nl;
            clear INN
            INN=regress_matrix(zeros(N,1),u_bar_old,0,3);
            u_bar=INN(:,1)+2*INN(:,2)+INN(:,3);
        
        else
            u_bar_old=vec(u);
            clear INN
            INN=regress_matrix(zeros(N,1),u_bar_old,0,3);
            u_bar=INN(:,1)+2*INN(:,2)+INN(:,3);
            sign_carrier=1;
            norm_carrier=1;
        end
        
        %%%%%%%%%%%%--------------------%%%%%%%%%%%%%
        %Here we estimate the linear component. At this point we already
        %computed the intermediate signal (the input to the linear
        %element), that signal and the output can be used  to estimate
        %the linear system. We do that using a TV-instrumental variables
        %approach. 
        u_bar_trials=reshape(u_bar,[],trials);

        b_est=b_est.*(sign_carrier.*norm_carrier);
      
        for iter_riv=1:Max_iter_riv    
            
            x_est=evaluate_filter_TV(u_bar_trials,b_est,a_est);

            XQX_f=[];
            XQY_f=[];
            Y_en=regress_matrix(vec(y),u_bar,na,0);
            X_en=regress_matrix(vec(x_est),u_bar,na,0);
            U_bar_en=regress_matrix(vec(x_est),u_bar,0,nb);
            for l=1:trials     
                Y=Y_en(N*(l-1)+1:N*(l),:);
                X=X_en(N*(l-1)+1:N*(l),:);
                U_bar=U_bar_en(N*(l-1)+1:N*(l),:);

                %this is where the Basis functions come in. Each colum of 
                %the regressor matrices is multiplied by the basis function 
                XQ=[];
                for i=1:na
                        %Remember that the first Basis function is 1. 
                        %So it does not affect the matrix of regressors
                        XQ=[XQ repmat(X(:,i),...
                            1,N_basis{1}(i)-1).*Basis{1}{i}(:,2:end)];  
                end
                UQ=[];
                for i=1:nb
                        UQ=[UQ ...
                          repmat(U_bar(:,i),1,N_basis{2}(i)).*Basis{2}{i}];
                end
    
                %creating the full regressor matrices, one for the 
                %noise-free estimated output and another for the measured 
                %output
                XQY=[Y XQ UQ];
                XQX=[X XQ UQ];
    
                %we now have to pre-filter the data, we are as the data
                %are periodic, we are using the final conditions of one
                %cycle as the initial conditions of the next one. This
                %should be corrected if the data are not periodic. 
                if l==1
                    [temp_xqx,init_xqx]=filter(1,F,XQX);
                    [temp_xqy,init_xqy]=filter(1,F,XQY);
                    [temp_y,init_y]=filter(1,F,y(:,l));
                else
                    [temp_xqx,init_xqx]=filter(1,F,XQX,init_xqx);
                    [temp_xqy,init_xqy]=filter(1,F,XQY,init_xqy);
                    [temp_y,init_y]=filter(1,F,y(:,l),init_y);
                end
                %filtering all this with the filter 1/F
                XQX_f=[XQX_f;temp_xqx];
                XQY_f=[XQY_f;temp_xqy];
                %filtering the measured output with the filter 1/F
                y_f(:,l)=temp_y;
            end
            %solving the least-squares problem with instrumental variables
            %to compute the parameters
            params_lin=(XQX_f'*XQY_f)\(XQX_f'*vec(y_f));
    
            %updating the filter F, A and B
            F=[1 params_lin(1:na)'];
            F=polystab(F);  %if F is not stable then force it to be. 
            

            old_a=a_est;
            old_b=b_est;
            %estimating the time-varying filter using the estimated 
            %parameters and the basis functions
    
            %index_den=zeros(na,N_basis);
            mse_a=zeros(na,1);
            aux=na+1;
            for i=1:na
                index_den=[i aux:(aux-1)+(N_basis{1}(i)-1)];
                aux=aux+N_basis{1}(i)-1;
                a_est(:,i)=Basis{1}{i}*params_lin(index_den);
                clear index_den
                %checking if there is a significant large between the 
                %old and new paramerer
                mse_a(i,1)=(1/N)*sum((a_est(:,i)-old_a(:,i)).^2);
            end
            %index_num=zeros(nb,N_basis);
            mse_b=zeros(nb,1);
            for i=1:nb
                index_num=aux:(aux-1)+N_basis{2}(i);
                aux=aux+N_basis{2}(i);
                b_est(:,i)=(Basis{2}{i}*params_lin(index_num));%*(sign_carrier*norm_lambda);
                clear index_num
                %checking if there is a large difference between the 
                %old and new paramerer
                mse_b(i,1)=(1/N)*sum((b_est(:,i)-old_b(:,i)).^2);
            end
    
            %checking if there was any significant change in the TV parameters
            %this is assuming that all the parameters are reaching its limit value.
            if max(max(mse_a),max(mse_b))<1e-3
                break;
            end
   

        end
        %computing the difference between the measured and predicted
        %output, the algorithm finalizes if this difference doesn't change
        %much between iterations
        mse(iter_alg,1)=(1/(N*trials))*sum((vec(y)-vec(x_est)).^2);
        if iter_alg>1
            relative_diff=(abs((mse(iter_alg,1)-mse(iter_alg-1,1))/mse(iter_alg-1,1)))*100;
            if relative_diff<0.5
                break;
            end
        end
  
    end
    
    if iter_alg==Max_iter
        disp('The algorithm did not converge. Increase the maximum number of iterations or modify the model')
    end
    if nn>0       
        %computing the non-linearity one last time with the final values of
        %b_est and a_est
        [total_size,elem_nol_size]=size(NLU_all);
        NLU_fil=zeros(total_size,elem_nol_size);
        for k=1:elem_nol_size
            clear TEMP_matrix TEMP_matrix2 INN new_NLU_all
            INN=regress_matrix(zeros(N,1),NLU_all(:,k),0,3);
            new_NLU_all=INN(:,1)+2*INN(:,2)+INN(:,3);
            TEMP_matrix=reshape(new_NLU_all,[],trials);
            TEMP_matrix2=evaluate_filter_TV(TEMP_matrix,b_est,a_est);
            NLU_fil(:,k)=vec(TEMP_matrix2);
        end
        %toc

        %params_nl=regress(y,NLU_fil);
        [~,params_nl,noise_var,lk,SIGMA]=TV_Bayes(vec(y),[],NLU_fil);

        sign_carrier=sign(params_nl(1));
        norm_carrier=norm(params_nl);
        params_nl=(params_nl.*sign_carrier)./norm_carrier;   
    

        u_bar_old=NLU_all*params_nl;
        clear INN
        INN=regress_matrix(zeros(N,1),u_bar_old,0,3);
        u_bar=INN(:,1)+2*INN(:,2)+INN(:,3);
        
        
        %normalizing the elements of the nonlinearity 
        aux=1;
        for i=1:nn+1
            index_params=aux:aux+(N_basis{3}(i)-1);        
            TV_Nonlinearity(:,i)=(Basis{3}{i}(:,:)*params_nl(index_params));%*(sign_carrier/norm_lambda);
            aux=aux+N_basis{3}(i);                                    
        end 
        
    else
        params_nl=[];
        TV_Nonlinearity=[];
        u_bar_old=vec(u);
        clear INN
        INN=regress_matrix(zeros(N,1),u_bar_old,0,3);
        u_bar=INN(:,1)+2*INN(:,2)+INN(:,3);
    end
    
    u_bar_trials=reshape(u_bar,[],trials);
    b_est=b_est*sign_carrier*norm_carrier;

    x_est=evaluate_filter_TV(u_bar_trials,b_est ,a_est);

end

function out = evaluate_filter_TV(u,B,A)

% N=length(u);
% out=zeros(N,1);
% for i=1:N
%     Atemp=[1 A(i,:)];
%     Atemp=polystab(Atemp);  %make sure that the filter at time i is stable
%     temp=filter(B(i,:),Atemp,u);  %filter the input with the current filter.
%     out(i,1)=temp(i); %extrac the output at the current time.
%     clear temp Atemp
% end
[N, N_trials]=size(u);

out=zeros(size(u));
for i=1:N
    A_e(i,:)=[1 A(i,:)];
    A_e(i,:)=polystab(A_e(i,:));  %make sure that the filter at time i is stable
end  
A_e(:,1)=[];
for n=1:N_trials
    if n==1
        out(1,n)=B(1)*u(1,n);
        out(2,n)=B(2)*u(2,n)-A_e(2,1)*out(1,n);
        for i=3:N
            out(i,n)=B(i)*u(i,n)-A_e(i,1)*out(i-1,n)-A_e(i,2)*out(i-2,n);
        end
    else
        out(1,n)=B(1)*u(1,n)-A_e(1,1)*out(end,n-1)-A_e(1,2)*out(end-1,n-1);
        out(2,n)=B(2)*u(2,n)-A_e(2,1)*out(1,n)-A_e(2,2)*out(end,n-1);
        for i=3:N
            out(i,n)=B(i)*u(i,n)-A_e(i,1)*out(i-1,n)-A_e(i,2)*out(i-2,n);
        end
    end
end

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
        %Y(:,j)=[0;Y(1:end-1,j-1)];
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
            %Utemp(:,j)=[0;Utemp(1:end-1,j-1)];
        end
    end
    U=[U Utemp];
end
X=[Y U]; %regressor matrix
end



function [W, flags] = multi_tcheb(V, max_order);
%
%  usage W = multi_herm2(V, max_order);
%
%  given a collection of column vectors, V, this function returns
%  a matrix of all of the Tcheb functions up to max_order applied
%  to all of the vectors in V.

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

if max_order==0,
   W=V*0 +1;
   return;
end




[nr,nc] = size(V);


%  create a matrix of flags, such that flags (i,j) points to the 
%  column of the i'th basis function raised to the j'th power.

flags = zeros(nc,max_order);
flags(:,1) = [2:nc+1]';
for order = 2:max_order
   % first column is offsetr by 1 from last value of previous order
   flags(1,order) = flags(nc,order-1)+1;
    for i = 2 : nc
        num_terms = flags(nc,order-1)-flags(i-1,order-1)+1;
        flags(i,order) = flags(i-1,order)+ num_terms;
      end
   end

   

% pre-allocate the W matrix
%W = zeros(nr,flags(nc,max_order));
%  generate the single basis function Hermite polynomials and place 
%  them  in the correct columns of W

Te = zeros(nr,max_order+1);
Te(:,1) = ones(nr,1);

%  generate all of the functions involving a Hermite polynomial
%  applied to a single basis vector

for i = 1:nc
    Te(:,2) = V(:,i);
    for j = 2:max_order
        Te(:,j+1) = 2*V(:,i).*Te(:,j) - Te(:,j-1);
      end
    for j = 1:max_order
        W(:,flags(i,j)) = Te(:,j+1);
      end
  end
W(:,1) = ones(nr,1);

%clear Te V


%  generate all of the functions involving a powers of a single input vector


  
%  Now, using the functions that we just created, fill in the 
%  rest of the matrix



for order = 2:max_order
    for v1 = 1 : nc-1
        index = flags(v1,order);
        for v1_order = order-1:-1:1
            term1 = W(:,flags(v1,v1_order));
            rem_order = order - v1_order;
            
%           find the terms of order rem_order, whose leading variable
%           is 'greater' than v1

            first_term = flags(v1+1,rem_order);
            last_term = flags(nc,rem_order);
            for j = first_term:last_term
                index = index+1;
                W(:,index)=term1.*W(:,j);
                disp([index flags(v1,v1_order) j])
              end
          end
      end
end
end


% function U = radial_basis_func(u,number)
% 
%     centers=linspace(min(u),max(u),number);
%     sd=1*(centers(2)-centers(1));
%     U=zeros(length(u),number);
%     %figure; hold on
%     for i=1:number
%         U(:,i)=(1/(2*pi*sd)^(1/2))*exp(-(0.5/(sd^2))*((u-centers(i)).*(u-centers(i))));
%         %plot(u,BF(:,i),'.');
%     end
% 
% end