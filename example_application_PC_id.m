%loading simulated data 

load('ENSEMBLE_DATA_new.mat')
%this file contains simulated noise free and noisy data. You can compare
%the prediction results vs the noise-free simulated data. 

POSITION_SIGNAL=ENSEMBLE_POSITION;
TORQUE_SIGNAL=ENSEMBLE_TORQUE;

%there are 1000 cycles, you need to define how many will be used for ID
trials=80;

%sampling rate was 1000Hz
T=1/1000;




%Extracting data
%we assume that data is periodic so it will be arrange as a single vector
%for preprosesing
POS_signal=vec(POSITION_SIGNAL(:,1:trials));
TOR_signal=vec(TORQUE_SIGNAL(:,1:trials));

%computing velocity signal
VELOCITY_signal=ddt_real(POS_signal,1/1000);


%data will be decimated to 100Hz
d_r=10;
Ts=T*d_r;
%decimation
pos_signal=decimate(POS_signal,d_r);
tor_signal=decimate(TOR_signal,d_r);
vel_signal=decimate(VELOCITY_signal,d_r);

%we will assume the reflex delay to be 40ms. This can be verified in the
%experimetal data with the provided EMG signal
delay=0.04;
tau=ceil(delay/(T*d_r)); %discrete delay
del_vel_signal=[zeros(tau,1);vel_signal(1:end-tau)];

%now we can re-arrange the decimated data into trials 
Ns=143;
pos=reshape(pos_signal,Ns,trials);
tor=reshape(tor_signal,Ns,trials);
delvel=reshape(del_vel_signal,Ns,trials);



%Generating basis functions for representation of TV-IRF
%B-splines 
%these are time-based basis functions. You need to define the numner 
%of basis to use, a time vector, a center vector and the stardar deviation
%of each spline. 
number=10;
q=linspace(0,1,Ns)';
centers=linspace(min(q),max(q),number);
sd=0.1;
B=generate_B_splines(q,centers,sd);


%approximation of reflex stiffness, 
%M = [elements_in_demonimator elements_in_numerator static_nonlinearity];
M=[2 1 5];

%Basis functions to approximate the TV, static-nonlinearity 
number=10;
q=linspace(0,1,Ns)';
centers=linspace(min(q),max(q),number);
sd=0.1;
BR_1=generate_B_splines(q,centers,sd);
N_basis=6;
BR_2=multi_tcheb(linspace(-1,1,Ns)',N_basis);

clear BASIS_NL
for i=1:M(3)+1  
        BASIS_NL{i}=ones(Ns,1);

end
%Basis functions to approximate the TV, linear component. In this example
%we will make the linear component time-invariant
clear BASIS
BASIS{1}={BR_2 BR_2}; %basis for the elements in the denominator
BASIS{2}={BR_1}; %basis for the elements in the numerator 
BASIS{3}=BASIS_NL;


%Generating basis functions for representation additional torque 
%Chebyshev polynomials
N_basis=3;
BASIS_TVMEAN=multi_tcheb(linspace(-1,1,Ns)',N_basis);


%initialization for iterative algorithm 
tqI=zeros(Ns,trials);
tqR=zeros(Ns,trials);
tqT=zeros(Ns,trials);
TVmean=zeros(Ns,trials);
VAFbest=-1000;
id_tolerance=0.01; 
max_iter=3;


disp('starts here')

for iter=1:max_iter
          
   %estimation of intrinsic component
   tqI=tor-(tqR+TVmean);
   [H_I, x_pred] = np_TV_ident(pos, tqI, B, 'nLags',4,'nSides',2,'periodic','yes','method','Bayes');
   tqI=x_pred;
   clear x_pred
    
   %estimation of reflex component 
   tqR=tor-(tqI+TVmean);
   %The results from previous iterations are used as initial conditions
   if iter==1
        [tqR,u_bar_old,a_est,b_est,TV_Nonlinearity,params_lin,params_nonlin,Basis]=Hammer_TV_rivbj_2ndorder_ens(tqR,delvel,M,BASIS,[]);
   else
        [tqR,u_bar_old,a_est,b_est,TV_Nonlinearity,params_lin,params_nonlin,Basis]=Hammer_TV_rivbj_2ndorder_ens(tqR,delvel,M,BASIS,params_lin);
   end
    
   %estimation of additional component
   TVmean=tor-(tqI+tqR);
   for i=1:trials
       [TVmean(:,i),param(:,i)]=TV_Bayes(TVmean(:,i),[],BASIS_TVMEAN);
   end

   
    

   tqT=tqI+tqR+TVmean;
   
   VAFT=VAFnl(vec(tor(:,2:end-1)),vec(tqT(:,2:end-1)));
   ss=sprintf('Iteration: %4i \tTotal VAF: %6.3f\n',iter, VAFT);
   disp(ss);
    
    if abs((VAFT-VAFbest)) < id_tolerance,
        break
    elseif VAFT<VAFbest
        break
    else
        VAFbest=VAFT;
    end
    
end

%%
%some plots
%%%%%%%%%%%%%%%-----------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Intrinsic component 
%Estimate static stiffness, viscosity and inertia (those were the ones
%simulated)
N_B=size(B,2);
BASS{1}=repmat(B,trials,1);
BASS{2}=repmat(B,trials,1);
BASS{3}=ones(Ns*trials,1);%repmat(B,trials,1);

ACCELERATION_signal=ddt_real(VELOCITY_signal,1/1000);
acceleration_signal=decimate(ACCELERATION_signal,10);

[x_est,mu,noise_var,lk,SIGMA,par,x_est1]=TV_Bayes(tqI(:),...
[pos_signal vel_signal acceleration_signal],BASS);
Estimated_K=B*mu(1:N_B);
Estimated_B=B*mu(N_B+1:2*N_B);
Estimated_I=mu(end)*ones(143,1);%B*mu(2*N_B+1:end);%*ones(Ns,1);



figure;
subplot(3,1,1)
plot((0:143-1)*(1/100),Estimated_K,'Color',[77 190 238]/255,'LineWidth',3)
title('Intrinsic static stiffness')
subplot(3,1,2)
plot((0:143-1)*(1/100),Estimated_B,'Color',[77 190 238]/255,'LineWidth',3)
title('Intrinsic viscosity')
subplot(3,1,3)
plot((0:143-1)*(1/100),Estimated_I,'Color',[77 190 238]/255,'LineWidth',3)
title('Intrinsic inertia')


%%%%%%%%%%%%%%%-----------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reflex component 
%computing the reflex linear elements: gain, natural frquency and damping

Gain=4*(b_est./(1+a_est(:,1)+a_est(:,2)));
nat_frequency=(2/Ts)*((1+a_est(:,1)+a_est(:,2))./(1-a_est(:,1)+a_est(:,2))).^(0.5);
damping=((1-a_est(:,2)))./sqrt((1+a_est(:,1)+a_est(:,2)).*(1-a_est(:,1)+a_est(:,2)));

figure;
subplot(3,1,1)
plot((0:143-1)*(1/100),Gain,'Color',[77 190 238]/255,'LineWidth',3)
title('Reflex gain')
subplot(3,1,2)
plot((0:143-1)*(1/100),nat_frequency,'Color',[77 190 238]/255,'LineWidth',3)
title('Reflex natural frequency')
subplot(3,1,3)
plot((0:143-1)*(1/100),damping,'Color',[77 190 238]/255,'LineWidth',3)
title('Damping')
%%%%%%%%%%%%%%%-----------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input and output of static-nonlinearity 
figure
subplot(2,1,1)
plot((0:Ns*trials-1)*(1/100),delvel(:),'Color',[77 190 238]/255,'LineWidth',3);
title('Joint velocity')
subplot(2,1,2)
plot((0:Ns*trials-1)*(1/100),u_bar_old(:),'Color',[77 190 238]/255,'LineWidth',3);
title('Non-linear joint velocity')

