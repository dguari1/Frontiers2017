This folder contains all the files necesary to try the identification method proposed in 
D.L Guarin and R.E. Kearney 'Estimation of Time-Varying, Intrinsic and Reflex Dynamic Joint Stiffness during Movement. 
Application to the Ankle Joint'

The folder includes a simulated data set consisting of 1000 cycles of joint position and torque data. 

The included funcions are:

- example_application_PC_id.m -> Code used to estimate the intrinsic and reflex dynamic stiffness

- np_TV_ident.m -> Implementation of proposed TV-IRF identification algorithm

- Hammer_TV_rivbj_2ndorder_ens.m -> Identification of time-varying, Hammerstein systems. Linear dynamic element has to be of second order

- TV_Bayes.m  -> Linear identification algorithm used for parameter estimation



Other tools included for completeness:
- arg_parse_c.m  -> needed to read to optional input arguments
- vec.m -> vectorize a matrix
- VAFnl.m -> compute VAF
- generate_B_splines.m -> Generate B-Splines
- multi-tcheb.m -> Generate Tchebichev polynomials
- ddt_real.m -> compute derivative 


This code is property of Diego L. Guarin, please email diego.guarinlopez at mail.mcgill.ca if you requiere further information.

