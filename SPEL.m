%
% SCRIPT_NAME: 
%   SPEL: State and Parameters Estimation Library
%
% Syntax:  
%   SPEL.m (user defined input)
%
% Inputs:
%   input1 - model:         (1) name, (2) parameters, (3) time step, (4)
%                           simulation period, (5) model error
%   input2 - assimilation:  (1) obs. variables, (2) frequency of obs., (3)
%                           ensembl size, (4) initial state var., (5)
%                           initial parameter var., (6) no. of experiments,
%                           (7) algorithm, (8) obs. error, (9)
%                           localization length scale, (10) inflation 
%                           factor, (11) index for state variables to 
%                           diagnose, (12) parameter index to diagnose
%                           
% Outputs:
%   output1 - Structure Array: Output( 1,2,...,Nexp ).{var1, var2, ...}
%               var01: Reference Trajectory for state
%               var02: Reference Parameters  
%               var03: Mean-Squared-Errors for the state
%               var04: Ensemble-Spread for the state
%               var05: Mean-Squared-Errors for the parameters
%               var06: Ensemble-Spread for the parameters
%               var07: Prediction-Skill for the state
%               var08: Prediction-Skill for the parameters
%               var09: Filter-Inbreeding-Factor for the state
%               var10: Ratio measure for state ensemble verification
%               var11: Skewness (3rd order moment)
%               var12: Kurtosis (4th order moment)
%               var13: diagnosis for state variable
%               var14: diagnosis for parameter
%
% Example: 
%   An example is readily available if you skip the input step; i.e., press
%   enter for all asked questions. 
%
% Major required m-files: 
%   1. enkf.m
%   2. Assimilate.m
%   3. Getie.m
%   4. Hmatrix
%   5. obs_operator
%
% Author:       Moha Gharamti
% Work address: Nansen Environmental and Remote Sensing Center
% Email:        mohamad.gharamti@nersc.no | new: gharamti@ucar.edu
% Website:      https://www.nersc.no | https://www.ucar.edu



fprintf( '\n============= ************* ============= \n')
fprintf(   '              SPEL STARTING               \n')
fprintf(   '============= ************* ============= \n\n')

GetModel= 'Option[1]-Lorenz63, Option[2]-Lorenz96? [1/2]: ';
model.no= input(GetModel); if isempty(model.no), model.no= 2; end

GetParm1= 'How many parameters you would like to estimate? [1/2/3]: ';
GetParm2= 'The paramter to estimate: F (atmospheric forcing term)';


switch model.no
    
    case 1
        DeclareModel= 'Selected Model: Lorenz63 (3 state variables)';
        fprintf( '%s\n\n',DeclareModel )
        
        model.Np= 3;
        model.Nx= 3; 
        model.Pn= zeros(model.Np,1);
                
        ParmNuml= input(GetParm1); if isempty(ParmNuml), ParmNuml= model.Np; end
        ParmName= cell(model.Np,1);
        if ( ParmNuml<3 )
            for iPAR = 1:ParmNuml
                ParmName{iPAR}= strcat( input( [ 'Parameter #' num2str(iPAR) ' is? [sigma/rho/beta]: ' ],'s' ) );
                if ( strcmp(ParmName{iPAR},'sigma')==0 && strcmp(ParmName{iPAR},'rho')==0 && strcmp(ParmName{iPAR},'beta')==0 )
                    error('Input parameter must be either: sigma, rho, or beta')
                else
                    parId= find( strcmp( ParmName{iPAR},{'sigma';'rho';'beta'} )~=0 );
                    model.Pn(parId)= 1;
                end
            end
        else
            model.Pn= ones(model.Np,1);
        end
        
    case 2
        DeclareModel= 'Selected Model: Lorenz96 (40 state variables)';
        fprintf( '%s\n',DeclareModel )
        fprintf( '%s\n\n',GetParm2 )
        
        model.Np= 1;
        model.Pn= 1;
        model.Nx= 40;
        
end


ModelStep= '- Model step in sec is: ';
SimPeriod= '- Total simulation period in sec: ';
ModelEror= '- Gaussian model error with mean zero and variance? [0.0001/0.001/0.01/0.1/1/...]: ';
ExpRepeat= '- Repeat the experiment how many times? [1/2/5/10/20/30/...]: ';

EnsSize  = '- Ensemble size? [10/20/50/100/...]: ';
InitVarSt= '- Variance of initial state ensemble? [1/2/3/...]: ';
InitVarPr= '- Variance of initial parameter ensemble? [0.1/0.5/1/...]: ';
ObsFreqcy= '- Observation avaiable every how many time steps? [1/2/3/4/...]: ';
ObsSparce= '- Observe which state variables? [all/half/quarter]: ';
Inflation= '- Inflation factor for State (1) and Parameters (2)? { [ 1.01,1.01 ]/[ 1.02,1.02 ]/... }: ';
Localize = '- Localization length scale? [1/5/10/20/...]: ';

DiagnosSt= '- State variable index to diagnose? [1/2/...]: ';
DiagnosPr= '- Parameter index to diagnose? [1/2/...]: ';
FilterSch= '- Ensemble filtering scheme? [EnKF/EnSRF]: ';
Algorithm= '- State-Parameters estimation algorithm? [joint/dual/conf-step/joint-osa/dual-osa]: ';
UpdateStr= '- All-at-once or serial assimilation of observations? [batch,ser]: ';


model.dt= input(ModelStep);      
if isempty(model.dt) 
    if (model.no==1), model.dt= 0.02; end
    if (model.no==2), model.dt= 0.05; end
end
model.t2= input(SimPeriod);      
if isempty(model.t2) 
    if (model.no==1), model.t2= 30; end
    if (model.no==2), model.t2= 73; end
end
model.me= input(ModelEror);      
if isempty(model.me)
    if (model.no==1), model.me= 0.001; model.na= 'L63'; end
    if (model.no==2), model.me= 0.001; model.na= 'L96'; end
end
model.Ns= input(ExpRepeat);      if isempty(model.Ns), model.Ns= 1;        end
model.Ne= input(EnsSize);        if isempty(model.Ne), model.Ne= 100;      end
model.Vx= input(InitVarSt);      if isempty(model.Vx), model.Vx= 3.0;      end
model.Vp= input(InitVarPr);      if isempty(model.Vp), model.Vp= 0.1;      end
model.do= input(ObsFreqcy);      if isempty(model.do), model.do= 4;        end
model.ob= input(ObsSparce,'s');  if isempty(model.ob), model.ob= 'half';   end
model.ga= input(Inflation);      
if isempty(model.ga)
    if (model.no==1), model.ga= [ 1.03,1.005 ]; end
    if (model.no==2), model.ga= [ 1.05,1.05 ];  end
end
model.lc= input(Localize);       
if isempty(model.lc)
    if (model.no==1), model.lc= 2; end           
    if (model.no==2), model.lc= 9; end  
end
model.fs= input(FilterSch,'s');  if isempty(model.fs), model.fs= 'EnKF' ; end
model.aa= input(Algorithm,'s');  if isempty(model.aa), model.aa= 'joint'; end
model.us= input(UpdateStr,'s');  if isempty(model.us), model.us= 'batch'; end
model.Ix= input(DiagnosSt);      
if isempty(model.Ix)
    if (model.no==1), model.Ix= 2;  end
    if (model.no==2), model.Ix= 20; end
end
model.Ip= input(DiagnosPr);      
if isempty(model.Ip)
    if (model.no==1), model.Ip= 2; end
    if (model.no==2), model.Ip= 1; end
end


fprintf( '\n\n============= ************** ============= \n')
fprintf(     '              :ASSIMILATION:               \n')
fprintf(     '============= ************** ============= \n\n')

Estimates= enkf(model);

fprintf( '\n\n========= ********************** ========= \n')
fprintf(     '          :OUTPUT VISUALIZATION:           \n')
fprintf(     '========= ********************** ========= \n\n')

VisualOutput(model,Estimates);
