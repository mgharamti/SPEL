function OutExp= enkf(model)

rng default


%% REFERENCE:
param.dt= model.dt ; param.Np= model.Np ;  
param.t1= param.dt ; param.t2= model.t2 ;        
param.Nx= model.Nx ; param.Ns= model.Ns ;
param.Vx= model.Vx ; param.Vp= model.Vp ;
param.Ix= model.Ix ; param.Ip= model.Ip ;
param.ob= model.ob ; param.Pn= model.Pn ;

param.Nt= length( param.t1:param.dt:param.t2 ) ;

if (model.no==1)
    param.fm= 'florenz63'; 
    param.Fr= [ 10,28,8/3 ];
    X0      = [ 20,15,30 ] ;
end
if (model.no==2)
    param.fm= 'florenz96'; 
    param.Fr= 8 ; 
    X0= param.Fr*ones( 1,param.Nx ) ; X0(20) = ( 1+0.001 )*param.Fr ;
    X0= modelrun( param,X0,param.Fr,param.Nt*0.25 );
end
XR= modelrun( param,X0,param.Fr,param.Nt,'YES' );


%% OBSERVATIONS
Yo = obs_operator( XR,param.ob ) ;
param.No = size( Yo,1 ) ;  

var = ones( param.No,1 ) ;       
sig = sqrt( var ) ;
Ro  = diag( var ) ;
[ H,Ol ] = Hmatrix( param.Nx,param.ob ) ;

Y = zeros( param.No,param.Nt ) ;
for t = 1:param.Nt
    Y( :,t ) = Yo( :,t ) + sig .* randn( param.No,1 ) ;
end


%% ASSIMILATION SETUP:
param.Ne = model.Ne; % Ensemble Size
param.do = model.do; % Observation frequency in time
param.ga = model.ga; % Inflation factor
param.lc = model.lc; % When lc is set to 30; zeros correlations (no more localization)
param.me = model.me; % Variance of gaussian model error 

param.fs = model.fs; % Ensemble Filter
param.aa = model.aa; % State-Parameters Algorithm
param.us = model.us; % Updating Strategy 

disp( [ '- Observe ' num2str( param.No ) ' variables' ] )
disp( [ '- Assimilation every ' num2str( param.do ) ' time steps' ] )
disp( [ '- Ensemble size is ' num2str( param.Ne ) ] )
disp( [ '- Model error is gaussian noise with 0 mean and variance ' num2str( param.me ) ] )
disp( [ '- Inflation factor is: State= ' num2str( param.ga(1) ) ', Parameters= ' num2str( param.ga(2) ) ] )
disp( [ '- Localization length scale is ' num2str( param.lc ) ] )
    
[ Xa,Fa ] = GetIE( param,mean( XR,2 )' ) ;

Xf = zeros( param.Nx,param.Ne ) ;
Ff = zeros( param.Np,param.Ne ) ; 

fprintf( '\n' )
OutExp(param.Ns)= struct;
for iEXP = 1:param.Ns 
    disp( [ '  > Experiment: ',num2str(iEXP) ] )
    Output= Assimilate( param,Xa,Fa,Xf,Ff,XR,H,Ol,Y,Ro,sig ) ;
    if ( iEXP==1 ), OutExp= Output; else OutExp(iEXP)= Output; end
end