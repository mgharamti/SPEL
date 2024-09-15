
function X1 = lorenz(param,X0,F)

    predictor= str2func(param.fm);
    
    K1 = predictor(param.Nx,X0                ,F);
    K2 = predictor(param.Nx,X0+0.5*param.dt*K1,F);
    K3 = predictor(param.Nx,X0+0.5*param.dt*K2,F);
    K4 = predictor(param.Nx,X0+1.0*param.dt*K3,F);

    X1 = X0 + (param.dt/6) * (K1 + 2*K2 + 2*K3 + K4);

end

