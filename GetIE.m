function [ Xa,Fa ] = GetIE( param,Xm )

    var0   = param.Vx*ones( param.Nx,1 );
    sigma0 = sqrt(var0);
    var1   = param.Vp*param.Pn;
    sigma1 = sqrt(var1);
    
    Fm = param.Fr + (0.10*param.Fr*param.Pn).*randn(1,param.Np);

    Xa = zeros( param.Nx,param.Ne ) ;
    Fa = zeros( param.Np,param.Ne ) ;
    for e = 1:param.Ne
        Xa( :,e ) = Xm' + sigma0 .* randn( param.Nx,1 ) ;
        Fa( :,e ) = Fm' + sigma1 .* randn( param.Np,1 ) ;
    end
    
end