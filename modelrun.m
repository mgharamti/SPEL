function XR = modelrun( param,Xn,F,Nt,varargin )

    % Sinmulate reference states from t1 to t2.
    XR= zeros( param.Nx,Nt ) ;
    for t = 1:Nt
      % Integrate states in time.
      Xn= lorenz( param,Xn,F ) ;

      if( nargin>4 ) % Save trajectory.
          XR( :,t )= Xn ;
      end
    end
    
    if( nargin<=4 )
        XR= Xn ;
    end
 
end
