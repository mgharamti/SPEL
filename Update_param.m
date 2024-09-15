function Fa = Update_param( param,Xf,Ff,PXf,PFf,CL2,H,Yp,R )


    Fa= zeros( param.Np,param.Ne );
    z = zeros( 1,param.Ne ); 
    
    
    if ( strcmp(param.us,'batch') == 1 )
        
        [ U,S,V ]= svd( CL2 .* ( H*PXf*H' ) + R ) ;
        Ti = find( cumsum( diag(S)/sum(S(:)) ) > 0.98,1 );
        Im = V( :,1:Ti )*diag( 1./diag( S(1:Ti,1:Ti) ) )*U( :,1:Ti )';
        
        Gain = ( PFf*H' ) * Im ;
        for e = 1:param.Ne
            Fa( :,e )= Ff( :,e ) + Gain * ( Yp( :,e ) - H*Xf( :,e ) ) ;
        end
        
        
    elseif ( strcmp(param.us,'ser') == 1 )
        
        Fa= Ff; Fam= mean( Fa,2 );
        Xa= Xf; Xam= mean( Xa,2 );
        
        for ip= 1:param.No
            
            zm= H(ip,:)*Xam;
            for ie= 1:param.Ne
                z(ie)= H(ip,:)*Xa(:,ie);
            end
            
            sumt= 0; sumo= 0;
            for ie= 1:param.Ne
                sumt= sumt+( Fa(:,ie)-Fam )*( z(ie)-zm );
                sumo= sumo+( z(ie)-zm )^2;
            end
            sumo= sumo+(param.Ne-1)*R(ip,ip);
            for ie= 1:param.Ne
                Fa(:,ie)= Fa( :,ie ) + sumt/sumo*( Yp(ip,ie)-z(ie) ); 
            end
            Fam= mean( Fa,2 );
            
        end
        
        
    else 
        warning('Update_param: "Update Scheme" is not defined, analysis is ignored.')
        
    end

    
end