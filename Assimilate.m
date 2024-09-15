function Output = Assimilate( param,Xa,Fa,Xf,Ff,XR,H,Ol,Y,R,sig )


    CL1 = Taper( param.Nx,Ol,param.lc, 'PH' ) ; 
    CL2 = Taper( param.Nx,Ol,param.lc,'HPH' ) ; 
    
    Fam = mean( Fa,2 ) ;
    Xam = mean( Xa,2 ) ;

    t = 0 ; 
    s = 0 ;
    
    while( t < param.Nt )

        t = t + param.do ;
        s = s + 1 ;

        % Propagate ensemble members to next observation time:
        for e = 1:param.Ne
            Ff( :,e )= Fa( :,e ) ;
            Xf( :,e )= modelrun( param,Xa( :,e )',Fa( :,e ),param.do )' + sqrt(param.me)*randn(param.Nx,1) ;
        end      
        Xfm= mean( Xf,2 ) ;
        Ffm= mean( Ff,2 ) ;

        % Perturb Observations for each ensemble member:
        Yp= zeros( param.No,param.Ne ) ;
        for e = 1:param.Ne
            Yp( :,e )= Y( :,t ) + sig.*randn( param.No,1 ) ;
        end

        % Inflate forecast members around the mean
        for e = 1:param.Ne
            Xf( :,e )= param.ga(1) * ( Xf( :,e ) - Xfm ) + Xfm ; 
            Ff( :,e )= param.ga(2) * ( Ff( :,e ) - Ffm ) + Ffm ; 
        end
        Xfp= Xf - repmat( Xfm,1,param.Ne ) ;
        Ffp= Ff - repmat( Ffm,1,param.Ne ) ;

        % Forecast statistics:
        PXf= 1/( param.Ne-1 ) * ( Xfp * Xfp' ) ;
        PFf= 1/( param.Ne-1 ) * ( Ffp * Xfp' ) ; 

        if ( strcmp(param.aa,'joint') == 1 )
            
            Fa = Update_param( param,Xf,Ff,PXf,PFf,CL2,H,Yp,R );        %Update param ensmeble
            [ GainX,Xa ] = Update_state( param,Xf,PXf,CL1,CL2,H,Yp,R ); %Update state ensemble
            

        elseif ( strcmp(param.aa,'dual') == 1 )

            Fa = Update_param( param,Xf,Ff,PXf,PFf,CL2,H,Yp,R );
            
            for e = 1:param.Ne
                Ff( :,e )= Fa( :,e ) ;
                Xf( :,e )= modelrun( param,Xa( :,e )',Fa( :,e ),param.do )' + sqrt(param.me)*randn(param.Nx,1) ;
            end      
            Xfm= mean( Xf,2 ) ;

            for e = 1:param.Ne
                Xf( :,e )= param.ga(1) * ( Xf( :,e ) - Xfm ) + Xfm ; 
            end
            Xfp= Xf - repmat( Xfm,1,param.Ne ) ;
            PXf= 1/( param.Ne-1 ) * ( Xfp * Xfp' ) ;

            [ GainX,Xa ] = Update_state( param,Xf,PXf,CL1,CL2,H,Yp,R );


        elseif ( strcmp(param.aa,'joint-osa') == 1 )

            % Most recent analysis ensemble:
            Xap= Xa - repmat( Xam,1,param.Ne ) ;
            Fap= Fa - repmat( Fam,1,param.Ne ) ;
            
            PFa= 1/( param.Ne-1 ) * ( Fap * Xfp' ) ;
            PXa= 1/( param.Ne-1 ) * ( Xap * Xfp' ) ;

            Fs = Update_param( param,Xf,Fa,PXf,PFa,CL2,H,Yp,R );
            [ GainX,Xs ] = Smooth_state( param,Xf,Xa,PXf,PXa,CL1,CL2,H,Yp,R );

            for e = 1:param.Ne
                Ff( :,e )= Fs( :,e ) ;
                Xf( :,e )= modelrun( param,Xs( :,e )',Fs( :,e ),param.do )' + sqrt(param.me)*randn(param.Nx,1) ;
            end  
            Ffm= mean( Ff,2 );
            Xfm= mean( Xf,2 );
            
            for e = 1:param.Ne
                Xf( :,e )= param.ga(1) * ( Xf( :,e ) - Xfm ) + Xfm ;
                Ff( :,e )= param.ga(2) * ( Ff( :,e ) - Ffm ) + Ffm ; 
            end
            
            Fa= Ff; Fam= mean( Fa,2 );
            Xa= Xf; Xam= mean( Xa,2 ) ;


        elseif ( strcmp(param.aa,'dual-osa') == 1 )

            % Most recent analysis ensemble:
            Xap= Xa - repmat( Xam,1,param.Ne ) ;
            Fap= Fa - repmat( Fam,1,param.Ne ) ;
            
            PFa= 1/( param.Ne-1 ) * ( Fap * Xfp' ) ;
            PXa= 1/( param.Ne-1 ) * ( Xap * Xfp' ) ;

            Fs = Update_param( param,Xf,Fa,PXf,PFa,CL2,H,Yp,R );
            [ GainX,Xs ] = Smooth_state( param,Xf,Xa,PXf,PXa,CL1,CL2,H,Yp,R );

            for e = 1:param.Ne
                Ff( :,e )= Fs( :,e ) ;
                Xf( :,e )= modelrun( param,Xs( :,e )',Fs( :,e ),param.do )' + sqrt(param.me)*randn(param.Nx,1) ;
            end    
            Ffm= mean( Ff,2 );
            Xfm= mean( Xf,2 );
            
            for e = 1:param.Ne
                Xf( :,e )= param.ga(1) * ( Xf( :,e ) - Xfm ) + Xfm ;
                Ff( :,e )= param.ga(2) * ( Ff( :,e ) - Ffm ) + Ffm ; 
            end
            
            Xfp= Xf - repmat( Xfm,1,param.Ne ) ;
            Ffp= Ff - repmat( Ffm,1,param.Ne ) ;

            PXf= 1/( param.Ne-1 ) * ( Xfp * Xfp' ) ;
            PFf= 1/( param.Ne-1 ) * ( Ffp * Xfp' ) ; 
            
            Fa = Update_param( param,Xf,Ff,PXf,PFf,CL2,H,Yp,R );       
            [ ~,Xa ] = Update_state( param,Xf,PXf,CL1,CL2,H,Yp,R ); 
            
            Fam= mean( Fa,2 );
            Xam= mean( Xa,2 );
            
        else
            error('Assimilation algorithm is wrong!')

        end


        Output.mseXf(s)= 1/( param.Nx*param.Ne ) * sum( sum( abs( XR( :,t )*ones(1,param.Ne)-Xf ) ) ) ;
        Output.mseFf(s)= 1/( param.Np*param.Ne ) * sum( sum( abs( param.Fr(:)*ones(1,param.Ne) -Ff ) ) ) ;

        Output.espXf(s)= 1/( param.Nx*param.Ne ) * sum( sum( abs( Xf-Xfm*ones(1,param.Ne) ) ) ) ;
        Output.espFf(s)= 1/( param.Np*param.Ne ) * sum( sum( abs( Ff-Ffm*ones(1,param.Ne) ) ) ) ;
        
        Output.skwXf(s)= 1/( param.Nx*param.Ne ) * sum( sum( ( ( Xf-Xfm*ones(1,param.Ne) )/ ...
                         sqrt( 1/( param.Nx*param.Ne ) * sum( sum( ( Xf-Xfm*ones(1,param.Ne) ).^2 ) ) ) ).^3 ) );
        Output.skwFf(s)= 1/( param.Np*param.Ne ) * sum( sum( ( ( Ff-Ffm*ones(1,param.Ne) )/ ...
                         sqrt( 1/( param.Np*param.Ne ) * sum( sum( ( Ff-Ffm*ones(1,param.Ne) ).^2 ) ) ) ).^3 ) );
        
        Output.kurXf(s)= 1/( param.Nx*param.Ne ) * sum( sum( ( ( Xf-Xfm*ones(1,param.Ne) )/ ...
                         sqrt( 1/( param.Nx*param.Ne ) * sum( sum( ( Xf-Xfm*ones(1,param.Ne) ).^2 ) ) ) ).^4 ) )-3;
        Output.kurFf(s)= 1/( param.Np*param.Ne ) * sum( sum( ( ( Ff-Ffm*ones(1,param.Ne) )/ ...
                         sqrt( 1/( param.Np*param.Ne ) * sum( sum( ( Ff-Ffm*ones(1,param.Ne) ).^2 ) ) ) ).^4 ) )-3;

        Output.dfsXf(s)= trace( GainX'*H' )/param.No;
        Output.corFf(s)= sum( PFf(param.Ip,:)*H' )/param.Fr( param.Ip );
        
        Output.eskXf(s)= 1/( param.Nx ) * sum( ( XR( :,t )-Xfm ).^2 ) ;
        Output.eskFf(s)= 1/( param.Np ) * sum( ( param.Fr(:)-Ffm ).^2 ) ;

        Output.inbXf(s)= mean(Output.eskXf)/mean(Output.espXf);
        Output.indXf(s)= mean(sqrt(Output.eskXf))/mean(sqrt(Output.mseXf));

        Output.dXf(s,:)= Xf(param.Ix,:); Output.refX(s)= XR( param.Ix,t );
        Output.dFf(s,:)= Ff(param.Ip,:); Output.refF(s)= param.Fr( param.Ip );

    end
    
end




%         elseif ( strcmp(algo,'w-EnKF') == 1 )
% 
%             % Apply Kalman-Correction to each state&parameter member:
%             GainX = ( PXf*H' ) / ( H*PXf*H' + R ) ;
%             GainF = ( PFf*H' ) / ( H*PXf*H' + R ) ;
% %             Vaprx = dia( diag(inv(eye(param.No)-H*GainX)) );
%             for e = 1:param.Ne
%                 Xa( :,e )= Xf( :,e ) + GainX * ( Yp( :,e ) - H*Xf( :,e ) ) ;
%             end     
%             for e = 1:param.Ne
%                 Fa( :,e )= Ff( :,e ) + GainF * ( Yp( :,e ) - ...
%                            ( eye(param.No)-H*GainX )\( H*Xa(:,e)-H*GainX*Yp(:,e) ) ) ;
%             end