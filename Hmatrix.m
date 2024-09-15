function [ H,ObsLoc ] = Hmatrix( Nx,cas )

    switch cas 

      case 'all'
        hdiag = zeros(Nx,1) ;
        hdiag(1:1:end) = 1  ;
        izeros =  ~hdiag  ;
        H = diag(hdiag) ;
        H( izeros,: ) = [] ; 

      case 'half'
        hdiag = zeros(Nx,1) ;
        hdiag(1:2:end) = 1  ;
        izeros =  ~hdiag  ;
        H = diag(hdiag) ;
        H( izeros,: ) = [] ;

      case 'quarter'
        hdiag = zeros(Nx,1) ;
        hdiag(1:3:end) = 1  ;
        izeros =  ~hdiag  ;
        H = diag(hdiag) ;
        H( izeros,: ) = [] ;

    end

    ObsLoc = find( hdiag>0 ) ; 
    
end