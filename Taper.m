% function coeffs = Taper( dist,lc ) 
% 
%     R = lc * sqrt(10/3) ;
% 
%     coeffs = zeros( size(dist) ) ;            
%     ind1 = find(dist <= R ) ;
%     r2 = ( dist(ind1) / R ) .^ 2 ;
%     r3 = ( dist(ind1) / R ) .^ 3 ;
%     coeffs( ind1 ) = 1 + r2 .* (- r3 / 4 + r2 / 2) + r3 * (5 / 8) - r2 * (5 / 3) ;
%     ind2 = find( dist > R & dist <= R * 2 ) ;
%     r1 = ( dist(ind2) / R ) ;
%     r2 = ( dist(ind2) / R ) .^ 2 ;
%     r3 = ( dist(ind2) / R ) .^ 3 ;
%     coeffs( ind2 ) = r2 .* (r3 / 12 - r2 / 2) + r3 * (5 / 8) + r2 * (5 / 3) - r1 * 5 + 4 - (2 / 3) ./ r1 ;
% 
% end


function Taper = Taper( Nx,P,lc,type ) 

    Ny = 1 ;
    p  = length( P ) ;
    [ I,J ] = ind2sub( [ Nx,Ny ],P ) ;
    pos( :,1 ) = I ; 
    pos( :,2 ) = J ;

    %%% Calculating the distances to the observation points
    R = lc * sqrt(10/3) ;

    switch type
        case 'PH'
            Taper = zeros( Nx*Ny,p ) ; 
            k = 0 ;

            for ii = 1:Nx
                for jj = 1:Ny
                    k = k + 1 ;
                    dist = hypot( ii - pos( :,1 ), jj - pos( :,2 ) ) ;

                    coeffs = zeros( size(dist) ) ;            
                    ind1 = find( dist <= R ) ;
                    r2 = ( dist(ind1) / R ) .^ 2 ;
                    r3 = ( dist(ind1) / R ) .^ 3 ;
                    coeffs( ind1 ) = 1 + r2 .* (- r3 / 4 + r2 / 2) + r3 * (5 / 8) - r2 * (5 / 3) ;
                    ind2 = find( dist > R & dist <= R * 2 ) ; 
                    r1 = ( dist(ind2) / R ) ;
                    r2 = ( dist(ind2) / R ) .^ 2 ;
                    r3 = ( dist(ind2) / R ) .^ 3 ;
                    coeffs( ind2 ) = r2 .* (r3 / 12 - r2 / 2) + r3 * (5 / 8) + r2 * (5 / 3) - r1 * 5 + 4 - (2 / 3) ./ r1 ;

                    tmp = coeffs < 1e-3 ;
                    coeffs( tmp )  = 0 ;
                    Taper( k,: ) = coeffs' ;
                end
            end

        case 'HPH'
            Taper = zeros( p,p ) ; 
            k = 0 ;
            
            for i = 1:p
                k = k + 1 ;
                dist = hypot( pos( i,1 ) - pos( :,1 ), pos( i,2 ) - pos( :,2 ) ) ;

                coeffs = zeros( size(dist) );            
                ind1 = find(dist <= R ) ;
                r2 = ( dist(ind1) / R ) .^ 2;
                r3 = ( dist(ind1) / R ) .^ 3;
                coeffs( ind1 ) = 1 + r2 .* (- r3 / 4 + r2 / 2) + r3 * (5 / 8) - r2 * (5 / 3);
                ind2 = find( dist > R & dist <= R * 2 ) ;
                r1 = ( dist(ind2) / R );
                r2 = ( dist(ind2) / R ) .^ 2;
                r3 = ( dist(ind2) / R ) .^ 3;
                coeffs( ind2 ) = r2 .* (r3 / 12 - r2 / 2) + r3 * (5 / 8) + r2 * (5 / 3) - r1 * 5 + 4 - (2 / 3) ./ r1 ;

                Taper( k,: ) = coeffs' ;
            end
    end
    
end