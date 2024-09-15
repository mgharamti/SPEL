
function fX = florenz63(Nx,X,F)

    % Set fX(1).
    fX(1) = F(1)*( X(2) - X(1) ) ;

    % Set fX(2).
    fX(2) = X(1)*( F(2) - X(3) ) - X(2) ;

    % Set fX(3).
    fX(Nx) = X(1)*X(2) - F(3)*X(3) ;

end
