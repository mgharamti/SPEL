
function HX = obs_operator( X,caso )

    % All variables are observed.
    if strcmp(caso,'all') == 1
      HX = X;
    end

    % Every other variable is observed.
    if strcmp(caso,'half') == 1
      HX = X(1:2:end,:);
    end

    % Every other 3 variables is observed.
    if strcmp(caso,'quarter') == 1
      HX = X(1:3:end,:);
    end

    % Nonlinear observation operator - 1.
    if strcmp(caso,'nl1') == 1
      HX = X(1:2:end,:).*X(2:2:end,:);
    end

    % Nonlinear observation operator - 2.
    if strcmp(caso,'nl2') == 1
      HX = X(1:2:end,:).*X(1:2:end,:);
    end

end