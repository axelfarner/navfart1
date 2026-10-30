function n = open_loop_speed_control(U_ref)
    Ja   = 0;       % advance number (bollard pull)
    PD   = 1.5;     % pitch/diameter ratio
    AEAO = 0.65;    % blade area ratio
    z    = 4;       % number of propeller blades
    KT = wageningen(Ja,PD,AEAO,z);
    D = 3.3;
    rho = 1025;
    m = 17.0677e6;          % mass (kg)
    T1 = 20;
    Xudot = -8.9830e5;
    Xu = -(m-Xudot)/T1;
    t = 0.05;
    T_d = U_ref * Xu / (t - 1);
    
    n = sqrt(T_d / (rho * KT * D ^4));

end