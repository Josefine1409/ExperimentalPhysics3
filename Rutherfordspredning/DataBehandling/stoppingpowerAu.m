function [S] = stoppingpowerAu(E)
    rho =19.30;
    A = load('AuStoppingPower.txt');
    Es = A(1:2:end,1);
    Ss = A(1:2:end,2)*rho*100;
    S = spline(Es,Ss,E);
end