function [S] = stoppingpowerC(E)
    rho = 2.267;
    A = load('CStoppingPower.txt');
    Es = A(1:2:end,1);
    Ss = A(1:2:end,2)*rho*100;
    S = spline(Es,Ss,E);
end