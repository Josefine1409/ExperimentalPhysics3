clear all; close all;clc;

syms phi theta p_i p_f p_t m1 m2 K 

assumeAlso(K>0)
assumeAlso(p_i>0)
assumeAlso(p_t>0)
assumeAlso(m1>0)
assumeAlso(m2>0)
assumeAlso(m2>m1)

eq3 = K*p_i*sin(theta)==p_t*sin(phi);
eq2 = p_i== K*p_i*cos(theta)+p_t*cos(phi);
eq1 = p_i^2/m1 == K^2*p_i^2/m1+p_t^2/m2;



eqns=[eq1,eq2,eq3]
S = solve(eqns,[K,p_t,phi])
