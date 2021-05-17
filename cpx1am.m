function C=cpx1am(U,C0,U0,m)
if (U<=0.9*U0)
    C=C0/(1-U/U0)^m;
else
    C=C0*10^m*(1+10*m*U/U0-9*m);
end