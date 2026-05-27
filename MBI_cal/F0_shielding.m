function f = F0_shielding(b)
% F0_shielding  Auxiliary function for parallel-plate CSR shielding
%
% Computes the F0(beta) function using Airy functions.
% Ref: Agoh & Yokoya, PRSTAB 7, 054403 (2004), Appendix

    if b > 10
        % Asymptotic expansion for large beta
        t1 = b/(2*pi) * exp(-4/3*b^3);
        t2 = 1 + 1/(24*b^3) + 1/(1152*b^6);
        t3 = -3/(16*pi*b^5);
        t4 = 1 + 105/32/b^6;
        f  = t1*t2 + 1i*t3*t4;
    else
        b2 = b^2;
        t1 = airy(1,b2).*airy(1,b2) - 1i*airy(1,b2).*airy(3,b2);
        t2 = b2.*airy(0,b2).*airy(0,b2) - 1i*b2.*airy(0,b2).*airy(2,b2);
        f  = t1 + t2;
    end
end
