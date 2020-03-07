function w = waldweight(V);
%
%   This Matlab program computes the weights for the one-sided Wald test
%   under the null hypothesis.
%
% From Raymond Kan, 16 Sep 2010
n = length(V);
w = zeros(n+1,1);
index = false(n,1);
for i=1:2^n-2
    j = 1;
    while index(j)
        j = j+1;
    end
    index(1:j-1) = false;
    index(j) = true;
    V22i = inv(V(~index,~index));
    V12 = V(index,~index);
    V11 = V(index,index)-V12*V22i*V12';
    k = sum(index)+1;
    w(k) = w(k)+posprob(V22i)*posprob(V11);
end
if rem(n,2)==0
    w(end) = posprob(V);
else
    w(end) = 0.5-sum(w(2:2:n-1));
end
w(1) = 1-sum(w(2:end));
end

function y = posprob(A)
%
%   This Matlab function computes the probability for X>0 when
%   X ~ N(0,A).  The dimension of A has to be less than or equal
%   to 11 for us to use the fast algorithm, which is based on
%   Sun, Hong-Jie(1988) 'A general reduction method for n-variate
%   normal orthant probability', Communications in Statistics -
%   Theory and Methods, 17: 11, 3913 - 3921
%   Sun, Hong-Jie(1998) 'A Fortran subroutine for computing normal
%   orthant probabilities of dimensions up to nine', Communications
%   in Statistics - Simulation and Computation, 17: 3, 1097 - 1111
%
p = length(A);
if p>3&&p<12
    % 12 point Gauss-Legendre quadrature, nodes and weights are
    % obtained using Mathematica
    nw = 12;
    xg = [9.90780317123359625e-1
        9.52058628185237428e-1
        8.84951337097152343e-1
        7.93658977143308724e-1
        6.83915749499090097e-1
        5.62616704255734458e-1
        4.37383295744265542e-1
        3.16084250500909903e-1
        2.06341022856691276e-1
        1.15048662902847656e-1
        4.79413718147625717e-2
        9.21968287664037465e-3];
    wg = [2.35876681932559136e-2
        5.34696629976592155e-2
        8.00391642716731132e-2
        1.01583713361532961e-1
        1.16746268269177404e-1
        1.24573522906701393e-1
        1.24573522906701393e-1
        1.16746268269177404e-1
        1.01583713361532961e-1
        8.00391642716731132e-2
        5.34696629976592155e-2
        2.35876681932559136e-2];
    %     [xg,wg] = lgwt(nw,0,1);
    x2 = xg.*xg;
end
%
%   Normalize the matrix
%
d = sqrt(diag(A));
A = A./(d*d');
%  Use analytical or integration formula for p<=9
if p==1
    y = 1/2;
elseif p==2
    y = 0.25+asin(A(2))/(2*pi);
elseif p==3
    y = 0.125+sum(asin(A([2 3 6])))/(4*pi);
elseif p==4
    r = [A(2:4,1); A(3:4,2); A(4,3)]';
    f = fi4(r);
    y = 0.0625+sum(asin(r))/(8*pi)+f/(4*pi^2);
elseif p==5
    r = [A(2:5,1); A(3:5,2); A(4:5,3); A(5,4)]';
    ind = setindex(p,4);
    f = sum(fi4(r(ind)));
    y = 1/32+sum(asin(r))/(16*pi)+f/(8*pi^2);
elseif p==6
    r = [A(2:6,1); A(3:6,2); A(4:6,3); A(5:6,4); A(6,5)]';
    ind = setindex(p,4);
    f1 = sum(fi4(r(ind)));
    f2 = fi6(r);
    y = 1/64+sum(asin(r))/(32*pi)+f1/(16*pi^2)+f2/(8*pi^3);
elseif p==7
    r = [A(2:7,1); A(3:7,2); A(4:7,3); A(5:7,4); A(6:7,5); A(7,6)]';
    ind = setindex(p,4);
    f1 = sum(fi4(r(ind)));
    ind = setindex(p,6);
    f2 = sum(fi6(r(ind)));
    y = 1/128+sum(asin(r))/(64*pi)+f1/(32*pi^2)+f2/(16*pi^3);
elseif p==8
    r = [A(2:8,1); A(3:8,2); A(4:8,3); A(5:8,4); A(6:8,5); A(7:8,6); A(8,7)]';
    ind = setindex(p,4);
    f1 = sum(fi4(r(ind)));
    ind = setindex(p,6);
    f2 = sum(fi6(r(ind)));
    f3 = fi8(r);
    y = 1/256+sum(asin(r))/(128*pi)+f1/(64*pi^2)+f2/(32*pi^3)+f3/(16*pi^4);
elseif p==9
    r = [A(2:9,1); A(3:9,2); A(4:9,3); A(5:9,4); A(6:9,5); A(7:9,6); A(8:9,7); A(9,8)]';
    ind = setindex(p,4);
    f1 = sum(fi4(r(ind)));
    ind = setindex(p,6);
    f2 = sum(fi6(r(ind)));
    ind = setindex(p,8);
    f3 = sum(fi8(r(ind)));
    y = 1/512+sum(asin(r))/(256*pi)+f1/(128*pi^2)+f2/(64*pi^3)+f3/(32*pi^4);
elseif p==10
    r = [A(2:10,1); A(3:10,2); A(4:10,3); A(5:10,4); A(6:10,5); A(7:10,6); A(8:10,7); A(9:10,8); A(10,9)]';
    ind = setindex(p,4);
    f1 = sum(fi4(r(ind)));
    ind = setindex(p,6);
    f2 = sum(fi6(r(ind)));
    ind = setindex(p,8);
    f3 = sum(fi8(r(ind)));
    f4 = fi10(r);
    y = 1/1024+sum(asin(r))/(512*pi)+f1/(256*pi^2)+f2/(128*pi^3)+f3/(64*pi^4)+f4/(32*pi^5);
elseif p==11
    r = [A(2:11,1); A(3:11,2); A(4:11,3); A(5:11,4); A(6:11,5); A(7:11,6); A(8:11,7); A(9:11,8); A(10:11,9); A(11,10)]';
    ind = setindex(p,4);
    f1 = sum(fi4(r(ind)));
    ind = setindex(p,6);
    f2 = sum(fi6(r(ind)));
    ind = setindex(p,8);
    f3 = sum(fi8(r(ind)));
    ind = setindex(p,10);
    f4 = sum(fi10(r(ind)));
    y = 1/2048+sum(asin(r))/(1024*pi)+f1/(512*pi^2)+f2/(256*pi^3)+f3/(128*pi^4)+f4/(64*pi^5);
else
    options = optimset('TolFun',1e-4);
    y = mvncdf(zeros(p,1),zeros(p,1),A,options);
end

    function z = fi4(rr)
        %
        %  This function computes the integral of I_4(S1,S2,S3,S4,S5,S6),
        %  where S1 to S6 are the correlation coefficients.  When Si's are
        %  vectors, the function returns a vector output of the integral.
        %
        z = zeros(size(rr,1),1);
        for i=1:size(rr,1)
            S1 = rr(i,1);  S2 = rr(i,2);  S3 = rr(i,3);  S4 = rr(i,4);  S5 = rr(i,5);  S6 = rr(i,6);
            c1 = 1-S4^2-(S1^2+S2^2-2*S1*S2*S4)*x2;
            c2 = 1-S5^2-(S1^2+S3^2-2*S1*S3*S5)*x2;
            c3 = 1-S6^2-(S2^2+S3^2-2*S2*S3*S6)*x2;
            R1 = (S6-S4*S5-(S2*S3+S1^2*S6-S1*S3*S4-S1*S2*S5)*x2)./sqrt(c1.*c2);
            R2 = (S5-S4*S6-(S1*S3+S2^2*S5-S2*S1*S6-S2*S3*S4)*x2)./sqrt(c1.*c3);
            R3 = (S4-S5*S6-(S1*S2+S3^2*S4-S1*S3*S6-S2*S3*S5)*x2)./sqrt(c2.*c3);
            z(i) = wg'*(S1./sqrt(1-S1^2*x2).*asin(R1)+S2./sqrt(1-S2^2*x2).*asin(R2)+ ...
                S3./sqrt(1-S3^2*x2).*asin(R3));
        end
    end

    function z = fi6(rr)
        %
        %  This function computes the integral for I_6(S1,...,S15),
        %  where S1 to S15 are the correlation coefficients.  When Si's are
        %  vectors, the function returns a vector output of the integral.
        %
        z = zeros(size(rr,1),1);
        for i=1:size(rr,1)
            S1 = rr(i,1);  S2 = rr(i,2);  S3 = rr(i,3);  S4 = rr(i,4);  S5 = rr(i,5);
            S6 = rr(i,6);  S7 = rr(i,7);  S8 = rr(i,8);  S9 = rr(i,9);  S10 = rr(i,10);
            S11 = rr(i,11);  S12 = rr(i,12);  S13 = rr(i,13);  S14 = rr(i,14);  S15 = rr(i,15);
            c1 = 1-S6^2-(S1^2+S2^2-2*S1*S2*S6)*x2;
            c2 = 1-S7^2-(S1^2+S3^2-2*S1*S3*S7)*x2;
            c3 = 1-S8^2-(S1^2+S4^2-2*S1*S4*S8)*x2;
            c4 = 1-S9^2-(S1^2+S5^2-2*S1*S5*S9)*x2;
            c5 = 1-S10^2-(S2^2+S3^2-2*S2*S3*S10)*x2;
            c6 = 1-S11^2-(S2^2+S4^2-2*S2*S4*S11)*x2;
            c7 = 1-S12^2-(S2^2+S5^2-2*S2*S5*S12)*x2;
            c8 = 1-S13^2-(S3^2+S4^2-2*S3*S4*S13)*x2;
            c9 = 1-S14^2-(S3^2+S5^2-2*S3*S5*S14)*x2;
            c10 = 1-S15^2-(S4^2+S5^2-2*S4*S5*S15)*x2;
            %
            %  i=2
            %
            R1 = (S10-S6*S7-(S2*S3+S1^2*S10-S1*S2*S7-S1*S3*S6)*x2)./sqrt(c1.*c2);
            R2 = (S11-S6*S8-(S2*S4+S1^2*S11-S1*S2*S8-S1*S4*S6)*x2)./sqrt(c1.*c3);
            R3 = (S12-S6*S9-(S2*S5+S1^2*S12-S1*S2*S9-S1*S5*S6)*x2)./sqrt(c1.*c4);
            R4 = (S13-S7*S8-(S3*S4+S1^2*S13-S1*S3*S8-S1*S4*S7)*x2)./sqrt(c2.*c3);
            R5 = (S14-S7*S9-(S3*S5+S1^2*S14-S1*S3*S9-S1*S5*S7)*x2)./sqrt(c2.*c4);
            R6 = (S15-S8*S9-(S4*S5+S1^2*S15-S1*S4*S9-S1*S5*S8)*x2)./sqrt(c3.*c4);
            z1 = S1*fi4([R1 R2 R3 R4 R5 R6])./sqrt(1-S1^2*x2);
            %
            %  i=3
            %
            R1 = (S7-S6*S10-(S1*S3+S2^2*S7-S2*S1*S10-S2*S3*S6)*x2)./sqrt(c1.*c5);
            R2 = (S8-S6*S11-(S1*S4+S2^2*S8-S2*S1*S11-S2*S4*S6)*x2)./sqrt(c1.*c6);
            R3 = (S9-S6*S12-(S1*S5+S2^2*S9-S2*S1*S12-S2*S5*S6)*x2)./sqrt(c1.*c7);
            R4 = (S13-S10*S11-(S3*S4+S2^2*S13-S2*S3*S11-S2*S4*S10)*x2)./sqrt(c5.*c6);
            R5 = (S14-S10*S12-(S3*S5+S2^2*S14-S2*S5*S10-S2*S3*S12)*x2)./sqrt(c5.*c7);
            R6 = (S15-S11*S12-(S4*S5+S2^2*S15-S2*S4*S12-S2*S5*S11)*x2)./sqrt(c6.*c7);
            z2 = S2*fi4([R1 R2 R3 R4 R5 R6])./sqrt(1-S2^2*x2);
            %
            %  i=4
            %
            R1 = (S6-S7*S10-(S1*S2+S3^2*S6-S3*S1*S10-S3*S2*S7)*x2)./sqrt(c2.*c5);
            R2 = (S8-S7*S13-(S1*S4+S3^2*S8-S3*S1*S13-S3*S4*S7)*x2)./sqrt(c2.*c8);
            R3 = (S9-S7*S14-(S1*S5+S3^2*S9-S3*S1*S14-S3*S5*S7)*x2)./sqrt(c2.*c9);
            R4 = (S11-S10*S13-(S2*S4+S3^2*S11-S3*S2*S13-S3*S4*S10)*x2)./sqrt(c5.*c8);
            R5 = (S12-S10*S14-(S2*S5+S3^2*S12-S3*S2*S14-S3*S5*S10)*x2)./sqrt(c5.*c9);
            R6 = (S15-S13*S14-(S4*S5+S3^2*S15-S3*S4*S14-S3*S5*S13)*x2)./sqrt(c8.*c9);
            z3 = S3*fi4([R1 R2 R3 R4 R5 R6])./sqrt(1-S3^2*x2);
            %
            %  i=5
            %
            R1 = (S6-S8*S11-(S1*S2+S4^2*S6-S4*S1*S11-S4*S2*S8)*x2)./sqrt(c3.*c6);
            R2 = (S7-S8*S13-(S1*S3+S4^2*S7-S4*S1*S13-S4*S3*S8)*x2)./sqrt(c3.*c8);
            R3 = (S9-S8*S15-(S1*S5+S4^2*S9-S4*S1*S15-S4*S5*S8)*x2)./sqrt(c3.*c10);
            R4 = (S10-S11*S13-(S2*S3+S4^2*S10-S4*S2*S13-S4*S3*S11)*x2)./sqrt(c6.*c8);
            R5 = (S12-S11*S15-(S2*S5+S4^2*S12-S4*S2*S15-S4*S5*S11)*x2)./sqrt(c6.*c10);
            R6 = (S14-S13*S15-(S3*S5+S4^2*S14-S4*S3*S15-S4*S5*S13)*x2)./sqrt(c8.*c10);
            z4 = S4*fi4([R1 R2 R3 R4 R5 R6])./sqrt(1-S4^2*x2);
            %
            %  i=6
            %
            R1 = (S6-S9*S12-(S1*S2+S5^2*S6-S5*S1*S12-S5*S2*S9)*x2)./sqrt(c4.*c7);
            R2 = (S7-S9*S14-(S1*S3+S5^2*S7-S5*S1*S14-S5*S3*S9)*x2)./sqrt(c4.*c9);
            R3 = (S8-S9*S15-(S1*S4+S5^2*S8-S5*S1*S15-S5*S4*S9)*x2)./sqrt(c4.*c10);
            R4 = (S10-S12*S14-(S2*S3+S5^2*S10-S5*S2*S14-S5*S3*S12)*x2)./sqrt(c7.*c9);
            R5 = (S11-S12*S15-(S2*S4+S5^2*S11-S5*S2*S15-S5*S4*S12)*x2)./sqrt(c7.*c10);
            R6 = (S13-S14*S15-(S3*S4+S5^2*S13-S5*S3*S15-S5*S4*S14)*x2)./sqrt(c9.*c10);
            z5 = S5*fi4([R1 R2 R3 R4 R5 R6])./sqrt(1-S5^2*x2);
            z(i) = wg'*(z1+z2+z3+z4+z5);
        end
    end

    function z = fi8(rr)
        %
        %  This function computes the integral for I_8(S1,...,S28),
        %  where S1 to S28 are the correlation coefficients.  When Si's are
        %  vectors, the function returns a vector output of the integral.
        %
        z = zeros(size(rr,1),1);
        for i=1:size(rr,1)
            S1 = rr(i,1);  S2 = rr(i,2);  S3 = rr(i,3);  S4 = rr(i,4);  S5 = rr(i,5);
            S6 = rr(i,6);  S7 = rr(i,7);  S8 = rr(i,8);  S9 = rr(i,9);  S10 = rr(i,10);
            S11 = rr(i,11);  S12 = rr(i,12);  S13 = rr(i,13);  S14 = rr(i,14);  S15 = rr(i,15);
            S16 = rr(i,16);  S17 = rr(i,17);  S18 = rr(i,18);  S19 = rr(i,19);  S20 = rr(i,20);
            S21 = rr(i,21);  S22 = rr(i,22);  S23 = rr(i,23);  S24 = rr(i,24);  S25 = rr(i,25);
            S26 = rr(i,26);  S27 = rr(i,27);  S28 = rr(i,28);
            c1 = 1-S8^2-(S1^2+S2^2-2*S1*S2*S8)*x2;
            c2 = 1-S9^2-(S1^2+S3^2-2*S1*S3*S9)*x2;
            c3 = 1-S10^2-(S1^2+S4^2-2*S1*S4*S10)*x2;
            c4 = 1-S11^2-(S1^2+S5^2-2*S1*S5*S11)*x2;
            c5 = 1-S12^2-(S1^2+S6^2-2*S1*S6*S12)*x2;
            c6 = 1-S13^2-(S1^2+S7^2-2*S1*S7*S13)*x2;
            c7 = 1-S14^2-(S2^2+S3^2-2*S2*S3*S14)*x2;
            c8 = 1-S15^2-(S2^2+S4^2-2*S2*S4*S15)*x2;
            c9 = 1-S16^2-(S2^2+S5^2-2*S2*S5*S16)*x2;
            c10 = 1-S17^2-(S2^2+S6^2-2*S2*S6*S17)*x2;
            c11 = 1-S18^2-(S2^2+S7^2-2*S2*S7*S18)*x2;
            c12 = 1-S19^2-(S3^2+S4^2-2*S3*S4*S19)*x2;
            c13 = 1-S20^2-(S3^2+S5^2-2*S3*S5*S20)*x2;
            c14 = 1-S21^2-(S3^2+S6^2-2*S3*S6*S21)*x2;
            c15 = 1-S22^2-(S3^2+S7^2-2*S3*S7*S22)*x2;
            c16 = 1-S23^2-(S4^2+S5^2-2*S4*S5*S23)*x2;
            c17 = 1-S24^2-(S4^2+S6^2-2*S4*S6*S24)*x2;
            c18 = 1-S25^2-(S4^2+S7^2-2*S4*S7*S25)*x2;
            c19 = 1-S26^2-(S5^2+S6^2-2*S5*S6*S26)*x2;
            c20 = 1-S27^2-(S5^2+S7^2-2*S5*S7*S27)*x2;
            c21 = 1-S28^2-(S6^2+S7^2-2*S6*S7*S28)*x2;
            %
            %  i=2
            %
            R1 = ((S1*S2*S9-S1^2*S14-S2*S3+S1*S3*S8)*x2+S14-S8*S9)./sqrt(c1.*c2);
            R2 = ((S1*S10*S2-S1^2*S15-S2*S4+S1*S4*S8)*x2+S15-S10*S8)./sqrt(c1.*c3);
            R3 = ((S1*S11*S2-S1^2*S16-S2*S5+S1*S5*S8)*x2+S16-S11*S8)./sqrt(c1.*c4);
            R4 = ((S1*S12*S2-S1^2*S17-S2*S6+S1*S6*S8)*x2+S17-S12*S8)./sqrt(c1.*c5);
            R5 = ((S1*S13*S2-S1^2*S18-S2*S7+S1*S7*S8)*x2+S18-S13*S8)./sqrt(c1.*c6);
            R6 = ((S1*S10*S3-S1^2*S19-S3*S4+S1*S4*S9)*x2+S19-S10*S9)./sqrt(c2.*c3);
            R7 = ((S1*S11*S3-S1^2*S20-S3*S5+S1*S5*S9)*x2+S20-S11*S9)./sqrt(c2.*c4);
            R8 = ((S1*S12*S3-S1^2*S21-S3*S6+S1*S6*S9)*x2+S21-S12*S9)./sqrt(c2.*c5);
            R9 = ((S1*S13*S3-S1^2*S22-S3*S7+S1*S7*S9)*x2+S22-S13*S9)./sqrt(c2.*c6);
            R10 = ((S1*S11*S4-S1^2*S23-S4*S5+S1*S10*S5)*x2+S23-S10*S11)./sqrt(c3.*c4);
            R11 = ((S1*S12*S4-S1^2*S24-S4*S6+S1*S10*S6)*x2+S24-S10*S12)./sqrt(c3.*c5);
            R12 = ((S1*S13*S4-S1^2*S25-S4*S7+S1*S10*S7)*x2+S25-S10*S13)./sqrt(c3.*c6);
            R13 = ((S1*S12*S5-S1^2*S26-S5*S6+S1*S11*S6)*x2+S26-S11*S12)./sqrt(c4.*c5);
            R14 = ((S1*S13*S5-S1^2*S27-S5*S7+S1*S11*S7)*x2+S27-S11*S13)./sqrt(c4.*c6);
            R15 = ((S1*S13*S6-S1^2*S28-S6*S7+S1*S12*S7)*x2+S28-S12*S13)./sqrt(c5.*c6);
            z1 = S1*fi6([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15])./sqrt(1-S1^2*x2);
            %
            %  i=3
            %
            R1 = ((S1*S14*S2-S2^2*S9-S1*S3+S2*S3*S8)*x2+S9-S14*S8)./sqrt(c1.*c7);
            R2 = ((S1*S15*S2-S10*S2^2-S1*S4+S2*S4*S8)*x2+S10-S15*S8)./sqrt(c1.*c8);
            R3 = ((S1*S16*S2-S11*S2^2-S1*S5+S2*S5*S8)*x2+S11-S16*S8)./sqrt(c1.*c9);
            R4 = ((S1*S17*S2-S12*S2^2-S1*S6+S2*S6*S8)*x2+S12-S17*S8)./sqrt(c1.*c10);
            R5 = ((S1*S18*S2-S13*S2^2-S1*S7+S2*S7*S8)*x2+S13-S18*S8)./sqrt(c1.*c11);
            R6 = ((S15*S2*S3-S19*S2^2-S3*S4+S14*S2*S4)*x2+S19-S14*S15)./sqrt(c7.*c8);
            R7 = ((S16*S2*S3-S2^2*S20-S3*S5+S14*S2*S5)*x2+S20-S14*S16)./sqrt(c7.*c9);
            R8 = ((S17*S2*S3-S2^2*S21-S3*S6+S14*S2*S6)*x2+S21-S14*S17)./sqrt(c7.*c10);
            R9 = ((S18*S2*S3-S2^2*S22-S3*S7+S14*S2*S7)*x2+S22-S14*S18)./sqrt(c7.*c11);
            R10 = ((S16*S2*S4-S2^2*S23-S4*S5+S15*S2*S5)*x2+S23-S15*S16)./sqrt(c8.*c9);
            R11 = ((S17*S2*S4-S2^2*S24-S4*S6+S15*S2*S6)*x2+S24-S15*S17)./sqrt(c8.*c10);
            R12 = ((S18*S2*S4-S2^2*S25-S4*S7+S15*S2*S7)*x2+S25-S15*S18)./sqrt(c8.*c11);
            R13 = ((S17*S2*S5-S2^2*S26-S5*S6+S16*S2*S6)*x2+S26-S16*S17)./sqrt(c9.*c10);
            R14 = ((S18*S2*S5-S2^2*S27-S5*S7+S16*S2*S7)*x2+S27-S16*S18)./sqrt(c9.*c11);
            R15 = ((S18*S2*S6-S2^2*S28-S6*S7+S17*S2*S7)*x2+S28-S17*S18)./sqrt(c10.*c11);
            z2 = S2*fi6([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15])./sqrt(1-S2^2*x2);
            %
            %  i=4
            %
            R1 = ((S1*S14*S3-S3^2*S8-S1*S2+S2*S3*S9)*x2+S8-S14*S9)./sqrt(c2.*c7);
            R2 = ((S1*S19*S3-S10*S3^2-S1*S4+S3*S4*S9)*x2+S10-S19*S9)./sqrt(c2.*c12);
            R3 = ((S1*S20*S3-S11*S3^2-S1*S5+S3*S5*S9)*x2+S11-S20*S9)./sqrt(c2.*c13);
            R4 = ((S1*S21*S3-S12*S3^2-S1*S6+S3*S6*S9)*x2+S12-S21*S9)./sqrt(c2.*c14);
            R5 = ((S1*S22*S3-S13*S3^2-S1*S7+S3*S7*S9)*x2+S13-S22*S9)./sqrt(c2.*c15);
            R6 = ((S19*S2*S3-S15*S3^2-S2*S4+S14*S3*S4)*x2+S15-S14*S19)./sqrt(c7.*c12);
            R7 = ((S2*S20*S3-S16*S3^2-S2*S5+S14*S3*S5)*x2+S16-S14*S20)./sqrt(c7.*c13);
            R8 = ((S2*S21*S3-S17*S3^2-S2*S6+S14*S3*S6)*x2+S17-S14*S21)./sqrt(c7.*c14);
            R9 = ((S2*S22*S3-S18*S3^2-S2*S7+S14*S3*S7)*x2+S18-S14*S22)./sqrt(c7.*c15);
            R10 = ((S20*S3*S4-S23*S3^2-S4*S5+S19*S3*S5)*x2+S23-S19*S20)./sqrt(c12.*c13);
            R11 = ((S21*S3*S4-S24*S3^2-S4*S6+S19*S3*S6)*x2+S24-S19*S21)./sqrt(c12.*c14);
            R12 = ((S22*S3*S4-S25*S3^2-S4*S7+S19*S3*S7)*x2+S25-S19*S22)./sqrt(c12.*c15);
            R13 = ((S21*S3*S5-S26*S3^2-S5*S6+S20*S3*S6)*x2+S26-S20*S21)./sqrt(c13.*c14);
            R14 = ((S22*S3*S5-S27*S3^2-S5*S7+S20*S3*S7)*x2+S27-S20*S22)./sqrt(c13.*c15);
            R15 = ((S22*S3*S6-S28*S3^2-S6*S7+S21*S3*S7)*x2+S28-S21*S22)./sqrt(c14.*c15);
            z3 = S3*fi6([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15])./sqrt(1-S3^2*x2);
            %
            %  i=5
            %
            R1 = ((S1*S15*S4-S4^2*S8-S1*S2+S10*S2*S4)*x2+S8-S10*S15)./sqrt(c3.*c8);
            R2 = ((S1*S19*S4-S4^2*S9-S1*S3+S10*S3*S4)*x2+S9-S10*S19)./sqrt(c3.*c12);
            R3 = ((S1*S23*S4-S11*S4^2-S1*S5+S10*S4*S5)*x2+S11-S10*S23)./sqrt(c3.*c16);
            R4 = ((S1*S24*S4-S12*S4^2-S1*S6+S10*S4*S6)*x2+S12-S10*S24)./sqrt(c3.*c17);
            R5 = ((S1*S25*S4-S13*S4^2-S1*S7+S10*S4*S7)*x2+S13-S10*S25)./sqrt(c3.*c18);
            R6 = ((S19*S2*S4-S14*S4^2-S2*S3+S15*S3*S4)*x2+S14-S15*S19)./sqrt(c8.*c12);
            R7 = ((S2*S23*S4-S16*S4^2-S2*S5+S15*S4*S5)*x2+S16-S15*S23)./sqrt(c8.*c16);
            R8 = ((S2*S24*S4-S17*S4^2-S2*S6+S15*S4*S6)*x2+S17-S15*S24)./sqrt(c8.*c17);
            R9 = ((S2*S25*S4-S18*S4^2-S2*S7+S15*S4*S7)*x2+S18-S15*S25)./sqrt(c8.*c18);
            R10 = ((S23*S3*S4-S20*S4^2-S3*S5+S19*S4*S5)*x2+S20-S19*S23)./sqrt(c12.*c16);
            R11 = ((S24*S3*S4-S21*S4^2-S3*S6+S19*S4*S6)*x2+S21-S19*S24)./sqrt(c12.*c17);
            R12 = ((S25*S3*S4-S22*S4^2-S3*S7+S19*S4*S7)*x2+S22-S19*S25)./sqrt(c12.*c18);
            R13 = ((S24*S4*S5-S26*S4^2-S5*S6+S23*S4*S6)*x2+S26-S23*S24)./sqrt(c16.*c17);
            R14 = ((S25*S4*S5-S27*S4^2-S5*S7+S23*S4*S7)*x2+S27-S23*S25)./sqrt(c16.*c18);
            R15 = ((S25*S4*S6-S28*S4^2-S6*S7+S24*S4*S7)*x2+S28-S24*S25)./sqrt(c17.*c18);
            z4 = S4*fi6([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15])./sqrt(1-S4^2*x2);
            %
            %  i=6
            %
            R1 = ((S1*S16*S5-S5^2*S8-S1*S2+S11*S2*S5)*x2+S8-S11*S16)./sqrt(c4.*c9);
            R2 = ((S1*S20*S5-S5^2*S9-S1*S3+S11*S3*S5)*x2+S9-S11*S20)./sqrt(c4.*c13);
            R3 = ((S1*S23*S5-S10*S5^2-S1*S4+S11*S4*S5)*x2+S10-S11*S23)./sqrt(c4.*c16);
            R4 = ((S1*S26*S5-S12*S5^2-S1*S6+S11*S5*S6)*x2+S12-S11*S26)./sqrt(c4.*c19);
            R5 = ((S1*S27*S5-S13*S5^2-S1*S7+S11*S5*S7)*x2+S13-S11*S27)./sqrt(c4.*c20);
            R6 = ((S2*S20*S5-S14*S5^2-S2*S3+S16*S3*S5)*x2+S14-S16*S20)./sqrt(c9.*c13);
            R7 = ((S2*S23*S5-S15*S5^2-S2*S4+S16*S4*S5)*x2+S15-S16*S23)./sqrt(c9.*c16);
            R8 = ((S2*S26*S5-S17*S5^2-S2*S6+S16*S5*S6)*x2+S17-S16*S26)./sqrt(c9.*c19);
            R9 = ((S2*S27*S5-S18*S5^2-S2*S7+S16*S5*S7)*x2+S18-S16*S27)./sqrt(c9.*c20);
            R10 = ((S23*S3*S5-S19*S5^2-S3*S4+S20*S4*S5)*x2+S19-S20*S23)./sqrt(c13.*c16);
            R11 = ((S26*S3*S5-S21*S5^2-S3*S6+S20*S5*S6)*x2+S21-S20*S26)./sqrt(c13.*c19);
            R12 = ((S27*S3*S5-S22*S5^2-S3*S7+S20*S5*S7)*x2+S22-S20*S27)./sqrt(c13.*c20);
            R13 = ((S26*S4*S5-S24*S5^2-S4*S6+S23*S5*S6)*x2+S24-S23*S26)./sqrt(c16.*c19);
            R14 = ((S27*S4*S5-S25*S5^2-S4*S7+S23*S5*S7)*x2+S25-S23*S27)./sqrt(c16.*c20);
            R15 = ((S27*S5*S6-S28*S5^2-S6*S7+S26*S5*S7)*x2+S28-S26*S27)./sqrt(c19.*c20);
            z5 = S5*fi6([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15])./sqrt(1-S5^2*x2);
            %
            %  i=7
            %
            R1 = ((S1*S17*S6-S6^2*S8-S1*S2+S12*S2*S6)*x2+S8-S12*S17)./sqrt(c5.*c10);
            R2 = ((S1*S21*S6-S6^2*S9-S1*S3+S12*S3*S6)*x2+S9-S12*S21)./sqrt(c5.*c14);
            R3 = ((S1*S24*S6-S10*S6^2-S1*S4+S12*S4*S6)*x2+S10-S12*S24)./sqrt(c5.*c17);
            R4 = ((S1*S26*S6-S11*S6^2-S1*S5+S12*S5*S6)*x2+S11-S12*S26)./sqrt(c5.*c19);
            R5 = ((S1*S28*S6-S13*S6^2-S1*S7+S12*S6*S7)*x2+S13-S12*S28)./sqrt(c5.*c21);
            R6 = ((S2*S21*S6-S14*S6^2-S2*S3+S17*S3*S6)*x2+S14-S17*S21)./sqrt(c10.*c14);
            R7 = ((S2*S24*S6-S15*S6^2-S2*S4+S17*S4*S6)*x2+S15-S17*S24)./sqrt(c10.*c17);
            R8 = ((S2*S26*S6-S16*S6^2-S2*S5+S17*S5*S6)*x2+S16-S17*S26)./sqrt(c10.*c19);
            R9 = ((S2*S28*S6-S18*S6^2-S2*S7+S17*S6*S7)*x2+S18-S17*S28)./sqrt(c10.*c21);
            R10 = ((S24*S3*S6-S19*S6^2-S3*S4+S21*S4*S6)*x2+S19-S21*S24)./sqrt(c14.*c17);
            R11 = ((S26*S3*S6-S20*S6^2-S3*S5+S21*S5*S6)*x2+S20-S21*S26)./sqrt(c14.*c19);
            R12 = ((S28*S3*S6-S22*S6^2-S3*S7+S21*S6*S7)*x2+S22-S21*S28)./sqrt(c14.*c21);
            R13 = ((S26*S4*S6-S23*S6^2-S4*S5+S24*S5*S6)*x2+S23-S24*S26)./sqrt(c17.*c19);
            R14 = ((S28*S4*S6-S25*S6^2-S4*S7+S24*S6*S7)*x2+S25-S24*S28)./sqrt(c17.*c21);
            R15 = ((S28*S5*S6-S27*S6^2-S5*S7+S26*S6*S7)*x2+S27-S26*S28)./sqrt(c19.*c21);
            z6 = S6*fi6([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15])./sqrt(1-S6^2*x2);
            %
            %  i=8
            %
            R1 = ((S1*S18*S7-S7^2*S8-S1*S2+S13*S2*S7)*x2+S8-S13*S18)./sqrt(c6.*c11);
            R2 = ((S1*S22*S7-S7^2*S9-S1*S3+S13*S3*S7)*x2+S9-S13*S22)./sqrt(c6.*c15);
            R3 = ((S1*S25*S7-S10*S7^2-S1*S4+S13*S4*S7)*x2+S10-S13*S25)./sqrt(c6.*c18);
            R4 = ((S1*S27*S7-S11*S7^2-S1*S5+S13*S5*S7)*x2+S11-S13*S27)./sqrt(c6.*c20);
            R5 = ((S1*S28*S7-S12*S7^2-S1*S6+S13*S6*S7)*x2+S12-S13*S28)./sqrt(c6.*c21);
            R6 = ((S2*S22*S7-S14*S7^2-S2*S3+S18*S3*S7)*x2+S14-S18*S22)./sqrt(c11.*c15);
            R7 = ((S2*S25*S7-S15*S7^2-S2*S4+S18*S4*S7)*x2+S15-S18*S25)./sqrt(c11.*c18);
            R8 = ((S2*S27*S7-S16*S7^2-S2*S5+S18*S5*S7)*x2+S16-S18*S27)./sqrt(c11.*c20);
            R9 = ((S2*S28*S7-S17*S7^2-S2*S6+S18*S6*S7)*x2+S17-S18*S28)./sqrt(c11.*c21);
            R10 = ((S25*S3*S7-S19*S7^2-S3*S4+S22*S4*S7)*x2+S19-S22*S25)./sqrt(c15.*c18);
            R11 = ((S27*S3*S7-S20*S7^2-S3*S5+S22*S5*S7)*x2+S20-S22*S27)./sqrt(c15.*c20);
            R12 = ((S28*S3*S7-S21*S7^2-S3*S6+S22*S6*S7)*x2+S21-S22*S28)./sqrt(c15.*c21);
            R13 = ((S27*S4*S7-S23*S7^2-S4*S5+S25*S5*S7)*x2+S23-S25*S27)./sqrt(c18.*c20);
            R14 = ((S28*S4*S7-S24*S7^2-S4*S6+S25*S6*S7)*x2+S24-S25*S28)./sqrt(c18.*c21);
            R15 = ((S28*S5*S7-S26*S7^2-S5*S6+S27*S6*S7)*x2+S26-S27*S28)./sqrt(c20.*c21);
            z7 = S7*fi6([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15])./sqrt(1-S7^2*x2);
            z(i) = wg'*(z1+z2+z3+z4+z5+z6+z7);
        end
    end

    function z = fi10(rr)
        %
        %  This function computes the integral for I_10(S1,...,S45),
        %  where S1 to S45 are the correlation coefficients.  When Si's are
        %  vectors, the function returns a vector output of the integral.
        %
        z = zeros(size(rr,1),1);
        for i=1:size(rr,1)
            S1 = rr(i,1);  S2 = rr(i,2);  S3 = rr(i,3);  S4 = rr(i,4);  S5 = rr(i,5);
            S6 = rr(i,6);  S7 = rr(i,7);  S8 = rr(i,8);  S9 = rr(i,9);  S10 = rr(i,10);
            S11 = rr(i,11);  S12 = rr(i,12);  S13 = rr(i,13);  S14 = rr(i,14);  S15 = rr(i,15);
            S16 = rr(i,16);  S17 = rr(i,17);  S18 = rr(i,18);  S19 = rr(i,19);  S20 = rr(i,20);
            S21 = rr(i,21);  S22 = rr(i,22);  S23 = rr(i,23);  S24 = rr(i,24);  S25 = rr(i,25);
            S26 = rr(i,26);  S27 = rr(i,27);  S28 = rr(i,28);  S29 = rr(i,29);  S30 = rr(i,30);
            S31 = rr(i,31);  S32 = rr(i,32);  S33 = rr(i,33);  S34 = rr(i,34);  S35 = rr(i,35);
            S36 = rr(i,36);  S37 = rr(i,37);  S38 = rr(i,38);  S39 = rr(i,39);  S40 = rr(i,40);
            S41 = rr(i,41);  S42 = rr(i,42);  S43 = rr(i,43);  S44 = rr(i,44);  S45 = rr(i,45);
            c1 = 1-S10^2-(S1^2+S2^2-2*S1*S2*S10)*x2;
            c2 = 1-S11^2-(S1^2+S3^2-2*S1*S3*S11)*x2;
            c3 = 1-S12^2-(S1^2+S4^2-2*S1*S4*S12)*x2;
            c4 = 1-S13^2-(S1^2+S5^2-2*S1*S5*S13)*x2;
            c5 = 1-S14^2-(S1^2+S6^2-2*S1*S6*S14)*x2;
            c6 = 1-S15^2-(S1^2+S7^2-2*S1*S7*S15)*x2;
            c7 = 1-S16^2-(S1^2+S8^2-2*S1*S8*S16)*x2;
            c8 = 1-S17^2-(S1^2+S9^2-2*S1*S9*S17)*x2;
            c9 = 1-S18^2-(S2^2+S3^2-2*S2*S3*S18)*x2;
            c10 = 1-S19^2-(S2^2+S4^2-2*S2*S4*S19)*x2;
            c11 = 1-S20^2-(S2^2+S5^2-2*S2*S5*S20)*x2;
            c12 = 1-S21^2-(S2^2+S6^2-2*S2*S6*S21)*x2;
            c13 = 1-S22^2-(S2^2+S7^2-2*S2*S7*S22)*x2;
            c14 = 1-S23^2-(S2^2+S8^2-2*S2*S8*S23)*x2;
            c15 = 1-S24^2-(S2^2+S9^2-2*S2*S9*S24)*x2;
            c16 = 1-S25^2-(S3^2+S4^2-2*S3*S4*S25)*x2;
            c17 = 1-S26^2-(S3^2+S5^2-2*S3*S5*S26)*x2;
            c18 = 1-S27^2-(S3^2+S6^2-2*S3*S6*S27)*x2;
            c19 = 1-S28^2-(S3^2+S7^2-2*S3*S7*S28)*x2;
            c20 = 1-S29^2-(S3^2+S8^2-2*S3*S8*S29)*x2;
            c21 = 1-S30^2-(S3^2+S9^2-2*S3*S9*S30)*x2;
            c22 = 1-S31^2-(S4^2+S5^2-2*S4*S5*S31)*x2;
            c23 = 1-S32^2-(S4^2+S6^2-2*S4*S6*S32)*x2;
            c24 = 1-S33^2-(S4^2+S7^2-2*S4*S7*S33)*x2;
            c25 = 1-S34^2-(S4^2+S8^2-2*S4*S8*S34)*x2;
            c26 = 1-S35^2-(S4^2+S9^2-2*S4*S9*S35)*x2;
            c27 = 1-S36^2-(S5^2+S6^2-2*S5*S6*S36)*x2;
            c28 = 1-S37^2-(S5^2+S7^2-2*S5*S7*S37)*x2;
            c29 = 1-S38^2-(S5^2+S8^2-2*S5*S8*S38)*x2;
            c30 = 1-S39^2-(S5^2+S9^2-2*S5*S9*S39)*x2;
            c31 = 1-S40^2-(S6^2+S7^2-2*S6*S7*S40)*x2;
            c32 = 1-S41^2-(S6^2+S8^2-2*S6*S8*S41)*x2;
            c33 = 1-S42^2-(S6^2+S9^2-2*S6*S9*S42)*x2;
            c34 = 1-S43^2-(S7^2+S8^2-2*S7*S8*S43)*x2;
            c35 = 1-S44^2-(S7^2+S9^2-2*S7*S9*S44)*x2;
            c36 = 1-S45^2-(S8^2+S9^2-2*S8*S9*S45)*x2;
            %
            %  i=2
            %
            R1 =  ((S1*S11*S2-S1^2*S18-S2*S3+S1*S10*S3)*x2+S18-S10*S11)./sqrt(c1.*c2);
            R2 =  ((S1*S12*S2-S1^2*S19-S2*S4+S1*S10*S4)*x2+S19-S10*S12)./sqrt(c1.*c3);
            R3 =  ((S1*S13*S2-S1^2*S20-S2*S5+S1*S10*S5)*x2+S20-S10*S13)./sqrt(c1.*c4);
            R4 =  ((S1*S14*S2-S1^2*S21-S2*S6+S1*S10*S6)*x2+S21-S10*S14)./sqrt(c1.*c5);
            R5 =  ((S1*S15*S2-S1^2*S22-S2*S7+S1*S10*S7)*x2+S22-S10*S15)./sqrt(c1.*c6);
            R6 =  ((S1*S16*S2-S1^2*S23-S2*S8+S1*S10*S8)*x2+S23-S10*S16)./sqrt(c1.*c7);
            R7 =  ((S1*S17*S2-S1^2*S24-S2*S9+S1*S10*S9)*x2+S24-S10*S17)./sqrt(c1.*c8);
            R8 =  ((S1*S12*S3-S1^2*S25-S3*S4+S1*S11*S4)*x2+S25-S11*S12)./sqrt(c2.*c3);
            R9 =  ((S1*S13*S3-S1^2*S26-S3*S5+S1*S11*S5)*x2+S26-S11*S13)./sqrt(c2.*c4);
            R10 = ((S1*S14*S3-S1^2*S27-S3*S6+S1*S11*S6)*x2+S27-S11*S14)./sqrt(c2.*c5);
            R11 = ((S1*S15*S3-S1^2*S28-S3*S7+S1*S11*S7)*x2+S28-S11*S15)./sqrt(c2.*c6);
            R12 = ((S1*S16*S3-S1^2*S29-S3*S8+S1*S11*S8)*x2+S29-S11*S16)./sqrt(c2.*c7);
            R13 = ((S1*S17*S3-S1^2*S30-S3*S9+S1*S11*S9)*x2+S30-S11*S17)./sqrt(c2.*c8);
            R14 = ((S1*S13*S4-S1^2*S31-S4*S5+S1*S12*S5)*x2+S31-S12*S13)./sqrt(c3.*c4);
            R15 = ((S1*S14*S4-S1^2*S32-S4*S6+S1*S12*S6)*x2+S32-S12*S14)./sqrt(c3.*c5);
            R16 = ((S1*S15*S4-S1^2*S33-S4*S7+S1*S12*S7)*x2+S33-S12*S15)./sqrt(c3.*c6);
            R17 = ((S1*S16*S4-S1^2*S34-S4*S8+S1*S12*S8)*x2+S34-S12*S16)./sqrt(c3.*c7);
            R18 = ((S1*S17*S4-S1^2*S35-S4*S9+S1*S12*S9)*x2+S35-S12*S17)./sqrt(c3.*c8);
            R19 = ((S1*S14*S5-S1^2*S36-S5*S6+S1*S13*S6)*x2+S36-S13*S14)./sqrt(c4.*c5);
            R20 = ((S1*S15*S5-S1^2*S37-S5*S7+S1*S13*S7)*x2+S37-S13*S15)./sqrt(c4.*c6);
            R21 = ((S1*S16*S5-S1^2*S38-S5*S8+S1*S13*S8)*x2+S38-S13*S16)./sqrt(c4.*c7);
            R22 = ((S1*S17*S5-S1^2*S39-S5*S9+S1*S13*S9)*x2+S39-S13*S17)./sqrt(c4.*c8);
            R23 = ((S1*S15*S6-S1^2*S40-S6*S7+S1*S14*S7)*x2+S40-S14*S15)./sqrt(c5.*c6);
            R24 = ((S1*S16*S6-S1^2*S41-S6*S8+S1*S14*S8)*x2+S41-S14*S16)./sqrt(c5.*c7);
            R25 = ((S1*S17*S6-S1^2*S42-S6*S9+S1*S14*S9)*x2+S42-S14*S17)./sqrt(c5.*c8);
            R26 = ((S1*S16*S7-S1^2*S43-S7*S8+S1*S15*S8)*x2+S43-S15*S16)./sqrt(c6.*c7);
            R27 = ((S1*S17*S7-S1^2*S44-S7*S9+S1*S15*S9)*x2+S44-S15*S17)./sqrt(c6.*c8);
            R28 = ((S1*S17*S8-S1^2*S45-S8*S9+S1*S16*S9)*x2+S45-S16*S17)./sqrt(c7.*c8);
            z1 = S1*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S1^2*x2);
            %
            %  i=3
            %
            R1 =  ((S1*S18*S2-S11*S2^2-S1*S3+S10*S2*S3)*x2+S11-S10*S18)./sqrt(c1.*c9);
            R2 =  ((S1*S19*S2-S12*S2^2-S1*S4+S10*S2*S4)*x2+S12-S10*S19)./sqrt(c1.*c10);
            R3 =  ((S1*S2*S20-S13*S2^2-S1*S5+S10*S2*S5)*x2+S13-S10*S20)./sqrt(c1.*c11);
            R4 =  ((S1*S2*S21-S14*S2^2-S1*S6+S10*S2*S6)*x2+S14-S10*S21)./sqrt(c1.*c12);
            R5 =  ((S1*S2*S22-S15*S2^2-S1*S7+S10*S2*S7)*x2+S15-S10*S22)./sqrt(c1.*c13);
            R6 =  ((S1*S2*S23-S16*S2^2-S1*S8+S10*S2*S8)*x2+S16-S10*S23)./sqrt(c1.*c14);
            R7 =  ((S1*S2*S24-S17*S2^2-S1*S9+S10*S2*S9)*x2+S17-S10*S24)./sqrt(c1.*c15);
            R8 =  ((S19*S2*S3-S2^2*S25-S3*S4+S18*S2*S4)*x2+S25-S18*S19)./sqrt(c9.*c10);
            R9 =  ((S2*S20*S3-S2^2*S26-S3*S5+S18*S2*S5)*x2+S26-S18*S20)./sqrt(c9.*c11);
            R10 = ((S2*S21*S3-S2^2*S27-S3*S6+S18*S2*S6)*x2+S27-S18*S21)./sqrt(c9.*c12);
            R11 = ((S2*S22*S3-S2^2*S28-S3*S7+S18*S2*S7)*x2+S28-S18*S22)./sqrt(c9.*c13);
            R12 = ((S2*S23*S3-S2^2*S29-S3*S8+S18*S2*S8)*x2+S29-S18*S23)./sqrt(c9.*c14);
            R13 = ((S2*S24*S3-S2^2*S30-S3*S9+S18*S2*S9)*x2+S30-S18*S24)./sqrt(c9.*c15);
            R14 = ((S2*S20*S4-S2^2*S31-S4*S5+S19*S2*S5)*x2+S31-S19*S20)./sqrt(c10.*c11);
            R15 = ((S2*S21*S4-S2^2*S32-S4*S6+S19*S2*S6)*x2+S32-S19*S21)./sqrt(c10.*c12);
            R16 = ((S2*S22*S4-S2^2*S33-S4*S7+S19*S2*S7)*x2+S33-S19*S22)./sqrt(c10.*c13);
            R17 = ((S2*S23*S4-S2^2*S34-S4*S8+S19*S2*S8)*x2+S34-S19*S23)./sqrt(c10.*c14);
            R18 = ((S2*S24*S4-S2^2*S35-S4*S9+S19*S2*S9)*x2+S35-S19*S24)./sqrt(c10.*c15);
            R19 = ((S2*S21*S5-S2^2*S36-S5*S6+S2*S20*S6)*x2+S36-S20*S21)./sqrt(c11.*c12);
            R20 = ((S2*S22*S5-S2^2*S37-S5*S7+S2*S20*S7)*x2+S37-S20*S22)./sqrt(c11.*c13);
            R21 = ((S2*S23*S5-S2^2*S38-S5*S8+S2*S20*S8)*x2+S38-S20*S23)./sqrt(c11.*c14);
            R22 = ((S2*S24*S5-S2^2*S39-S5*S9+S2*S20*S9)*x2+S39-S20*S24)./sqrt(c11.*c15);
            R23 = ((S2*S22*S6-S2^2*S40-S6*S7+S2*S21*S7)*x2+S40-S21*S22)./sqrt(c12.*c13);
            R24 = ((S2*S23*S6-S2^2*S41-S6*S8+S2*S21*S8)*x2+S41-S21*S23)./sqrt(c12.*c14);
            R25 = ((S2*S24*S6-S2^2*S42-S6*S9+S2*S21*S9)*x2+S42-S21*S24)./sqrt(c12.*c15);
            R26 = ((S2*S23*S7-S2^2*S43-S7*S8+S2*S22*S8)*x2+S43-S22*S23)./sqrt(c13.*c14);
            R27 = ((S2*S24*S7-S2^2*S44-S7*S9+S2*S22*S9)*x2+S44-S22*S24)./sqrt(c13.*c15);
            R28 = ((S2*S24*S8-S2^2*S45-S8*S9+S2*S23*S9)*x2+S45-S23*S24)./sqrt(c14.*c15);
            z2 = S2*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S2^2*x2);
            %
            %  i=4
            %
            R1 =  ((S1*S18*S3-S10*S3^2-S1*S2+S11*S2*S3)*x2+S10-S11*S18)./sqrt(c2.*c9);
            R2 =  ((S1*S25*S3-S12*S3^2-S1*S4+S11*S3*S4)*x2+S12-S11*S25)./sqrt(c2.*c16);
            R3 =  ((S1*S26*S3-S13*S3^2-S1*S5+S11*S3*S5)*x2+S13-S11*S26)./sqrt(c2.*c17);
            R4 =  ((S1*S27*S3-S14*S3^2-S1*S6+S11*S3*S6)*x2+S14-S11*S27)./sqrt(c2.*c18);
            R5 =  ((S1*S28*S3-S15*S3^2-S1*S7+S11*S3*S7)*x2+S15-S11*S28)./sqrt(c2.*c19);
            R6 =  ((S1*S29*S3-S16*S3^2-S1*S8+S11*S3*S8)*x2+S16-S11*S29)./sqrt(c2.*c20);
            R7 =  ((S1*S3*S30-S17*S3^2-S1*S9+S11*S3*S9)*x2+S17-S11*S30)./sqrt(c2.*c21);
            R8 =  ((S2*S25*S3-S19*S3^2-S2*S4+S18*S3*S4)*x2+S19-S18*S25)./sqrt(c9.*c16);
            R9 =  ((S2*S26*S3-S20*S3^2-S2*S5+S18*S3*S5)*x2+S20-S18*S26)./sqrt(c9.*c17);
            R10 = ((S2*S27*S3-S21*S3^2-S2*S6+S18*S3*S6)*x2+S21-S18*S27)./sqrt(c9.*c18);
            R11 = ((S2*S28*S3-S22*S3^2-S2*S7+S18*S3*S7)*x2+S22-S18*S28)./sqrt(c9.*c19);
            R12 = ((S2*S29*S3-S23*S3^2-S2*S8+S18*S3*S8)*x2+S23-S18*S29)./sqrt(c9.*c20);
            R13 = ((S2*S3*S30-S24*S3^2-S2*S9+S18*S3*S9)*x2+S24-S18*S30)./sqrt(c9.*c21);
            R14 = ((S26*S3*S4-S3^2*S31-S4*S5+S25*S3*S5)*x2+S31-S25*S26)./sqrt(c16.*c17);
            R15 = ((S27*S3*S4-S3^2*S32-S4*S6+S25*S3*S6)*x2+S32-S25*S27)./sqrt(c16.*c18);
            R16 = ((S28*S3*S4-S3^2*S33-S4*S7+S25*S3*S7)*x2+S33-S25*S28)./sqrt(c16.*c19);
            R17 = ((S29*S3*S4-S3^2*S34-S4*S8+S25*S3*S8)*x2+S34-S25*S29)./sqrt(c16.*c20);
            R18 = ((S3*S30*S4-S3^2*S35-S4*S9+S25*S3*S9)*x2+S35-S25*S30)./sqrt(c16.*c21);
            R19 = ((S27*S3*S5-S3^2*S36-S5*S6+S26*S3*S6)*x2+S36-S26*S27)./sqrt(c17.*c18);
            R20 = ((S28*S3*S5-S3^2*S37-S5*S7+S26*S3*S7)*x2+S37-S26*S28)./sqrt(c17.*c19);
            R21 = ((S29*S3*S5-S3^2*S38-S5*S8+S26*S3*S8)*x2+S38-S26*S29)./sqrt(c17.*c20);
            R22 = ((S3*S30*S5-S3^2*S39-S5*S9+S26*S3*S9)*x2+S39-S26*S30)./sqrt(c17.*c21);
            R23 = ((S28*S3*S6-S3^2*S40-S6*S7+S27*S3*S7)*x2+S40-S27*S28)./sqrt(c18.*c19);
            R24 = ((S29*S3*S6-S3^2*S41-S6*S8+S27*S3*S8)*x2+S41-S27*S29)./sqrt(c18.*c20);
            R25 = ((S3*S30*S6-S3^2*S42-S6*S9+S27*S3*S9)*x2+S42-S27*S30)./sqrt(c18.*c21);
            R26 = ((S29*S3*S7-S3^2*S43-S7*S8+S28*S3*S8)*x2+S43-S28*S29)./sqrt(c19.*c20);
            R27 = ((S3*S30*S7-S3^2*S44-S7*S9+S28*S3*S9)*x2+S44-S28*S30)./sqrt(c19.*c21);
            R28 = ((S3*S30*S8-S3^2*S45-S8*S9+S29*S3*S9)*x2+S45-S29*S30)./sqrt(c20.*c21);
            z3 = S3*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S3^2*x2);
            %
            %  i=5
            %
            R1 =  ((S1*S19*S4-S10*S4^2-S1*S2+S12*S2*S4)*x2+S10-S12*S19)./sqrt(c3.*c10);
            R2 =  ((S1*S25*S4-S11*S4^2-S1*S3+S12*S3*S4)*x2+S11-S12*S25)./sqrt(c3.*c16);
            R3 =  ((S1*S31*S4-S13*S4^2-S1*S5+S12*S4*S5)*x2+S13-S12*S31)./sqrt(c3.*c22);
            R4 =  ((S1*S32*S4-S14*S4^2-S1*S6+S12*S4*S6)*x2+S14-S12*S32)./sqrt(c3.*c23);
            R5 =  ((S1*S33*S4-S15*S4^2-S1*S7+S12*S4*S7)*x2+S15-S12*S33)./sqrt(c3.*c24);
            R6 =  ((S1*S34*S4-S16*S4^2-S1*S8+S12*S4*S8)*x2+S16-S12*S34)./sqrt(c3.*c25);
            R7 =  ((S1*S35*S4-S17*S4^2-S1*S9+S12*S4*S9)*x2+S17-S12*S35)./sqrt(c3.*c26);
            R8 =  ((S2*S25*S4-S18*S4^2-S2*S3+S19*S3*S4)*x2+S18-S19*S25)./sqrt(c10.*c16);
            R9 =  ((S2*S31*S4-S20*S4^2-S2*S5+S19*S4*S5)*x2+S20-S19*S31)./sqrt(c10.*c22);
            R10 = ((S2*S32*S4-S21*S4^2-S2*S6+S19*S4*S6)*x2+S21-S19*S32)./sqrt(c10.*c23);
            R11 = ((S2*S33*S4-S22*S4^2-S2*S7+S19*S4*S7)*x2+S22-S19*S33)./sqrt(c10.*c24);
            R12 = ((S2*S34*S4-S23*S4^2-S2*S8+S19*S4*S8)*x2+S23-S19*S34)./sqrt(c10.*c25);
            R13 = ((S2*S35*S4-S24*S4^2-S2*S9+S19*S4*S9)*x2+S24-S19*S35)./sqrt(c10.*c26);
            R14 = ((S3*S31*S4-S26*S4^2-S3*S5+S25*S4*S5)*x2+S26-S25*S31)./sqrt(c16.*c22);
            R15 = ((S3*S32*S4-S27*S4^2-S3*S6+S25*S4*S6)*x2+S27-S25*S32)./sqrt(c16.*c23);
            R16 = ((S3*S33*S4-S28*S4^2-S3*S7+S25*S4*S7)*x2+S28-S25*S33)./sqrt(c16.*c24);
            R17 = ((S3*S34*S4-S29*S4^2-S3*S8+S25*S4*S8)*x2+S29-S25*S34)./sqrt(c16.*c25);
            R18 = ((S3*S35*S4-S30*S4^2-S3*S9+S25*S4*S9)*x2+S30-S25*S35)./sqrt(c16.*c26);
            R19 = ((S32*S4*S5-S36*S4^2-S5*S6+S31*S4*S6)*x2+S36-S31*S32)./sqrt(c22.*c23);
            R20 = ((S33*S4*S5-S37*S4^2-S5*S7+S31*S4*S7)*x2+S37-S31*S33)./sqrt(c22.*c24);
            R21 = ((S34*S4*S5-S38*S4^2-S5*S8+S31*S4*S8)*x2+S38-S31*S34)./sqrt(c22.*c25);
            R22 = ((S35*S4*S5-S39*S4^2-S5*S9+S31*S4*S9)*x2+S39-S31*S35)./sqrt(c22.*c26);
            R23 = ((S33*S4*S6-S4^2*S40-S6*S7+S32*S4*S7)*x2+S40-S32*S33)./sqrt(c23.*c24);
            R24 = ((S34*S4*S6-S4^2*S41-S6*S8+S32*S4*S8)*x2+S41-S32*S34)./sqrt(c23.*c25);
            R25 = ((S35*S4*S6-S4^2*S42-S6*S9+S32*S4*S9)*x2+S42-S32*S35)./sqrt(c23.*c26);
            R26 = ((S34*S4*S7-S4^2*S43-S7*S8+S33*S4*S8)*x2+S43-S33*S34)./sqrt(c24.*c25);
            R27 = ((S35*S4*S7-S4^2*S44-S7*S9+S33*S4*S9)*x2+S44-S33*S35)./sqrt(c24.*c26);
            R28 = ((S35*S4*S8-S4^2*S45-S8*S9+S34*S4*S9)*x2+S45-S34*S35)./sqrt(c25.*c26);
            z4 = S4*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S4^2*x2);
            %
            %  i=6
            %
            R1 =  ((S1*S20*S5-S10*S5^2-S1*S2+S13*S2*S5)*x2+S10-S13*S20)./sqrt(c4.*c11);
            R2 =  ((S1*S26*S5-S11*S5^2-S1*S3+S13*S3*S5)*x2+S11-S13*S26)./sqrt(c4.*c17);
            R3 =  ((S1*S31*S5-S12*S5^2-S1*S4+S13*S4*S5)*x2+S12-S13*S31)./sqrt(c4.*c22);
            R4 =  ((S1*S36*S5-S14*S5^2-S1*S6+S13*S5*S6)*x2+S14-S13*S36)./sqrt(c4.*c27);
            R5 =  ((S1*S37*S5-S15*S5^2-S1*S7+S13*S5*S7)*x2+S15-S13*S37)./sqrt(c4.*c28);
            R6 =  ((S1*S38*S5-S16*S5^2-S1*S8+S13*S5*S8)*x2+S16-S13*S38)./sqrt(c4.*c29);
            R7 =  ((S1*S39*S5-S17*S5^2-S1*S9+S13*S5*S9)*x2+S17-S13*S39)./sqrt(c4.*c30);
            R8 =  ((S2*S26*S5-S18*S5^2-S2*S3+S20*S3*S5)*x2+S18-S20*S26)./sqrt(c11.*c17);
            R9 =  ((S2*S31*S5-S19*S5^2-S2*S4+S20*S4*S5)*x2+S19-S20*S31)./sqrt(c11.*c22);
            R10 = ((S2*S36*S5-S21*S5^2-S2*S6+S20*S5*S6)*x2+S21-S20*S36)./sqrt(c11.*c27);
            R11 = ((S2*S37*S5-S22*S5^2-S2*S7+S20*S5*S7)*x2+S22-S20*S37)./sqrt(c11.*c28);
            R12 = ((S2*S38*S5-S23*S5^2-S2*S8+S20*S5*S8)*x2+S23-S20*S38)./sqrt(c11.*c29);
            R13 = ((S2*S39*S5-S24*S5^2-S2*S9+S20*S5*S9)*x2+S24-S20*S39)./sqrt(c11.*c30);
            R14 = ((S3*S31*S5-S25*S5^2-S3*S4+S26*S4*S5)*x2+S25-S26*S31)./sqrt(c17.*c22);
            R15 = ((S3*S36*S5-S27*S5^2-S3*S6+S26*S5*S6)*x2+S27-S26*S36)./sqrt(c17.*c27);
            R16 = ((S3*S37*S5-S28*S5^2-S3*S7+S26*S5*S7)*x2+S28-S26*S37)./sqrt(c17.*c28);
            R17 = ((S3*S38*S5-S29*S5^2-S3*S8+S26*S5*S8)*x2+S29-S26*S38)./sqrt(c17.*c29);
            R18 = ((S3*S39*S5-S30*S5^2-S3*S9+S26*S5*S9)*x2+S30-S26*S39)./sqrt(c17.*c30);
            R19 = ((S36*S4*S5-S32*S5^2-S4*S6+S31*S5*S6)*x2+S32-S31*S36)./sqrt(c22.*c27);
            R20 = ((S37*S4*S5-S33*S5^2-S4*S7+S31*S5*S7)*x2+S33-S31*S37)./sqrt(c22.*c28);
            R21 = ((S38*S4*S5-S34*S5^2-S4*S8+S31*S5*S8)*x2+S34-S31*S38)./sqrt(c22.*c29);
            R22 = ((S39*S4*S5-S35*S5^2-S4*S9+S31*S5*S9)*x2+S35-S31*S39)./sqrt(c22.*c30);
            R23 = ((S37*S5*S6-S40*S5^2-S6*S7+S36*S5*S7)*x2+S40-S36*S37)./sqrt(c27.*c28);
            R24 = ((S38*S5*S6-S41*S5^2-S6*S8+S36*S5*S8)*x2+S41-S36*S38)./sqrt(c27.*c29);
            R25 = ((S39*S5*S6-S42*S5^2-S6*S9+S36*S5*S9)*x2+S42-S36*S39)./sqrt(c27.*c30);
            R26 = ((S38*S5*S7-S43*S5^2-S7*S8+S37*S5*S8)*x2+S43-S37*S38)./sqrt(c28.*c29);
            R27 = ((S39*S5*S7-S44*S5^2-S7*S9+S37*S5*S9)*x2+S44-S37*S39)./sqrt(c28.*c30);
            R28 = ((S39*S5*S8-S45*S5^2-S8*S9+S38*S5*S9)*x2+S45-S38*S39)./sqrt(c29.*c30);
            z5 = S5*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S5^2*x2);
            %
            %  i=7
            %
            R1 =  ((S1*S21*S6-S10*S6^2-S1*S2+S14*S2*S6)*x2+S10-S14*S21)./sqrt(c5.*c12);
            R2 =  ((S1*S27*S6-S11*S6^2-S1*S3+S14*S3*S6)*x2+S11-S14*S27)./sqrt(c5.*c18);
            R3 =  ((S1*S32*S6-S12*S6^2-S1*S4+S14*S4*S6)*x2+S12-S14*S32)./sqrt(c5.*c23);
            R4 =  ((S1*S36*S6-S13*S6^2-S1*S5+S14*S5*S6)*x2+S13-S14*S36)./sqrt(c5.*c27);
            R5 =  ((S1*S40*S6-S15*S6^2-S1*S7+S14*S6*S7)*x2+S15-S14*S40)./sqrt(c5.*c31);
            R6 =  ((S1*S41*S6-S16*S6^2-S1*S8+S14*S6*S8)*x2+S16-S14*S41)./sqrt(c5.*c32);
            R7 =  ((S1*S42*S6-S17*S6^2-S1*S9+S14*S6*S9)*x2+S17-S14*S42)./sqrt(c5.*c33);
            R8 =  ((S2*S27*S6-S18*S6^2-S2*S3+S21*S3*S6)*x2+S18-S21*S27)./sqrt(c12.*c18);
            R9 =  ((S2*S32*S6-S19*S6^2-S2*S4+S21*S4*S6)*x2+S19-S21*S32)./sqrt(c12.*c23);
            R10 = ((S2*S36*S6-S20*S6^2-S2*S5+S21*S5*S6)*x2+S20-S21*S36)./sqrt(c12.*c27);
            R11 = ((S2*S40*S6-S22*S6^2-S2*S7+S21*S6*S7)*x2+S22-S21*S40)./sqrt(c12.*c31);
            R12 = ((S2*S41*S6-S23*S6^2-S2*S8+S21*S6*S8)*x2+S23-S21*S41)./sqrt(c12.*c32);
            R13 = ((S2*S42*S6-S24*S6^2-S2*S9+S21*S6*S9)*x2+S24-S21*S42)./sqrt(c12.*c33);
            R14 = ((S3*S32*S6-S25*S6^2-S3*S4+S27*S4*S6)*x2+S25-S27*S32)./sqrt(c18.*c23);
            R15 = ((S3*S36*S6-S26*S6^2-S3*S5+S27*S5*S6)*x2+S26-S27*S36)./sqrt(c18.*c27);
            R16 = ((S3*S40*S6-S28*S6^2-S3*S7+S27*S6*S7)*x2+S28-S27*S40)./sqrt(c18.*c31);
            R17 = ((S3*S41*S6-S29*S6^2-S3*S8+S27*S6*S8)*x2+S29-S27*S41)./sqrt(c18.*c32);
            R18 = ((S3*S42*S6-S30*S6^2-S3*S9+S27*S6*S9)*x2+S30-S27*S42)./sqrt(c18.*c33);
            R19 = ((S36*S4*S6-S31*S6^2-S4*S5+S32*S5*S6)*x2+S31-S32*S36)./sqrt(c23.*c27);
            R20 = ((S4*S40*S6-S33*S6^2-S4*S7+S32*S6*S7)*x2+S33-S32*S40)./sqrt(c23.*c31);
            R21 = ((S4*S41*S6-S34*S6^2-S4*S8+S32*S6*S8)*x2+S34-S32*S41)./sqrt(c23.*c32);
            R22 = ((S4*S42*S6-S35*S6^2-S4*S9+S32*S6*S9)*x2+S35-S32*S42)./sqrt(c23.*c33);
            R23 = ((S40*S5*S6-S37*S6^2-S5*S7+S36*S6*S7)*x2+S37-S36*S40)./sqrt(c27.*c31);
            R24 = ((S41*S5*S6-S38*S6^2-S5*S8+S36*S6*S8)*x2+S38-S36*S41)./sqrt(c27.*c32);
            R25 = ((S42*S5*S6-S39*S6^2-S5*S9+S36*S6*S9)*x2+S39-S36*S42)./sqrt(c27.*c33);
            R26 = ((S41*S6*S7-S43*S6^2-S7*S8+S40*S6*S8)*x2+S43-S40*S41)./sqrt(c31.*c32);
            R27 = ((S42*S6*S7-S44*S6^2-S7*S9+S40*S6*S9)*x2+S44-S40*S42)./sqrt(c31.*c33);
            R28 = ((S42*S6*S8-S45*S6^2-S8*S9+S41*S6*S9)*x2+S45-S41*S42)./sqrt(c32.*c33);
            z6 = S6*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S6^2*x2);
            %
            %  i=8
            %
            R1 =  ((S1*S22*S7-S10*S7^2-S1*S2+S15*S2*S7)*x2+S10-S15*S22)./sqrt(c6.*c13);
            R2 =  ((S1*S28*S7-S11*S7^2-S1*S3+S15*S3*S7)*x2+S11-S15*S28)./sqrt(c6.*c19);
            R3 =  ((S1*S33*S7-S12*S7^2-S1*S4+S15*S4*S7)*x2+S12-S15*S33)./sqrt(c6.*c24);
            R4 =  ((S1*S37*S7-S13*S7^2-S1*S5+S15*S5*S7)*x2+S13-S15*S37)./sqrt(c6.*c28);
            R5 =  ((S1*S40*S7-S14*S7^2-S1*S6+S15*S6*S7)*x2+S14-S15*S40)./sqrt(c6.*c31);
            R6 =  ((S1*S43*S7-S16*S7^2-S1*S8+S15*S7*S8)*x2+S16-S15*S43)./sqrt(c6.*c34);
            R7 =  ((S1*S44*S7-S17*S7^2-S1*S9+S15*S7*S9)*x2+S17-S15*S44)./sqrt(c6.*c35);
            R8 =  ((S2*S28*S7-S18*S7^2-S2*S3+S22*S3*S7)*x2+S18-S22*S28)./sqrt(c13.*c19);
            R9 =  ((S2*S33*S7-S19*S7^2-S2*S4+S22*S4*S7)*x2+S19-S22*S33)./sqrt(c13.*c24);
            R10 = ((S2*S37*S7-S20*S7^2-S2*S5+S22*S5*S7)*x2+S20-S22*S37)./sqrt(c13.*c28);
            R11 = ((S2*S40*S7-S21*S7^2-S2*S6+S22*S6*S7)*x2+S21-S22*S40)./sqrt(c13.*c31);
            R12 = ((S2*S43*S7-S23*S7^2-S2*S8+S22*S7*S8)*x2+S23-S22*S43)./sqrt(c13.*c34);
            R13 = ((S2*S44*S7-S24*S7^2-S2*S9+S22*S7*S9)*x2+S24-S22*S44)./sqrt(c13.*c35);
            R14 = ((S3*S33*S7-S25*S7^2-S3*S4+S28*S4*S7)*x2+S25-S28*S33)./sqrt(c19.*c24);
            R15 = ((S3*S37*S7-S26*S7^2-S3*S5+S28*S5*S7)*x2+S26-S28*S37)./sqrt(c19.*c28);
            R16 = ((S3*S40*S7-S27*S7^2-S3*S6+S28*S6*S7)*x2+S27-S28*S40)./sqrt(c19.*c31);
            R17 = ((S3*S43*S7-S29*S7^2-S3*S8+S28*S7*S8)*x2+S29-S28*S43)./sqrt(c19.*c34);
            R18 = ((S3*S44*S7-S30*S7^2-S3*S9+S28*S7*S9)*x2+S30-S28*S44)./sqrt(c19.*c35);
            R19 = ((S37*S4*S7-S31*S7^2-S4*S5+S33*S5*S7)*x2+S31-S33*S37)./sqrt(c24.*c28);
            R20 = ((S4*S40*S7-S32*S7^2-S4*S6+S33*S6*S7)*x2+S32-S33*S40)./sqrt(c24.*c31);
            R21 = ((S4*S43*S7-S34*S7^2-S4*S8+S33*S7*S8)*x2+S34-S33*S43)./sqrt(c24.*c34);
            R22 = ((S4*S44*S7-S35*S7^2-S4*S9+S33*S7*S9)*x2+S35-S33*S44)./sqrt(c24.*c35);
            R23 = ((S40*S5*S7-S36*S7^2-S5*S6+S37*S6*S7)*x2+S36-S37*S40)./sqrt(c28.*c31);
            R24 = ((S43*S5*S7-S38*S7^2-S5*S8+S37*S7*S8)*x2+S38-S37*S43)./sqrt(c28.*c34);
            R25 = ((S44*S5*S7-S39*S7^2-S5*S9+S37*S7*S9)*x2+S39-S37*S44)./sqrt(c28.*c35);
            R26 = ((S43*S6*S7-S41*S7^2-S6*S8+S40*S7*S8)*x2+S41-S40*S43)./sqrt(c31.*c34);
            R27 = ((S44*S6*S7-S42*S7^2-S6*S9+S40*S7*S9)*x2+S42-S40*S44)./sqrt(c31.*c35);
            R28 = ((S44*S7*S8-S45*S7^2-S8*S9+S43*S7*S9)*x2+S45-S43*S44)./sqrt(c34.*c35);
            z7 = S7*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S7^2*x2);
            %
            %  i=9
            %
            R1 =  ((S1*S23*S8-S10*S8^2-S1*S2+S16*S2*S8)*x2+S10-S16*S23)./sqrt(c7.*c14);
            R2 =  ((S1*S29*S8-S11*S8^2-S1*S3+S16*S3*S8)*x2+S11-S16*S29)./sqrt(c7.*c20);
            R3 =  ((S1*S34*S8-S12*S8^2-S1*S4+S16*S4*S8)*x2+S12-S16*S34)./sqrt(c7.*c25);
            R4 =  ((S1*S38*S8-S13*S8^2-S1*S5+S16*S5*S8)*x2+S13-S16*S38)./sqrt(c7.*c20);
            R5 =  ((S1*S41*S8-S14*S8^2-S1*S6+S16*S6*S8)*x2+S14-S16*S41)./sqrt(c7.*c25);
            R6 =  ((S1*S43*S8-S15*S8^2-S1*S7+S16*S7*S8)*x2+S15-S16*S43)./sqrt(c7.*c29);
            R7 =  ((S1*S45*S8-S17*S8^2-S1*S9+S16*S8*S9)*x2+S17-S16*S45)./sqrt(c7.*c32);
            R8 =  ((S2*S29*S8-S18*S8^2-S2*S3+S23*S3*S8)*x2+S18-S23*S29)./sqrt(c14.*c34);
            R9 =  ((S2*S34*S8-S19*S8^2-S2*S4+S23*S4*S8)*x2+S19-S23*S34)./sqrt(c14.*c36);
            R10 = ((S2*S38*S8-S20*S8^2-S2*S5+S23*S5*S8)*x2+S20-S23*S38)./sqrt(c14.*c20);
            R11 = ((S2*S41*S8-S21*S8^2-S2*S6+S23*S6*S8)*x2+S21-S23*S41)./sqrt(c14.*c25);
            R12 = ((S2*S43*S8-S22*S8^2-S2*S7+S23*S7*S8)*x2+S22-S23*S43)./sqrt(c14.*c29);
            R13 = ((S2*S45*S8-S24*S8^2-S2*S9+S23*S8*S9)*x2+S24-S23*S45)./sqrt(c14.*c32);
            R14 = ((S3*S34*S8-S25*S8^2-S3*S4+S29*S4*S8)*x2+S25-S29*S34)./sqrt(c20.*c34);
            R15 = ((S3*S38*S8-S26*S8^2-S3*S5+S29*S5*S8)*x2+S26-S29*S38)./sqrt(c20.*c36);
            R16 = ((S3*S41*S8-S27*S8^2-S3*S6+S29*S6*S8)*x2+S27-S29*S41)./sqrt(c20.*c25);
            R17 = ((S3*S43*S8-S28*S8^2-S3*S7+S29*S7*S8)*x2+S28-S29*S43)./sqrt(c20.*c29);
            R18 = ((S3*S45*S8-S30*S8^2-S3*S9+S29*S8*S9)*x2+S30-S29*S45)./sqrt(c20.*c32);
            R19 = ((S38*S4*S8-S31*S8^2-S4*S5+S34*S5*S8)*x2+S31-S34*S38)./sqrt(c25.*c34);
            R20 = ((S4*S41*S8-S32*S8^2-S4*S6+S34*S6*S8)*x2+S32-S34*S41)./sqrt(c25.*c36);
            R21 = ((S4*S43*S8-S33*S8^2-S4*S7+S34*S7*S8)*x2+S33-S34*S43)./sqrt(c25.*c29);
            R22 = ((S4*S45*S8-S35*S8^2-S4*S9+S34*S8*S9)*x2+S35-S34*S45)./sqrt(c25.*c32);
            R23 = ((S41*S5*S8-S36*S8^2-S5*S6+S38*S6*S8)*x2+S36-S38*S41)./sqrt(c29.*c34);
            R24 = ((S43*S5*S8-S37*S8^2-S5*S7+S38*S7*S8)*x2+S37-S38*S43)./sqrt(c29.*c36);
            R25 = ((S45*S5*S8-S39*S8^2-S5*S9+S38*S8*S9)*x2+S39-S38*S45)./sqrt(c29.*c32);
            R26 = ((S43*S6*S8-S40*S8^2-S6*S7+S41*S7*S8)*x2+S40-S41*S43)./sqrt(c32.*c34);
            R27 = ((S45*S6*S8-S42*S8^2-S6*S9+S41*S8*S9)*x2+S42-S41*S45)./sqrt(c32.*c36);
            R28 = ((S45*S7*S8-S44*S8^2-S7*S9+S43*S8*S9)*x2+S44-S43*S45)./sqrt(c34.*c36);
            z8 = S8*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S8^2*x2);
            %
            %  i=10
            %
            R1 =  ((S1*S24*S9-S10*S9^2-S1*S2+S17*S2*S9)*x2+S10-S17*S24)./sqrt(c8.*c15);
            R2 =  ((S1*S30*S9-S11*S9^2-S1*S3+S17*S3*S9)*x2+S11-S17*S30)./sqrt(c8.*c21);
            R3 =  ((S1*S35*S9-S12*S9^2-S1*S4+S17*S4*S9)*x2+S12-S17*S35)./sqrt(c8.*c26);
            R4 =  ((S1*S39*S9-S13*S9^2-S1*S5+S17*S5*S9)*x2+S13-S17*S39)./sqrt(c8.*c30);
            R5 =  ((S1*S42*S9-S14*S9^2-S1*S6+S17*S6*S9)*x2+S14-S17*S42)./sqrt(c8.*c33);
            R6 =  ((S1*S44*S9-S15*S9^2-S1*S7+S17*S7*S9)*x2+S15-S17*S44)./sqrt(c8.*c35);
            R7 =  ((S1*S45*S9-S16*S9^2-S1*S8+S17*S8*S9)*x2+S16-S17*S45)./sqrt(c8.*c36);
            R8 =  ((S2*S30*S9-S18*S9^2-S2*S3+S24*S3*S9)*x2+S18-S24*S30)./sqrt(c15.*c21);
            R9 =  ((S2*S35*S9-S19*S9^2-S2*S4+S24*S4*S9)*x2+S19-S24*S35)./sqrt(c15.*c26);
            R10 = ((S2*S39*S9-S20*S9^2-S2*S5+S24*S5*S9)*x2+S20-S24*S39)./sqrt(c15.*c30);
            R11 = ((S2*S42*S9-S21*S9^2-S2*S6+S24*S6*S9)*x2+S21-S24*S42)./sqrt(c15.*c33);
            R12 = ((S2*S44*S9-S22*S9^2-S2*S7+S24*S7*S9)*x2+S22-S24*S44)./sqrt(c15.*c35);
            R13 = ((S2*S45*S9-S23*S9^2-S2*S8+S24*S8*S9)*x2+S23-S24*S45)./sqrt(c15.*c36);
            R14 = ((S3*S35*S9-S25*S9^2-S3*S4+S30*S4*S9)*x2+S25-S30*S35)./sqrt(c21.*c26);
            R15 = ((S3*S39*S9-S26*S9^2-S3*S5+S30*S5*S9)*x2+S26-S30*S39)./sqrt(c21.*c30);
            R16 = ((S3*S42*S9-S27*S9^2-S3*S6+S30*S6*S9)*x2+S27-S30*S42)./sqrt(c21.*c33);
            R17 = ((S3*S44*S9-S28*S9^2-S3*S7+S30*S7*S9)*x2+S28-S30*S44)./sqrt(c21.*c35);
            R18 = ((S3*S45*S9-S29*S9^2-S3*S8+S30*S8*S9)*x2+S29-S30*S45)./sqrt(c21.*c36);
            R19 = ((S39*S4*S9-S31*S9^2-S4*S5+S35*S5*S9)*x2+S31-S35*S39)./sqrt(c26.*c30);
            R20 = ((S4*S42*S9-S32*S9^2-S4*S6+S35*S6*S9)*x2+S32-S35*S42)./sqrt(c26.*c33);
            R21 = ((S4*S44*S9-S33*S9^2-S4*S7+S35*S7*S9)*x2+S33-S35*S44)./sqrt(c26.*c35);
            R22 = ((S4*S45*S9-S34*S9^2-S4*S8+S35*S8*S9)*x2+S34-S35*S45)./sqrt(c26.*c36);
            R23 = ((S42*S5*S9-S36*S9^2-S5*S6+S39*S6*S9)*x2+S36-S39*S42)./sqrt(c30.*c33);
            R24 = ((S44*S5*S9-S37*S9^2-S5*S7+S39*S7*S9)*x2+S37-S39*S44)./sqrt(c30.*c35);
            R25 = ((S45*S5*S9-S38*S9^2-S5*S8+S39*S8*S9)*x2+S38-S39*S45)./sqrt(c30.*c36);
            R26 = ((S44*S6*S9-S40*S9^2-S6*S7+S42*S7*S9)*x2+S40-S42*S44)./sqrt(c33.*c35);
            R27 = ((S45*S6*S9-S41*S9^2-S6*S8+S42*S8*S9)*x2+S41-S42*S45)./sqrt(c33.*c36);
            R28 = ((S45*S7*S9-S43*S9^2-S7*S8+S44*S8*S9)*x2+S43-S44*S45)./sqrt(c35.*c36);
            z9 = S9*fi8([R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 R16 ...
                R17 R18 R19 R20 R21 R22 R23 R24 R25 R26 R27 R28])./sqrt(1-S9^2*x2);
            z(i) = wg'*(z1+z2+z3+z4+z5+z6+z7+z8+z9);
        end
    end
end

function y = setindex(n,k)
%
%   This Matlab program generates the index for different combinations
%   of submatrices of order kxk from an nxn matrix.
%
S = zeros(n);
count = 0;
for j=1:n-1
    for i=j+1:n
        count = count+1;
        S(i,j) = count;
    end
end
ind = nchoosek(1:n,k);
m = size(ind,1);
y = zeros(m,k*(k-1)/2);
for i=1:m
    ind1 = ind(i,:);
    S1 = S(ind1,ind1);
    S1 = S1(:);
    y(i,:) = S1(S1>0);
end
end
