
%% Junk I just testet smth
a = 1;
del = 1;
L = 5;
n = 2;

alpha(10, a, del, n, L)


function T = T1(k, a, del, n)

    A = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];
    B = [1 + n/del, 1 - n/del; 1 - n/del, 1 + n/del];
    C = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];

    T =  0.5 * A * B * C;
end

function T = T2(k, a, del, n)

    A = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];
    B = [1 + del/n, 1 - del/n; 1 - del/n, 1 + del/n];
    C = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];

    T =  0.5 * A * B * C;
end

function F = FL(k, L)

    F = [exp(-1i * k * L), 0 ; 0, exp(1i * k * L)];
end 

function M = M(k, a, del, n, L)

    M = T1(k, a, del, n) * T2(k, a, del, n) * FL(k,L);
end

function b = beta(k, a, del, n, L)

    eigenvalues = eig(M(k, a, del, n, L));
    %b = 1 / L * log(abs(eigenvalues(1)));
    b = [1 / L * log(abs(eigenvalues(1))), 1 / L * log(abs(eigenvalues(2)))];
end 

function a = alpha(k, a, del, n, L)

    eigenvalues = eig(M(k, a, del, n, L));
    %a = angle(eigenvalues(1));
    a = [angle(eigenvalues(1)), angle(eigenvalues(2))];
end
