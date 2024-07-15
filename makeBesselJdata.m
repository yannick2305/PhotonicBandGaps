function Jdata = makeBesselJdata2(N_multi,z)


Jdata=zeros(1,2*N_multi+1);

%Jdata(N_multi+1+0)=besselj(0,z);
%Jdata(N_multi+1+1)=besselj(1,z);

for n=-N_multi:N_multi
   
    Jdata(N_multi+1+n) = besselj(n,z);
    
end


%for n=2:N_multi
   
%    Jdata(N_multi+1+n)=-Jdata(N_multi+1+n-2)+2*(n-1)/z*Jdata(N_multi+1+n-1);
    
%end

%for n=-1:-1:-N_multi
   
%    Jdata(N_multi+1+n)=-Jdata(N_multi+1+n+2)+2*(n+1)/z*Jdata(N_multi+1+n+1);
    
%end

% -besselj(n-2,x)+2*(n-1)/x*besselj(n-1,x)
% besselj(n,x)
% 2*(n+1)/z*besselj(n+1,z)+besselj(n+2,z)

end