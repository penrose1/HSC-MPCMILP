clc, clear;
np =12;
t=24;
ndays = 2;
nh = ndays*t;
c1 =1:24;
c2=25:48;
c3=49:72;
c4=73:96;
c5=97:120;
c6=121:144;

f=[c1 c2 c3 c4 c5 c6];
nvars=3;

for n =1:ndays
    for hday = 1:t/np
        for k = 1:np
            ff=[];
            for kf=1:nvars
                disp('Day: ')
                disp(n)
                disp('Var: ')
                disp(kf)
                disp('Starting index:')
                disp((kf-1)*nh+np*(n-1)+n-kf+1)
                ff=[ff f((kf-1)*nh+np*(n-1)+n-k+1:n*np+(kf-1)*nh)]
    
            end
        end
    end
    
end
% (hday-1)*np +1:np*hday