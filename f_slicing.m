clc, clear ,close all;
x=1:47*21;
for kk= 1:24
    disp(kk)
    for k =  1:21
        disp(x((k-1)*47+kk:24+kk-1+(k-1)*47)')
    end
end