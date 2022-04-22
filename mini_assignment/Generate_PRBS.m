function [r,P,w] = Generate_PRBS(probability, N, mag, plt_flag)
%GENERATE_PRBS Summary of this function goes here
%  generate PRBS based on:
%  u(t) = u(t-1) with probability p
%  u(t) = -u(t-1) with probability 1-p
%  and may plot spectral density based on setting
%INPUT:
%  P : prbability
%  N : total input points
%  mag: output magnitude
%  plt_flg: if true, plot sepctral density
%OUTPUT:
%  u : generated input sequence
%  P : spectral density of P
%  w : spectral density points 

    % determine the initial entry of u, that is u(1)
    if rand(1,1)>=0.5
        r(1) = mag;
    else
        r(1) = -mag;
    end
    
    % generate sequence
    for i = 2:1:N
        if rand(1,1) <= probability
            r(i) = r(i-1);
        else
            r(i) = -r(i-1);
        end
    end
    
    % get spectrum of generated input sequence
    [P,w] = cpsd(r,r,N);
    
    if (plt_flag == true)
        cpsd(r, r);
        title('spectrum of generated input sequence')
    end
    
end

