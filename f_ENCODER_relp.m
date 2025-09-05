%RELP ENCODER portion:

function [aCoeff, b_LTopt, Topt, e_prime] = f_ENCODER_relp(x, fs);


M = 8; %prediction order for LP analysis

%INITIALIZATION;
b=1;        %index no. of starting data point of current frame
fsize = 20e-3;    %frame size (in milisec)   [chap10, page6, Pract Handbook of Speech Coders]
frame_length = round(fs .* fsize);   %=number data points in each framesize 
                                %of "x"
N= frame_length - 1;        %N+1 = frame length = number of data points in 
                            %each framesize

y_proc = filter([1 -1], [1 -0.999], x);  %pre-processing
%FRAME SEGMENTATION
for b=1 : frame_length : (length(x) - N),
% for b=1 : frame_length : frame_length.*8, %(length(x) - frame_length),  %temporary
% y_f=y_proc(1280:1280+159);   %temporary

y_f = y_proc(b:b+N);     %"b+N" denotes the end point of current frame.
                %"y" denotes an array of the data points of the current 
                %frame
        
    %LP ANALYSIS [lev-durb] & PREDICTION ERROR (short-term) FILTER;
    [a, tcount_of_aCoeff, e_s] = func_lev_durb_relp (y_f, M); %e=error signal from lev-durb proc
    aCoeff(b: (b + tcount_of_aCoeff - 1)) = a;  %aCoeff is array of "a" for whole "x"
    
    %LONG-TERM LP ANALYSIS, FILTERING, AND CODING
    %analysis:
    T_min = round (fs .* 5e-3);  %=total data points in 5ms of "x"      +++++++++++++++++++++++++++++ (pg291)
                                            %[chap10, page6, Pract Handbook of Speech Coders]
    T_max = round (fs .* 15e-3);   %=15ms+++++++++++++++++++  (pg291)
%   b = 15841;   %temporary
    c1 = 1;
    for bs = b : 40 : b+length(y_f)-40,  %subframing
%       bs = 1281;   %temp
        if bs < T_max,
            break;
        end
    
        Jmin(bs) = 10^9;
    
        for T = T_min : T_max,      %within 1 (current) frame
%               T = 40;   %temp
            
                %b loop:
                for c = 1:40,    %data points of current subframe
%                     c=1;  %temporary
                    sm1(c) = ( y_proc(bs+(c-1)) .* y_proc(bs-T+(c-1)) );  %"y_proc(bs+(c-1))" =es(n) of pg 121
                    sm2(c) = y_proc(bs-T+(c-1)); %=es(n-T) of pg 121
                    sm22(c) = sm2(c).^2;
                end
                q1 = sum(sm1);
                q2 = sum(sm22);
                b_LT(T) = -(q1./q2);
            
                %J loop:
                for c = 1:40,    %data points of current subframe
%                     c=1;  %temporary
                    smJ1(c) = y_proc(bs+(c-1));
                    smJ2(c) = b_LT(T) .* y_proc(bs-T+(c-1));            
                end            
                smJ = smJ1 + smJ2;
                qJ = smJ.^2;
                J(T) = sum(qJ);
            
                if J(T) < Jmin(bs),
                    Jmin(bs) = J(T);
                    Topt(bs) = T;
                    if b_LT(T)>=1,
                        b_LTopt(bs) = 0.9999; %trancation
                    else
                        b_LTopt(bs) = b_LT(T);
                    end
                end
        end     %T loop ends
        %predictor:
        LT_gain = [zeros(1, Topt(bs)-1), b_LTopt(bs)];     %as it says z^-T in page 121
         e_s_padded_for_delay_removal = [e_s(c1:c1+39); zeros(Topt(bs), 1)]; %it is padded with zeros to remove the effect of delay in filter. %Topt(bs) no. of 'z's and one '1' results in total 'Topt(bs)' amount of delay
        e_with_dummy_pad = filter([1 LT_gain], 1,  e_s_padded_for_delay_removal);   % = 1 + 0*z^-1 + 0*z^-2 + ... + b*z^-T
         e_LT(bs:bs+39,1) = e_with_dummy_pad(Topt(bs)+1 : Topt(bs)+1+39);   %LT predicted "e"
        e(bs:bs+39, 1) = e_s(c1 : c1+39) - e_LT(bs : bs+39);
        
        %WEIGHTING FILTER:
        w = flattopwin(11);     %11 point flattop window is temporarily chosen for page 292
                            %use wvtool(w) in command window to see its
                            %response
        wndd = conv(w, e(bs:bs+39));     %outputs total 50 samples
        x_n(bs:bs+39) = wndd(6:45);  %middle 40 samples are taken

        %POSITION SELECTION & EXCITATION GENERATOR:
        for i1 = 0:3,  % = 0:3 of page 293
            for i = i1+bs : 3 : bs+i1+38;
                x_m(i1+1,i) = x_n(i);
            end
            E_m(i1+1,1) = sum((x_m(i1+1, bs:end)).^2);
        end
        [E_m_max, index_max] = sort(E_m);
        e_prime(bs : bs+39) = x_m(index_max(4), bs:bs+39); %=e_prime of fig 10.5

        c1 = c1 + 40;
    end
end

    
    
   