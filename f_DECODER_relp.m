%RELP DECODER portion:

function [synth_speech, synth_speech1, LT_gain, e_prime_padded_for_delay_removal, e_prime_op_with_dummy_pad,e_prime_op, e_prime_op_padded_for_delay_removal,synth_speech_with_dummy_pad] = f_DECODER_relp(aCoeff, b_LTopt, Topt, e_prime);

%re-calculating frame_length for this decoder,
frame_length=9; %initial value for calculation
for i=10:length(aCoeff)
    if aCoeff(i) == 0,
        frame_length = frame_length + 1;
    else break;
    end
end     %=160?

e_prime = e_prime';     %making it a column matrix for convenience

for b=1 : frame_length : length(aCoeff),    %length(aCoeff) should be very close 
                                            %(i.e less than a frame_length error) to length(x)
    for bs = b : 40 : b+frame_length-40,  %subframing
        
        %EXCITATION GENERATOR:
        %not done yet. because e_prime has been sent to this decoder directly. without quantization.
        
        %PITCH SYNTHESIS FILTER:    %has to be done per subframe
        LT_gain = [zeros(1, Topt(bs)-1), b_LTopt(bs)];     %as it says z^-T in page 121
         e_prime_padded_for_delay_removal = [e_prime(bs:bs+39); zeros(Topt(bs), 1)]; %it is padded with zeros to remove the effect of delay in filter. %Topt(bs) no. of 'z's and one '1' results in total 'Topt(bs)' amount of delay
        e_prime_op_with_dummy_pad = filter(1, [1 LT_gain],  e_prime_padded_for_delay_removal);  % = 1 / (1 + 0*z^-1 + 0*z^-2 + ... + b*z^-T)
         e_prime_op(bs:bs+39,1) = e_prime_op_with_dummy_pad(Topt(bs)+1 : Topt(bs)+1+39);   %pitch-synthesis filter output
    end

        %FORMANT SYNTHESIS FILTER:
         e_prime_op_padded_for_delay_removal = [e_prime_op(b : b+159); zeros(1,1)]; %it is padded with zeros to remove the effect of delay in filter
        synth_speech_with_dummy_pad = filter(1, [1 aCoeff(b+1 : b+8)], e_prime_op_padded_for_delay_removal);    
         synth_speech1(b : b+159) = synth_speech_with_dummy_pad (2:end);  % = s^(n) with a cap on page 129 of the book          
        
        %DE-EMPHASIS (de-proprocessing):
        synth_speech(b : b+159) = filter([1 -0.999], [1 -1], synth_speech1(b : b+159));  %De-processing

end


