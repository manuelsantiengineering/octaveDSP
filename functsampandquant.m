function [ s_out, sq_out, sqh_out, Delta, SQNR ] = functsampandquant(sig_in, L, td, ts )
  % The function executes sampling and uniform quantization simultaneously.
  % Input:
  % L - number of uniform quantization levels
  % sig_in - input signal vector
  % td - original signal sampling period of sig_in
  % ts - new sampling Period
  % Output:
  % s_out - sampled Output
  % sq_out - sampled and quantized Output
  % sqh_out - sampled, quantized, and zero order hold Output
  % Delta - quantization interval
  % SQNR - actual signal to quantization noise ratio
  if(rem(ts/td,1) == 0) %rem returns the reminder of the division
    nfac = round(ts/td);
    p_zoh = ones(1, nfac);
    s_out = downsample(sig_in, nfac);
    [sq_out, Delta, SQNR] = functuniquan(s_out, L);
    s_out = upsample(s_out, nfac);
    sqh_out = kron(sq_out, p_zoh);
    sq_out = upsample(sq_out, nfac);
  else
    s_out = []; sq_out = []; sqh_out=[]; Delta=[]; SQNR=[];
  end
end  % function
