
function [ q_out, Delta, SQNR] = functuniquan(sig_in, L)
  % The function executes sampling and uniform quantization simultaneously.
  % Input:
  % L - number of uniform quantization levels
  % sig_in - input signal vector
  % Output:
  % q_out - quantized Output
  % Delta - quantization interval
  % SQNR - actual signal to quantization noise ratio

sig_pmax=max(sig_in); %Positive peak
sig_nmax=min(sig_in); %Negative peak
Delta =(sig_pmax-sig_nmax)/L;
q_level=sig_nmax+Delta/2:Delta:sig_pmax-Delta/2; %Defines the quantization levels
% L_sig=length(sig_in);
sigp=(sig_in-sig_nmax)/Delta+1/2; %Converts into 1/2 to L+1/2 range
qindex=round(sigp); %Round to 1,2,... L levels
qindex=min(qindex,L); %eliminate L+1 as rate possibility
q_out=q_level(qindex); %Use index vector to generate output
SQNR = 20*log10(norm(sig_in)/norm(sig_in-q_out));
end
