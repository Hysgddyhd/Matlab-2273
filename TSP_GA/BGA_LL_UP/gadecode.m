function f=gadecode(chrom,lo,hi,bits)
[M,N]=size(chrom);
npar=N/bits; % number of variables
quant=(0.5.^[1:bits]'); % quantization levels
quant=quant/sum(quant); % quantization levels normalized
ct=reshape(chrom',bits,npar*M)';% each column contains one variable
par=((ct*quant)*(hi-lo)+lo); % DA conversion and unnormalize variables
f=reshape(par,npar,M)'; % reassemble population