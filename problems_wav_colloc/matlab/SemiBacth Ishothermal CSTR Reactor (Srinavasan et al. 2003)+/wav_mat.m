function [DWT iDWT] = wav_mat(PH, Hi ,J)

% Integral do produto de Hs
    dim=     2^(J+1);
    DWT =    inv(Hi)*inv(eye(dim,dim))*inv(PH)*inv(eye(dim,dim));
    iDWT  =  PH*Hi;
    
end 