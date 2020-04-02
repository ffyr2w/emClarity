function IMAGE = ifft(obj, SPECTRUM)
%
% IMAGE = obj.ifft(SPECTRUM)
% N-D inverse Fast Fourier Transform of a complex 2d/3d non-redundant/redundant SPECTRUM.
%
% Input:
%   SPECTRUM (numeric): 2d or 3d, non-redundant or redundant, Fast Fourier Transform (complex).
%
% Output:
%   IMAGE (numeric):    Inverse Discrete Fourier Transform of SPECTRUM (real).
%
% Property used:
%   obj.half
%   obj.half_wrap
%   obj.size_real
%   obj.fftshift
%   obj.index_ifftshift
%   obj.index_half2full
%
% Method used:
%   EMC_Fourier.half2full
%
% Note:
%   - Recomputing the redundant SPECTRUM has a non negligeable overload.
%     This is only worth it if a lot of operation are done on the non-redundant.
%
% Example:
%   - >> ft = EMC_Fourier([128,128,128], 'gpu', {});
%     >> img = rand(128,128,128);
%     >> dft = ft.fft(img);
%     >> img_back = ft.ifft(dft);  % img_back == img
%

% Created:  3Feb2020
% Version:  v.1.0.  unittest (TF, 8Feb2020).
%           v.1.1.  use EMC_maskIndex and a persistent wrapper to speed up
%                   Fourier coefficients wrapping (half to full grid) (TF, 14Feb2020).
%           v.1.2.  EMC_irfftn is integrated into EMC_Fourier.ifft; can used
%                   half/full not-centered/centered DFTs (TF, 8Mar2020).
%

if obj.half
    IMAGE  = ifftn(EMC_Fourier.half2full(obj.wrap, SPECTRUM, obj.size_real, obj.index_half2full), 'symmetric');
elseif obj.centered
    IMAGE  = ifftn(SPECTRUM(obj.index_ifftshift), 'symmetric');
else
    IMAGE  = ifftn(SPECTRUM, 'symmetric');
end

end  % ifft
