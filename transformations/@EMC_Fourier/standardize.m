function SPECTRUM = standardize(obj, SPECTRUM)
%
% SPECTRUM = obj.standardize(SPECTRUM)
% Set the real space mean to ~0 and real space variance to ~1 of the ARRAY.
% The mean and variance are changed in frequency space.
%
% Input:
%   SPECTRUM (numeric):  2d/3d Fast Fourier Transform.
%
% Output:
%   SPECTRUM (numeric):  Standardize SPECTRUM.
%

if ~obj.half
    SPECTRUM(1) = 0;
   	SPECTRUM = SPECTRUM ./ (sqrt(sum(abs(SPECTRUM).^2, 'all')) / numel(SPECTRUM));
if obj.half
    % Capture every independant chunk (same for 2d or 3d)
    if obj.centered
        error('not finished')
        
    else
        SPECTRUM(1) = 0;
        cD = ceil(obj.size_real(1)/2);  % center of donor
        factor = sum(abs(SPECTRUM(1,:,:)).^2, 'all');  % unique row/plane
        factor = factor + 2*sum(abs(SPECTRUM(2:cD,:,:)).^2, 'all');  % common chunk
        if obj.isEven(1); factor = factor + sum(abs(SPECTRUM(cD+1,:,:)).^2, 'all'); end  % unique row/plane
    end

  	SPECTRUM = SPECTRUM ./ (sqrt(factor) / prod(obj.size_real, 'native'));
end

end  % obj.standardize
