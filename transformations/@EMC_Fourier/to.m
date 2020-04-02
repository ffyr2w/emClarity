function obj = to(obj, STR)
%
% obj = obj.to(STR)
% Change the precision or the method of the EMC_Fourier object.
%
% Input:
%   STR (str|char): If 'single' or 'double': cast obj properties to the desired precision.
%                   If 'cpu' or 'gpu': push object properties to the desired device/host.
%
% Output:
%   obj (handle):   EMC_Fourier instance.
%


switch STR
    case 'cpu'
        if ~strcmpi(obj.method, 'cpu')
         	obj.method = 'cpu';
            obj.isOnGpu = false;
            if ~isemtpy(obj.bandpass); obj.bandpass        = gpuArray(obj.bandpass);        end
            if index_fftshift_isSet;   obj.index_fftshift  = gpuArray(obj.index_fftshift);  end
            if index_ifftshift_isSet;  obj.index_ifftshift = gpuArray(obj.index_ifftshift); end
            if index_half2full_isSet;  obj.index_half2full = gpuArray(obj.index_half2full); end
        end

    case 'gpu'
        if ~strcmpi(obj.method, 'gpu')
            obj.method = 'gpu';
            obj.isOnGpu = true;
            if ~isemtpy(obj.bandpass); obj.bandpass        = gather(obj.bandpass);        end
            if index_fftshift_isSet;   obj.index_fftshift  = gather(obj.index_fftshift);  end
            if index_ifftshift_isSet;  obj.index_ifftshift = gather(obj.index_ifftshift); end
            if index_half2full_isSet;  obj.index_half2full = gather(obj.index_half2full); end
        end

    case 'double'
        if ~strcmpi(obj.precision, 'double')
            obj.precision = 'double';
            if ~isemtpy(obj.bandpass); obj.bandpass = cast(obj.bandpass, 'double'); end
        end

    case 'single'
        if ~strcmpi(obj.precision, 'single')
            obj.precision = 'single';
            if ~isemtpy(obj.bandpass); obj.bandpass = cast(obj.bandpass, 'single'); end
        end

    otherwise
        error('EMC:Fourier:to', "STR should be 'gpu', 'cpu', 'double' or 'single'")
end

end  % to
