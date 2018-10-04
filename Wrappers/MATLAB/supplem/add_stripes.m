function [stripes_sino] = add_stripes(b, percentage, maxthickness)
%  function to add stripes (constant offsets) to sinogram which results in rings in the reconstructed image
%        - percentage defines the density of stripes
%        - maxthickness defines maximal thickness of a stripe

if ((percentage <= 0) || (percentage > 100))            
error('Percentage must be larger than zero but smaller than 100');
end
if ((maxthickness < 0) || (maxthickness > 5))    
error('Maximum thickness must be in [0,5] range');
end

DetectorsDim = size(b,2);
range_detect = round((DetectorsDim)*(percentage/100));
stripes_sino = b;
max_intensity = max(b(:));

for jj = 1:range_detect
    randind = randi([1 DetectorsDim],1); % generate random index of a faulty detector
    randthickness = randi([0 maxthickness],1); % generate random thickness    
    randintens = -1 + (0.5 + 1).*rand(1,1); % generate random multiplier
    intensity = max_intensity*randintens;       
    for ll = -randthickness:randthickness
        if (((randind+ll) > 0) && ((randind+ll) < DetectorsDim))   
            stripes_sino(:,randind+ll) = stripes_sino(:,randind+ll) + intensity;
        end
    end
end
stripes_sino(stripes_sino < 0) = 0;
return;
