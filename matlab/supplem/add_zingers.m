function [zingers_sino] = add_zingers(b, percentage, modulus)
% adding zingers (zero pixels or small 4 pixels clusters) to data 
%        - percentage - the amount of zingers to be added
%        - modulus controls the amount of 4 pixel clusters to be added

if ((percentage <= 0) || (percentage > 100))            
error('Percentage must be larger than zero but smaller than 100');
end
if (modulus < 0)
error('Modulus integer must be positive');
end

length_sino = length(b(:));
num_values = round((length_sino)*(percentage/100));
zingers_sino = b;
     
for jj = 1:num_values
    randind = randi([1 length_sino],1); % generate random indeces 
    zingers_sino(randind) = 0;            
    if (mod(jj,modulus) == 0)
        %make a cluster of 4 pixels
        if ((randind > size(b,1)) && (randind < length_sino-size(b,1)))
           zingers_sino(randind+1) = 0;          
           zingers_sino(randind-1) = 0;          
           zingers_sino(randind+size(b,1)) = 0;          
           zingers_sino(randind-size(b,1)) = 0;          
        end
    end
end

return;
