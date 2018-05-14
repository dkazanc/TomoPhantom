function [noisy_sino,N_noise] = add_noise(b, sigma, noisetype)
% function which adds Gaussian or Poisson noise to data
% return data after negative log and also raw data

if (strcmp(noisetype, 'Gaussian') == 1)
    E = randn(size(b));
    noisy_sino = b + sigma*norm(b,'fro')*E/norm(E,'fro');  % adding normal noise to the sinogram
    noisy_sino(noisy_sino<0) = 0;
elseif (strcmp(noisetype, 'Poisson') == 1)
    maxVal = max(b(:));
    b = b./maxVal;
    dataExp = sigma.*exp(-b); % noiseless raw data
    N_noise = poissrnd(dataExp);  % adding Poisson noise to the sinogram
    noisy_sino = -log(N_noise./max(N_noise(:)))*maxVal; %log corrected data -> sinogram
    noisy_sino(noisy_sino<0) = 0;
else
    error('Please select Gaussian or Poisson noise');
end

return;
