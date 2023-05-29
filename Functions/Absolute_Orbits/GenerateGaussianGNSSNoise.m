function [r_noise, v_noise] = GenerateGaussianGNSSNoise(r_noise_std, v_noise_std, n)
    r_sigma = eye(3) * r_noise_std; 
    v_sigma = eye(3) * v_noise_std; 

    r_noise = sqrtm(r_sigma) * randn(3,n);
    v_noise = sqrt(v_sigma) * randn(3,n);
end
