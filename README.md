Steps Involved to model a 4x4 UPA and implement DOA estimation cum Visualization using MATLAB
4x4 Uniform Planar Array (UPA): The array is configured as a 4x4 grid, resulting in 16 antenna elements. The element spacing d is half the wavelength (lambda/2).
Steering Vector Calculation: The steering vector now accounts for two dimensions:
m and n iterate over the rows and columns of the UPA, respectively.
The phase shift is calculated for each element based on its position (m, n) in the planar array.
MUSIC Spectrum: The MUSIC pseudo-spectrum is calculated over both azimuth (theta) and elevation (phi) angles, resulting in a 2D spectrum. Peaks in this spectrum indicate the estimated directions of arrival.
Visualization: The MUSIC spectrum is plotted using surf to show the estimated DOA in both azimuth and elevation angles.


# DOAEstimation

% DOA estimation using for far fields using matlab dataset generation

DOA=[-2,3];
T = 500;
K= length(DOA); 
Nr =4; 
Nc=4;
c=3e8;
fc=3.3e9;
lambda=c/fc;
d = lambda/2;
SNR = 1;    
A = zeros(Nr*Nc,K); 
 for k=1:K 
     steering_vector = zeros(Nr, Nc);
    for m=1:Nr
        for n=1:Nc
            steering_vector(m,n)=exp(-1j*2*pi*d*(m*sind(DOA(k))+n*cosd(DOA(k)))/lambda);
        end

    end
    A(:,k)=reshape(steering_vector,[],1);

end 
Vj = diag(sqrt((10.^(SNR/10))/2));
s = Vj * (randn(K, T) + 1j * randn(K, T));
noise = sqrt(1/2) * (randn(Nr * Nc, T) + 1j * randn(Nr * Nc, T));
X = A * s;
X = X + noise;       
Rx = cov(X');                    
[eigenVec,eigenVal] = eig(Rx);    
Vn = eigenVec(:,1:Nr*Nc-K);         
theta = -90:0.05:90;
Pi=-90:0.05:90;
for i=1:length(theta) 
    for j=1:length(Pi)
        SS=zeros(Nr,Nc);
        for m=1:Nr
            for n=1:Nc
                SS(m,n)=exp(-1j*2*pi*d*(m*sind(theta(i))+n*cosd(Pi(j)))/lambda);
            end
        end
        SS=reshape(SS,[],1);
        PP=SS'*(Vn*Vn')*SS;
        Pmusic(i,j)=1/PP;

    end
end
Pmusic = real(10*log10(Pmusic));
% Find peaks in the spectrum
[pks, locs] = findpeaks(Pmusic(:), 'SortStr', 'descend', 'Annotate', 'extents');
[est_i, est_j] = ind2sub(size(Pmusic), locs(1:K));
estimated_angles = [theta(est_i)', Pi(est_j)'];

% Display estimated angles
disp('Estimated DOA angles (Azimuth, Elevation):');
disp(estimated_angles);

% Plot MUSIC Spectrum
figure;
surf(Pi, theta, Pmusic, 'EdgeColor', 'none');
xlabel('Elevation Angle \phi (degree)');
ylabel('Azimuth Angle \theta (degree)');
zlabel('Spatial Power Spectrum P(\theta, \phi) (dB)');
title('DOA estimation based on MUSIC algorithm for 4x4 UPA');
axis tight;
view(2);
colorbar;
grid on;
