% ECE490 Lab 3 
%% Clears all the workspace
clc
clear
%%

load('molecular_absorption.mat') % load the file which contains frequency absorption data

% The file above loads 'f' and 'k'


%% Section 1 - Finding the frequency of operation 
dummy_variable = 1;
for i = 1:length(f)
    if (f >= 100e8) && (f <= 1.5e12)
       freq_2_use(dummy_variable) = f;
       k_2_use(dummy_variable) = k;
       dummy_variable = dummy_variable + 1;
    end
end 


%% Finding the distances
r = [1  10  100]; % distances
%% Section 2 - Absorption Loss 
for i = 1:length(r)
    for j = 1:length(freq_2_use)
    ab_loss(i,j) = 1./(exp(-k_2_use(j).*r(i)));
    ab_loss_db (i,j) = 10*log(ab_loss(i,j));
    end
end


%% Section 3 - plotting the absorption loss
figure(1)
box on; 
h1 = plot(freq_2_use, ab_loss_db(1,:)); % plots the absorption loss for the first link 
hold on; 
h2 = ;% plot the absorption loss in dB for the second link distance 
h3 = ;% plot the absorption loss in dB for the third link distance 
legend([h1 h2 h3], {'Legend 1 ', 'Legend 2 ', 'Legend 3 '},'Location','northeast') % give your own names which are self explanatory
ylim ([0 300]) % a limit on the range of the absorption loss 
xlim([% the desired range])  % another limit on the range of the frequencies 
xlabel('Frequency (Hz)') 
ylabel('Absorption Loss (dB)')
%% Section 4 - Spreading Loss
c_const = 3e8; % speed of light 
wavelength = ; % find the wavelength from the frequencies.. 
A_eff = wavelength.^2./(4.*pi); % effective aperture of a dipole
 for i = 1:length(r)
     for j = 1:length(A_eff)
         spr_loss(i,j) = 4.*pi.*r(i).^2./(A_eff(j));
         spr_loss_db(i,j) = ; % find the spreading loss in dB
         path_loss_db(i,j) = ; % total path loss in dB
     end
 end

 %% Section 5 - Plotting the Spreading and Path Loss
 figure(2)
    h4 = ; % spreading loss in dB for first link distance 
    hold on 
    h5 = ; % spreading loss in dB for second link distance 
    h6 = ; % spreading loss in dB for third link distance  
    legend([h4 h5 h6], {' Legend titles'},'Location','northeast')
title('Spreading loss at different distances')
xlabel('Frequency (Hz)')
ylabel('Spreading loss(db)')
hold off;

 figure(3)
    h7 = ; %path loss in dB for first link distance 
    hold on 
    h8 = ; %path loss in dB for second link distance 
    h9 = ; %path loss in dB for third link distance  
    hold off;
legend([h7 h8 h9], {'legend titles'},'Location','northeast')
title('Path loss at different distances')
xlabel('Frequency (Hz)')
ylabel('Path loss(db)')
 
%% Section 6 - Comparison between the measured and simulated path losses
    %Q2
    load('channel.mat')
    %%
   %  dummy_variable = 1;
    for i = 1:length(freq_2_use)
    if % new range to plot data in 
       
        % your code here
    end
    end
% add your code for the new range of wavelength and affective aperture

r_new = [16e-2]; % new link distance 

    for j = 1:length(new_freq_range_2_use)
        % absorption loss
        % absorption loss in dB 
        % spreading loss 
        % spreading loss in dB
        % path loss in dB
    end

%% Section 7 - Plot and compare the two 
 figure(4)
    h10 = plot(new_freq_range_2_use, new_path_loss_in_dB);
    hold on 
    h11 = plot(f,H, '--*'); % the measured values 
    ylim([65 85])
    legend([h10 h11], {'Simulated data title', 'Measured path loss'},'Location','northeast')
    hold off;
title('A comparison between measured and simulated losses')
xlabel('Frequency (Hz)')
ylabel('Path loss(db)')

   
