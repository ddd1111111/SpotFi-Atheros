% sample CSI trace is a 90x1 vector where first 30 elements correspond to subcarriers for first rx antenna, second 30 correspond to CSI from next 30 subcarriers and so on.
% replace sample_csi_trace with CSI from Intel 5300 converted to 90x1 vector
calculate_number=0;
MUSIC_iter = 1;
format short g;
sanitization=0;
sample_csi_traceTmp_r(1,168)=0;
particle_filter=1;
d_noise= 3; % standard deviations of the dynamic noise.
measurement_noise = 0;% means for the measurement noise.
m_noise= 3;% standard deviations for the measurement noise.
init=1;
sub_carrier_number=1;
