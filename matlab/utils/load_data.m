%==========================================================================
% function [u_l,u_r,q_l,q_r] = load_data(path)
%==========================================================================
% @descirption: Extracts data from the APDM sensors and stores them in a
%               matrix. Also extracts processed quaternion information from 
%               the APMD sensors.
% @author     : Prateek Gundannavar
% @date       : 06/07/18
%
% @input      
%             - path        location of the APDM .h5 files
% @output
%             - u_l         9xN matrix with acc, gyro, and mag measurements
%             - u_r         9xN matrix with acc, gyro, and mag measurements
%             - q_l         4xN matrix with processed quaternions
%             - q_r         4xN matrix with processed quaternions
%
%               (subscripts _l and _r indicate *left* & *right* foot, resp)
% @copyright   : Copyright(c) 2019 Prateek Gundannavar
%==========================================================================
function [u_l,u_r,q_l,q_r,t_l,t_r] = load_data(path)

l_foot = '/Sensors/1516';
u_al = h5read(path,strcat(l_foot,'/Accelerometer'));
u_gl = h5read(path,strcat(l_foot,'/Gyroscope'));
u_ml = h5read(path,strcat(l_foot,'/Magnetometer'));
u_bl = h5read(path,strcat(l_foot,'/Barometer'));
u_l = [u_al;u_gl;u_ml;u_bl'];
t_l = h5read(path,strcat(l_foot,'/Time'));
t_l = convert_unix_to_real_time(t_l);
l_foot = '/Processed/1516';
q_l = h5read(path,strcat(l_foot,'/Orientation'));

r_foot = '/Sensors/1989';
u_ar = h5read(path,strcat(r_foot,'/Accelerometer'));
u_gr = h5read(path,strcat(r_foot,'/Gyroscope'));
u_mr = h5read(path,strcat(r_foot,'/Magnetometer'));
u_br = h5read(path,strcat(r_foot,'/Barometer'));
u_r = [u_ar;u_gr;u_mr;u_br'];
t_r = h5read(path,strcat(r_foot,'/Time'));
t_r = convert_unix_to_real_time(t_r);
r_foot = '/Processed/1989';
q_r = h5read(path,strcat(r_foot,'/Orientation'));


end

function real_time = convert_unix_to_real_time(unix_time)

unix_time = double(unix_time);
time_reference = datenum('1970', 'yyyy'); 
time_matlab = time_reference + unix_time / 8.64e10;         % micro seconds
real_time = datestr(time_matlab, 'yyyy-mm-dd HH:MM:SS.FFF');
real_time = datetime(real_time,'Format','yyyy-MM-dd HH:mm:ss.SSS','TimeZone','America/Chicago') + hours(-5) + seconds(-11.35);

end

