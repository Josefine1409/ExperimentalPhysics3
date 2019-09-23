A=importdata('Compton_110deg_15min_ch001.dat')
fid = fopen('Compton_110deg_15min_ch001.dat')
b = fread(fid);
A.data;
 
 
data = fread(fid, '*uint8');
temp = fread(fid, [2 inf], 'double');
