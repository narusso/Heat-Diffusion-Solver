fid=fopen('soln.dat','r');   % open data file
fseek(fid,0,'eof');          % go to end of file
position = ftell(fid);       % determine size of file
frewind(fid);                % go back to the beginning
n = sqrt(position/8/10);     % assume file contains 2x2 matrix
for i=1:40
  i
  f=fread(fid,[n,n],'double'); % of 8-byte doubles
  surf(f);                     % draw without contour lines
  %zlim([0 1e-3]);                 % fix the z-axis
  shading flat;                % disable black edges
  colormap(hsv(128));          % use more and different colors
  pause(1);
end
fclose(fid);                 % close data file
