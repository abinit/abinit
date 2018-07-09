function plot_band
%Anton Bokhanchuk January 6, 2015
%Plot band structure
data=load('H6o_DS2_EIG.dat');
DatatoPlot = [data(2:end,1) data(2:end,6:end)]; %Create plot matrix
hold on
for k = 1:size(DatatoPlot,1)
plot(k,DatatoPlot(k,2:end),'*-') %plot data
end
xlabel('K-point path (Y-G-Y; Z-G-Z)');
ylabel('Energy (eV)');
xlim([0 k+1]) %// To see data more clearly
set(gca,'Xtick');
end

