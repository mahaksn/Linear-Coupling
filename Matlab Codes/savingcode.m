function [filename,folder] = savingcode(RD,dim,suffix,layer)
model = strcat(RD,'_',dim,'_',layer,'_',suffix);
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
filename = strcat(model,'_',time);

folder = ['D:\MATLAB Output\Linear Coupling\',model,'\',];
end
