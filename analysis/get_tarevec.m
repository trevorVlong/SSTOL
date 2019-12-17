function [tarevec] = get_tarevec(filename)
%GET_TAREVEC Returns the tare vector from a given tare file
%   INPUTS
%
%   filename    -   String giving the name of the tare file. Tare file 
%                   must be in .xlsx format. .xlsx should not be included in 
%                   the file name e.g. filename = 'tare_0001' implies the
%                   opening of the file tare_0001.xlsx
%
%
%
%   OUTPUTS
%
%   tarevec     - 3x1 matrix of wind-off tare values
%
%   Aug 26. 2019
%   courtin@mit.edu

opts = spreadsheetImportOptions("NumVariables", 12);
opts.Sheet = 1;

data = readtable(strcat(filename,'.xlsx'), opts, "UseExcel", false);

[nr,nc] = size(data);

L_vals = table2array(data(2:end,2));
X_vals = table2array(data(2:end,3));
M_vals = table2array(data(2:end,4));

L_sum = 0;
X_sum = 0;
M_sum = 0;

for i = 1:nr-1
    %Kinda gross; can be more elegant with better understanding of data
    %type coversions
    L_sum = L_sum + str2num(L_vals{i});
    X_sum = X_sum + str2num(X_vals{i});
    M_sum = M_sum + str2num(M_vals{i});
end

L = L_sum/(nr-1);
X = X_sum/(nr-1);
M = M_sum/(nr-1);

tarevec = [L;X;M];



end

