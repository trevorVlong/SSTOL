function testlog = implog(workbookFile, sheetName, dataLines)
%IMPORTFILE1 Import data from a spreadsheet
%  TESTLOG = IMPORTFILE1(FILE) reads data from the first worksheet in
%  the Microsoft Excel spreadsheet file named FILE.  Returns the data as
%  a table.
%
%  TESTLOG = IMPORTFILE1(FILE, SHEET) reads from the specified worksheet.
%
%  TESTLOG = IMPORTFILE1(FILE, SHEET, DATALINES) reads from the
%  specified worksheet for the specified row interval(s). Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  testlog = importfile1("/Users/trevorlong/Dropbox (MIT)/Wind Tunnel Testing/Testing_Summer2019/testlog.xlsx", "Sheet1", [2, 108]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 03-Sep-2019 13:59:15

%% Input handling
workbookFile = 'testlog.xlsx';
% If no sheet is specified, read first sheet
if nargin == 0 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 200];
end

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":K" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["filename", "DateCollected", "FlapAngledegrees", "AoAlowdegrees", "AoAhigh", "DynamicPressurePa", "AmbientPressurePa", "TunnelTempK", "CalibrationFile", "TareFile", "GoodPoint"];
opts.SelectedVariableNames = ["filename", "DateCollected", "FlapAngledegrees", "AoAlowdegrees", "AoAhigh", "DynamicPressurePa", "AmbientPressurePa", "TunnelTempK", "CalibrationFile", "TareFile", "GoodPoint"];
opts.VariableTypes = ["string", "datetime", "double", "double", "double", "double", "double", "double", "string", "string", "categorical"];
opts = setvaropts(opts, 2, "InputFormat", "");
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 9, 10, 11], "EmptyFieldRule", "auto");

% Import the data
testlog = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":K" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    testlog = [testlog; tb]; %#ok<AGROW>
end

    testlog.Properties.VariableNames{'FlapAngledegrees'} = 'FA';
    testlog.Properties.VariableNames{'AoAlowdegrees'}    = 'AoAlo';
    testlog.Properties.VariableNames{'AoAhigh'}          = 'AoAhi';
    testlog.Properties.VariableNames{'DynamicPressurePa'} = 'DP';
    testlog.Properties.VariableNames{'AmbientPressurePa'} = 'AP';
    testlog.Properties.VariableNames{'TunnelTempK'} = 'TT';
    testlog.Properties.VariableNames{'CalibrationFile'} = 'Cfile';
    testlog.Properties.VariableNames{'TareFile'} = 'Tfile';
    testlog.Properties.VariableNames{'GoodPoint'} = 'GP';
end