function [q_tunnel, T_tunnel, p_amb, ...
            calfile, tarefile, ...
            AoA_high, AoA_low, dF, goodData] = get_testCondition(runName)
    %GET_TESTCONDITION Return the info from the testlog.xlsx spreadsheet
    %corresponding to a given test run. Returns values as given in sheet -
    %units and outputs described below may not be valid if spreadsheet
    %changes. 
    %
    %   INPUTS
    %
    %   runName - name of the run e.g. "run_1001". Must be in row 1/
    %
    %   OUTPUTS
    %
    %   q_tunnel    Tunnel dynamic pressure, [Pa]. Spreadsheet column 6
    %   T_tunnel    Tunnel temp, [K] Spreadsheet column 8
    %   p_amb       Ambient pressure, [Pa]. Spreadsheet column 7
    %   dF          Flap deflection, [deg]. Spreadsheet column 3
    %   AoA_high    Angle of attack, high range [deg]. Spreadsheet column 5
    %   AoA_low     Angle of attack, low range [deg]. Spreadsheet column 4
    %   calfile     Name of file with correct calibration matrix.
    %               Spreadsheet column 9
    %   tarefile    Name of file with correct tare values. Spreadsheet
    %               column 10
    %   goodData    Y/N, whether datapoint is good. Spreadsheet column 11
    opts = spreadsheetImportOptions("NumVariables", 11);
    opts.Sheet = 1;

    data = readtable("testlog.xlsx", opts, "UseExcel", false);
    
    idx = any(ismember(data{:,1},{char(runName)}),2);
    
    dF          = str2num(data.Var3{idx});
    AoA_low     = str2num(data.Var4{idx});
    AoA_high    = str2num(data.Var5{idx});
    q_tunnel    = str2num(data.Var6{idx});
    p_amb       = str2num(data.Var7{idx});
    T_tunnel    = str2num(data.Var8{idx});
    calfile     = data.Var9{idx};
    tarefile    = data.Var10{idx};
    goodData    = data.Var11{idx};
    
    
    

end

