function Homework2inputSLE020524 = importfile1(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  HOMEWORK2INPUTSLE020524 = IMPORTFILE(FILE) reads data from the first
%  worksheet in the Microsoft Excel spreadsheet file named FILE.
%  Returns the data as a table.
%
%  HOMEWORK2INPUTSLE020524 = IMPORTFILE(FILE, SHEET) reads from the
%  specified worksheet.
%
%  HOMEWORK2INPUTSLE020524 = IMPORTFILE(FILE, SHEET, DATALINES) reads
%  from the specified worksheet for the specified row interval(s).
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  Homework2inputSLE020524 = importfile("C:\Users\spenc\Desktop\school shit 2\Spring_2024\upper_atmosphere\hw2\Homework_2_input_SLE_020524.xlsx", "Sheet1", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 10-Feb-2024 17:56:46

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 9);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["AltitudeKm", "TemperatureK", "AltitudeKm1", "ArCm3", "AltitudeKm2", "N2Cm3", "EarthRadiusAtThisLocation", "VarName8", "km"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "string", "double", "string"];

% Specify variable properties
opts = setvaropts(opts, ["EarthRadiusAtThisLocation", "km"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["EarthRadiusAtThisLocation", "km"], "EmptyFieldRule", "auto");

% Import the data
Homework2inputSLE020524 = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    Homework2inputSLE020524 = [Homework2inputSLE020524; tb]; %#ok<AGROW>
end

end