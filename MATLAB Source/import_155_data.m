function [ output ] = import_155_data( filename )
%IMPORT_155_DATA Summary of this function goes here
%   Detailed explanation goes here
fprintf('Reading data from file...\n');
%
% Import LVM File
%
lvm_data = lvm_import(filename, false);
fprintf('Finished reading. Consolidating data...\n');
%
% Consolidate data into arrays
%
fields = fieldnames(lvm_data);
n = length(fields);
running_data = [];
for i = 6:n
    f = fields{i};
    new_data = lvm_data.(f).('data');
    running_data = [running_data; new_data];
end
output = running_data;
fprintf('Operation complete.\n');
end