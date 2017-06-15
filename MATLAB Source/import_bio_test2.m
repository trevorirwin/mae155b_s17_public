function [adata_out, cdata_out] = import_bio_test2()
%IMPORT_RP1_TEST Summary of this function goes here
%   Detailed explanation goes here
adata = import_155_data('analog-test-data_2.lvm');
cdata = import_155_data('counter_2.lvm');
%
% Correct counter data fuel K factor
%
cdata(:,7) = cdata(:,3)/2222;
cdata(:,9) = cdata(:,9).*735./2222;
%
% Perform timestamp correction
% Here we take the difference between the data recording timestamps
% And add it to one set of time columns so that T = 0 is synchronized.
%
analog_timestamp = datevec('16:12:49.85871');
counter_timestamp = datevec('16:12:49.95715');
e_offset = etime(counter_timestamp, analog_timestamp);
adata(:,1) = adata(:,1) + e_offset;
%
% Duplicate counter data so we have the same number of indices with the
% same timestamps
%
cdata_long = zeros(size(adata));
crow = 1;
for arow = 1:length(adata)
    if (adata(arow,1) > cdata(crow))
        crow = crow + 1;
    end
    cdata_long(arow,:) = cdata(crow,:);
end
%
% Truncate the data to just the useful rows
%
ia = find(adata(:,10) ~= Inf & adata(:,10) ~= -Inf & adata(:,10) ~= 0);
ic = find(cdata_long(:,10) ~= 0);

imin = min(ia(1), ic(1));
imax = max(ia(end), ic(end));

iboth = [imin:1:imax];

adata_trunc = adata(iboth,:);
cdata_trunc = cdata_long(iboth,:);

adata_out = adata_trunc;
cdata_out = cdata_trunc;
%
% End of function
%
end

