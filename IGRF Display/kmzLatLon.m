function [LAT, LON] = kmzLatLon(filename)
% returns kmz or kml file lat lon points

KS = kmz2struct(filename);
TS = struct2table(KS); clear KS

LAT = []; LON = [];
for idx = 1:height(TS)
    lat = TS.Lat{idx};
    lon = TS.Lon{idx};
    LAT = [LAT lat];
    LON = [LON lon];
end

IDX = or(LAT == -90, LAT == 90);
LAT(IDX) = []; 
LON(IDX) = [];

end