function [cellNames, colors] = assignNamesAndColors
cellNames = [{'HS_578T'} {'RPMI_8226'} {'HT29'} {'MALME_3M'} {'SR'}...
           {'UO_31'} {'MDMAMB_231'} {'HOP62'} {'NCI_H226'} {'HOP92'}...
           {'O_786'}];

colors    = [0    0    128             %dark blue
             0    0    255             %blue
             0    170  255             %light blue
             0    128  0               %forest green
             128  255  0               %lime green
             255  105  180             %pink
             255  200  050             %yellow
             191  0    255             %purple
             255  128  0               %orange
             0    0    0               %black
             255  0    0]./255;        %red
end     