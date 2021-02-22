function [colors, colors2] = getColorPalettes()
    colors = zeros(0, 3, 'double');
    colors(size(colors, 1)+1,:) = [201, 114, 92]/255; % [217, 77, 77]/255; % red
    colors(size(colors, 1)+1,:) = [229, 195, 115]/255; %[232, 193, 74]/255; % mustard
    colors(size(colors, 1)+1,:) = [132, 168, 122]/255; %[113, 209, 82]/255; % green
    colors(size(colors, 1)+1,:)= [120, 144, 190]/255; %[111, 157, 237]/255; % cornflower blue
    colors(size(colors, 1)+1,:)= [166, 136, 170]/255; %[169, 121, 212]/255; % purple
    colors(size(colors, 1)+1,:) = [226, 179, 161]/255; %[237, 168, 168]/255; % pink-red
    colors(size(colors, 1)+1,:)= [195, 213, 159]/255; %[186, 245, 128]/255; % light green
    colors(size(colors, 1)+1,:) = [222, 235, 206]/255; %[219, 255, 191]/255; % lightest green
    colors(size(colors, 1)+1,:) = [240, 221, 177]/255; %[242, 220, 160]/255; % lighter mustard
    colors(size(colors, 1)+1,:)= [253, 239, 211]/255; %[255, 240, 204]/255; % lightest mustard
    colors(size(colors, 1)+1,:)= [173, 190, 224]/255; %[173, 209, 255]/255; % blue
    colors(size(colors, 1)+1,:)=[204, 215, 236]/255; %[207, 230, 255]/255; % lightest blue
    colors(size(colors, 1)+1,:) = [238, 215, 203]/255;%[245, 212, 212]/255; % pink
    colors(size(colors, 1)+1,:) = [211, 197, 216]/255;%[218, 194, 242]/255; % light purple
    colors(size(colors, 1)+1,:)= [192, 181, 158]/255;%[191, 191, 153]/255; % khaki
    colors(size(colors, 1)+1,:) = [226, 221, 210]/255;%[230, 230, 214]/255; % beige
    colors(size(colors, 1)+1,:) = [132, 156, 144]/255;%[118, 180, 135]/255; % green 2
    colors(size(colors, 1)+1,:) = [209, 215, 211]/255;%[204, 217, 204]/255; % lighter green 2
    colors(size(colors, 1)+1,:) = [137, 86, 48]/255;%[132, 62, 14]/255; % brown
    colors(size(colors, 1)+1,:) = [204, 180, 152]/255; %[207, 177, 143]/255; % tan


    colors2(1,:) = [0.85, 0.325, 0.0980];    % red
    colors2(2,:) = [0.929, 0.694, 0.1250];% yellow
    colors2(3,:) = [0.4660, 0.674, 0.188]; % green
    colors2(4,:) = [0, 0.4470, 0.7410]; % blue
    colors2(5,:) = [0.3010, 0.745, 0.9330]; % light blue
    colors2(6,:) = [0.4940, 0.1840, 0.5560]; % purple
    colors2(7,:) = [0.6350, 0.078, 0.1840]; % burgundy?

end