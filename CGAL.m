%clear;
clc
close all;

figure;
lines = readFromTxt('./cmake-build-debug/CGAL_Line_Output.txt');
edges = readFromTxt('./cmake-build-debug/CGAL_Wall_Output.txt');



for i = 1:length(lines(1,:))
    line = lines(:,i);
    plot([line(1),line(3)],[line(2),line(4)],'LineWidth',0.1);
    hold on;
end

for i = 1:length(edges(1,:))
    edge = edges(:,i);
    plot([edge(1),edge(3)],[edge(2),edge(4)], 'b','LineWidth',2);
    hold on;
end
minX = -3.4815;
maxX = 7.24;
minY = -8.47;
maxY = 0.6095;
plot([minX ,minX],[minY ,maxY],'k','LineWidth',2); hold on;
plot([minX ,maxX],[maxY ,maxY],'k','LineWidth',2); hold on;
plot([maxX ,maxX],[maxY ,minY],'k','LineWidth',2); hold on;
plot([maxX ,minX],[minY ,minY],'k','LineWidth',2); hold on;
ylim([-10 2])
xlim([-5 10])
readFromTxtCells('./cmake-build-debug/CGAL_Face_Output.txt');
%% CGAL PART
cgal_vertex = readFromTxt('./cmake-build-debug/CGAL_Vertex_Output.txt');
for i = 1:length(cgal_vertex(1,:))
    vertex = cgal_vertex(:,i);
    scatter(vertex(1,:),vertex(2,:), 5, 'r','filled');
    hold on;
end

function cords = readFromTxt(path)
    fid = fopen(path);
    cords = [];
    tline = fgetl(fid);
    while ischar(tline)
        if tline == -1
            break;
        end
        line = cellfun(@str2num,split(tline));
        cords = [cords line];
        tline = fgetl(fid);
    end
    fclose(fid);
end

function readFromTxtCells(path)
    k = 1;
    n = 1;
    fid = fopen(path);
    tline = fgetl(fid);
    %figure;
    while ischar(tline)
        k = k + 1;
        if tline == -1
            break;
        end
        splitLine = split(tline);
        line = cellfun(@str2num, splitLine(1:length(splitLine)));
        
        if (length(line) < 6)
            length(line);
        end
        warning('');
        shape = polyshape(line(1:2:end-1),line(2:2:end-1));
        if (line(end) == 1)
            plot(shape, 'FaceColor', 'red'); hold on;
        else
            plot(shape, 'FaceColor', 'yellow'); hold on;
        end
        
        [warnmsg, msgid] = lastwarn;
        if strcmp(msgid,'MATLAB:polyshape:repairedBySimplify')
            n = n+1;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    n
end