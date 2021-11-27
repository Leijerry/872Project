% Besed on "3D-face-nose-tip-recognition"
% URL="https://github.com/RainMen1277/3D-face-nose-tip-recognition"
% Modified by Lei Jiang

clc;
input_path = 'D:\facedataset1\pending\';  
output_path = 'D:\facedataset2\result\';  
fileExt = '*.xyz';
files = dir(fullfile(input_path,fileExt));
len = size(files,1);
for file_num = 1:len
    fileName = strcat(input_path,files(file_num,1).name);
    read_file = importdata(fileName);
    file_x = read_file(:,1);
    file_y = read_file(:,2);
    file_z = read_file(:,3);
read_file = importdata('D:\facedata1'); 
file_x = read_file(:,1);
file_y = read_file(:,2);
file_z = read_file(:,3);
set(figure,'name','pointcloud');
plot3(file_x,file_y,file_z,'.')
[X,Y,Z] = griddata(file_x,file_y,file_z,linspace(min(file_x),max(file_x),200)',linspace(min(file_y),max(file_y),20),'cubic'); 
set(figure,'name','grid');
surf(Z);  
nose = []; 
x_gaos = []; 
II = [];
for j_new = 1:20 
    xz = Z(j_new,:);
    I = find(~isnan(xz)); 
    [x_s,~] = min(I); 
    [x_e,~] = max(I); 
    %r = round((x_e-x_s)/5); 
    r = 30; 
    h_gaos = []; 
    for x1 = x_s+r:x_e-r 
        if x1>1 
            z1 = xz(1,x1); 
            i = x1-r:x1+r; 
            y2 = sqrt(r.^2-(i-x1).^2)+z1; 
            y = xz(1,i);
            cy = y2-y;
            pos = cy>0;
            neg = cy<=0;
            fro = diff(pos)~=0; 
            rel = diff(neg)~=0; 
            zpf = find(fro==1); 
            zpr = find(rel==1)+1; 
            zpfr = [zpf,zpr];
            x0 = (i(zpr).*(y2(zpf)-y(zpf))-i(zpf).*(y2(zpr)-y(zpr)))./(y(zpr)+y2(zpf)-y(zpf)-y2(zpr));
            y0 = y(zpf)+(x0-i(zpf)).*(y(zpr)-y(zpf))./(i(zpr)-i(zpr)-i(zpf));
            x0 = [x0 x0].';
            y0 = [y0 y0].';
            jie = unique([x0,y0],'rows'); 
            [mm,nn] = size(jie);
            %fprintf('number of solution：');
            %disp(mm);
            if mm == 2;
                jie1 = jie(1,:);
                jie2 = jie(2,:);
                if jie1(:,1)<x1 && jie2(:,1)>x1 
                    l1 = ((jie1(:,1)-x1)^2+(jie1(:,2)-z1)^2)^(1/2);
                    l2 = ((jie2(:,1)-x1)^2+(jie2(:,2)-z1)^2)^(1/2);
                    l3 = ((jie1(:,1)-jie2(:,1))^2+(jie1(:,2)-jie2(:,2))^2)^(1/2); 
                    p = (l1+l2+l3)/2;
                    s = sqrt(p*(p-l1)*(p-l2)*(p-l3)); 
                    h_gao = 2*s/l3;  
                    h_gaos = [h_gaos;[x1,h_gao]]; 
                end
            end
        end
    end
    if length(h_gaos)>0
        [m,p] = max(h_gaos(:,2));
        x_gao = h_gaos(p,1); 
        x_gaos = [x_gaos;[j_new,x_gao,m]]; 
    end
end
[max_h_1,max_p_1] = max(x_gaos(:,3));
j_max_1 = x_gaos(max_p_1,1);
j_2_s = Y((j_max_1-1),1);
j_2_e = Y((j_max_1+1),1);
number = size(read_file,1);
new_file = [];
for i = 1:number
    if j_2_s < read_file(i,2) && read_file(i,2) < j_2_e
        new_file = [new_file;read_file(i,:)];
    end
end
new_x = new_file(:,1);
new_y = new_file(:,2);
new_z = new_file(:,3);
[Xx,Yy,Zz] = griddata(new_x,new_y,new_z,linspace(min(new_x),max(new_x),200)',linspace(min(new_y),max(new_y),40),'cubic'); 
 
 for j_new = 1:40 
    xz = Zz(j_new,:); 
    I = find(~isnan(xz)); 
    [x_s,~] = min(I); 
    [x_e,~] = max(I); 
    %r = round((x_e-x_s)/5); 
    r = 30; 
    h_gaos = []; 
    for x1 = x_s+r:x_e-r 
        if x1>1 
            z1 = xz(1,x1); 
            i = x1-r:x1+r; 
            y2 = sqrt(r.^2-(i-x1).^2)+z1; 
            y = xz(1,i);
            cy = y2-y;
            pos = cy>0;
            neg = cy<=0;
            fro = diff(pos)~=0; 
            rel = diff(neg)~=0; 
            zpf = find(fro==1); 
            zpr = find(rel==1)+1; 
            zpfr = [zpf,zpr];
            x0 = (i(zpr).*(y2(zpf)-y(zpf))-i(zpf).*(y2(zpr)-y(zpr)))./(y(zpr)+y2(zpf)-y(zpf)-y2(zpr));
            y0 = y(zpf)+(x0-i(zpf)).*(y(zpr)-y(zpf))./(i(zpr)-i(zpr)-i(zpf));
            x0 = [x0 x0].';
            y0 = [y0 y0].';
            jie = unique([x0,y0],'rows'); 
            [mm,nn] = size(jie);
            %fprintf('number of solution：');
            %disp(mm);
            if mm == 2;
                jie1 = jie(1,:);
                jie2 = jie(2,:);
                if jie1(:,1)<x1 && jie2(:,1)>x1 
                    l1 = ((jie1(:,1)-x1)^2+(jie1(:,2)-z1)^2)^(1/2);
                    l2 = ((jie2(:,1)-x1)^2+(jie2(:,2)-z1)^2)^(1/2);
                    l3 = ((jie1(:,1)-jie2(:,1))^2+(jie1(:,2)-jie2(:,2))^2)^(1/2); 
                    p = (l1+l2+l3)/2;
                    s = sqrt(p*(p-l1)*(p-l2)*(p-l3)); 
                    h_gao = 2*s/l3;  
                    h_gaos = [h_gaos;[x1,h_gao]]; 
                end
            end
        end
    end
    if length(h_gaos)>0
        [m,p] = max(h_gaos(:,2));
        x_gao = h_gaos(p,1); 
        x_gaos = [x_gaos;[j_new,x_gao,m]]; 
    end
end
ans = sortrows(x_gaos,-3); 
nose_cut = ans(1:20,:); 
for mq = 1:20
    nose = [nose;[Xx(1,nose_cut(mq,2)),Yy(nose_cut(mq,1),1),Zz(nose_cut(mq,1),nose_cut(mq,2))]];
end
% fprintf('dimension：');
% disp(x_gaos);
% fprintf('point：');
% disp(nose_cut);
% fprintf('original location：');
% disp(nose);
hold on;
nose_x = nose(:,1);
nose_y = nose(:,2);
nose_z = nose(:,3);
plot3(nose_x,nose_y,nose_z,'o');
data=[nose(:,1)';nose(:,2)'];
iter = 500; 
number = size(data,2); 
bestParameter1=0; bestParameter2=0; 
sigma = 1;
pretotal=0;     
for i=1:iter
    idx = randperm(number,2); 
    sample = data(:,idx); 
    line = zeros(1,3);
    x = sample(1, :);
    y = sample(2, :);
    if x(1)==x(2)
        a = 1;
        b = 0;
        c = -x(1);
        line = [a b c];
    else
        k=(y(1)-y(2))/(x(1)-x(2)); 
        a = sqrt(1-1/(1+k^2));
        b = sqrt(1-a^2);
        if k > 0
            b = -b;
        end
        c = -a*x(1)-b*y(1);
        line = [a b c];
    end
    mask=abs(line*[data; ones(1,size(data,2))]);    
    total=sum(mask<sigma);              
    if total>pretotal            
        pretotal=total;
        bestline=line;          
    end  
end
mask=abs(bestline*[data; ones(1,size(data,2))])<sigma;    
hold on;
k=1;
nose_z = 1000;
for i=1:length(mask)
    if mask(i)
        if nose(i,3) < nose_z
            nose_z = nose(i,3);
            nose_i = i;
        end
    end
end
hold on;
nose_max_x = nose(nose_i,1);
nose_max_y = nose(nose_i,2);
nose_max_z = nose(nose_i,3);
plot3(nose_max_x,nose_max_y,nose_max_z,'g*');
output_file = strcat(output_path,files(file_num,1).name);
radius = 90;
[file_row,file_list] = size(read_file);
face = [];
for cut_i = 1:file_row
    face_x = read_file(cut_i,1);
    face_y = read_file(cut_i,2);
    face_z = read_file(cut_i,3);
    distance = ((face_x-nose_max_x)^2+(face_y-nose_max_y)^2+(face_z-nose_max_z)^2)^0.5;
    if distance < radius
        face = [face;read_file(cut_i,:)];
    end
end
xx = face(:,1);
yy = face(:,2);
zz = face(:,3);
set(figure,'name','人脸');
plot3(xx,yy,zz,'.');
%[face_row,face_list] = size(face);
output_file = 'C:\myfiles\xyz\test1\test2\face_1_1_1.xyz';
fid = fopen(output_file,'w');
for i = 1:length(face)
%     for j = 1:face_list
%         fprintf(fid,'%f\n',face(i,j));
%     end
    fprintf(fid,'%.3f %.3f %.3f\n',face(i,:));
%     fprintf(fid,'\r\n');
end
fclose(fid);
end

