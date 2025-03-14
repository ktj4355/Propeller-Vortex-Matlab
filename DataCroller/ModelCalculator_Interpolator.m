%Create Interpolation Table 

close all
clear
clc
%Section1 Model
listing=dir("Section1");
DB_Aero={}; %Alpha Cl Cd
DBind=[]; %RE Mach
Total_DB=[];
% matlab dir 함수에서는 3번째 줄부터 실제 폴더이름
for idx=3:size(listing,1)
    

    % 데이터 예시 (아래줄부터 시작)
    % 파일 이름 : E63_4_27_T1_Re0.050_M0.20_N7.0.txt

    %
    % xflr5 v6.61
    %
    % Calculated polar for: E63_4_27
    %
    % 1 1 Reynolds number fixed          Mach number fixed
    %
    % xtrf =   1.000 (top)        1.000 (bottom)
    % Mach =   0.200     Re =     0.050 e 6     Ncrit =   7.000
    %
    %  alpha     CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge    XCp
    % ------- -------- --------- --------- -------- ------- ------- -------- --------- ---------
    % -24.000  -0.7585   0.29481   0.28669   0.0419  1.0000  0.0225  -1.5489   0.0000   0.0000   0.0000   0.2964
    %    
    fileName=listing(idx).name;
    filePath=listing(idx).folder;
    fullName=[filePath '\' fileName];
    NameSplit=split(fileName,'_');
    AirfoilName=NameSplit{1};
    thick_major=str2num(NameSplit{2});
    thick_minor=str2num(NameSplit{3});
    thick=thick_major+thick_minor*0.01;

    
    AeroDB_Local_line=readlines(fullName);
    Setting=double(split(AeroDB_Local_line(8)));
    Local_M=Setting(4);
    Local_Re=Setting(7)*(10^Setting(9));
    AeroDB_Local=double(split(AeroDB_Local_line(12:end-3)));
    AeroDB_Local(:,1)=[];
    DBind=[DBind;Local_Re,Local_M];
    DB_AeroLocal=[ones(size(AeroDB_Local,1),1)*thick ones(size(AeroDB_Local,1),1)*Local_Re  AeroDB_Local(:,1:3)];
    Total_DB=[Total_DB;DB_AeroLocal];
    DB_Aero=[DB_Aero;DB_AeroLocal];

    

end
Total_DB;
% Total_DB : [thick, RE, AoA, cl, cd]
Section1_CL=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,4),"linear","nearest");
Section1_CD=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,5),"linear","nearest");

sampleThk=10;
sampleRE=50000;
sampleAOA=0;

cl=Section1_CL(sampleThk,sampleRE,sampleAOA);

%% Section2
listing=dir("Section2");
listing(1);
DB_Aero={}; %Alpha Cl Cd
DBind=[]; %RE Mach
Total_DB=[];
for idx=3:size(listing,1)
    fileName=listing(idx).name;
    filePath=listing(idx).folder;
    fullName=[filePath '\' fileName];
    NameSplit=split(fileName,'_');
    AirfoilName=NameSplit{1};
    thick_major=str2num(NameSplit{2});
    thick_minor=str2num(NameSplit{3});
    thick=thick_major+thick_minor*0.01;


    %  alpha     CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge    XCp
    AeroDB_Local_line=readlines(fullName);
    Setting=double(split(AeroDB_Local_line(8)));
    Local_M=Setting(4);
    Local_Re=Setting(7)*(10^Setting(9));
    AeroDB_Local=double(split(AeroDB_Local_line(12:end-3)));
    AeroDB_Local(:,1)=[];
    DBind=[DBind;Local_Re,Local_M];
    DB_AeroLocal=[ones(size(AeroDB_Local,1),1)*thick ones(size(AeroDB_Local,1),1)*Local_Re  AeroDB_Local(:,1:3)];
    Total_DB=[Total_DB;DB_AeroLocal];
    DB_Aero=[DB_Aero;DB_AeroLocal];
end

Section2_CL=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,4),"linear","nearest");
Section2_CD=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,5),"linear","nearest");




%% blend Fucntion

listing=dir("Blended");
listing(1);
DB_Aero={}; %Alpha Cl Cd
DBind=[]; %RE Mach
Total_DB=[];
for idx=3:size(listing,1)
    fileName=listing(idx).name;
    filePath=listing(idx).folder;
    fullName=[filePath '\' fileName];
    NameSplit=split(fileName,'_');
    AirfoilName=NameSplit{1};
    BlendPosition=str2num(NameSplit{2});
    BlendPosition=BlendPosition./100;

    %  alpha     CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge    XCp
    AeroDB_Local_line=readlines(fullName);
    Setting=double(split(AeroDB_Local_line(8)));
    Local_M=Setting(4);
    Local_Re=Setting(7)*(10^Setting(9));
    AeroDB_Local=double(split(AeroDB_Local_line(12:end-3)));
    AeroDB_Local(:,1)=[];
    DBind=[DBind;Local_Re,Local_M];
    DB_AeroLocal=[ones(size(AeroDB_Local,1),1)*BlendPosition ones(size(AeroDB_Local,1),1)*Local_Re  AeroDB_Local(:,1:3)];
    Total_DB=[Total_DB;DB_AeroLocal];
    DB_Aero=[DB_Aero;DB_AeroLocal];
end


BLD_CL=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,4),"linear","nearest");
BLD_CD=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,5),"linear","nearest");

ModelTitle=sprintf("InterpolatedModel");
save(ModelTitle,'Section1_CL','Section1_CD','Section2_CL','Section2_CD','BLD_CL','BLD_CD')
