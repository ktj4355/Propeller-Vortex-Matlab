clc
clear 
close all
%unloadlibrary('VortexDLL')
loadlibrary('VortexDLL.dll','Vortex_DLL.h')
% Define input values
A = [1.0, 0.0, 0.0];
B = [1.0, 1.0, 0.0];
ColocationPoint = [0.5, 0.5, 0.0];
vortexStrength = 1.0;
rc = 0.1;

% Create pointers
A_ptr = libpointer('doublePtr', A);
B_ptr = libpointer('doublePtr', B);
ColocationPoint_ptr = libpointer('doublePtr', ColocationPoint);
Vout_ptr = libpointer('doublePtr', zeros(1, 3));

% Call the C function
calllib('VortexDLL', 'Vortex_Scully', A_ptr, B_ptr, ColocationPoint_ptr, vortexStrength, rc, Vout_ptr);
Vout=Vout_ptr.Value;
% Get the output
VoutRef = Vortex_Scully(A,B,ColocationPoint,vortexStrength,rc);
disp('Vout:');
disp(Vout);

disp('VoutRef:');
disp(VoutRef);
%unloadlibrary('VortexDLL')