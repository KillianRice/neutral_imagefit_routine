function [ output_args ] = Red_Spec_Power(  )
%Red_Spec_Power Takes Picoscope trace of the pulsed 689 spectroscopy beam,
%and determines it's amplitude.
%  

intenFilename = ['data.pico'];
intenDataFile = [intenFilename];
tmpIntVec     = dlmread(intenDataFile);
       
scopeTime  = tmpIntVec(1,:);
inten      = tmpIntVec(2,:);
           
pulseDur         = 200; % Time in ms for the PAS beam to be on. This comes in from the scan parameters every shot.

intenMeanBG      = mean(inten(scopeTime < pulseDur/1000));
intenMeanPulse   = mean(inten(scopeTime > pulseDur/1000));
       
       

%% Plotting

plot(scopeTime,inten,'b-','LineWidth',0.5)
%plot(time(time>0.000001),inten(time>0.000001),'b-','LineWidth',0.5)
str   = sprintf('Scope Trace');
title(str);
ylabel('Voltage (V)','FontSize',15,'FontWeight','bold');
xlabel('Time (s)','FontSize',15,'FontWeight','bold');
grid on
set(gca,'FontSize',23)
set(gca,'LineWidth',1)

hold on

upperRefline = refline([0 intenMeanBG]);
upperRefline.Color = 'r';

hold on

upperRefline = refline([0 intenMeanPulse ]);
upperRefline.Color = 'b';





end

