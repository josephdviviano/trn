%% Polar Angle Tremotopy for PTB 3 %% Makes Eyes Bleed and Brains Sing %%%%%%%%%
% Retinotopy - K Schneider 12/13/05, KbWait(-1) and KbCheck(-1) 2/3/10
% Flicker    - J Viviano   29/11/12 w/ frame-based timing

try

clc;
AssertOpenGL;            % Check OpenGL is avaliable
Priority = 100;          % Set High Priority
KbName('UnifyKeyNames'); % Portable keyboard handling (win / linux / OSX)
ivx=iViewXInitDefaults;  % Init eyetracker

%% Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timing
opt.RF_period = 21;	     % period (sec) for rotating hemifield (cosine)
opt.FF_period = 30;      % period (sec) for flicker rate (staircase)
    
opt.RF_cycles = 10;      % number of cycles for rotating hemifield
opt.FF_cycles = 7;

opt.rewind = 10;         % rewind of flickulus in secs (allow for steady state)
opt.rest = 0;           % blank period at end of stimulus in secs

% flicker frequencies
opt.HzProjector = 120;           % refresh rate of projector (theoretical)
opt.HzList = [2,5,10,15,20,24,30,40,60,120]; % frequencies to test

% eye tracker network settings
ivx.host = '192.168.0.1';
ivx.port = 4444;
ivx.udp = '127.0.0.1';
ivx.localport = 5555;

% misc options
opt.whichScreen = 0;     % normally 0 for the primary display
opt.shrinkFactor = 6;    % scales the image for unusual projection situations
opt.print = 1;           % print a rendering of the stimulus (.png)

% NB: All periods should not divide evenly into eachother, but should still 
%     allow for an interger number of cycles in the full scan time. 
% NB: iViex network settings start on line 185
  
%% Check Settings for Human Errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure user defined the right number of valid frequencies for the projector
Hz.all = 1:opt.HzProjector;
Hz.listAll = find(rem(opt.HzProjector, Hz.all) == 0);
Hz.test = setdiff(opt.HzList, Hz.listAll);

% abort the program if user makes a mistake
if isempty(Hz.test) == 0;
    badFreqs = sprintf('invalid frequency: %d \n', Hz.test);
    disp(badFreqs);
    break
end

if rem(opt.FF_period, length(opt.HzList)) > 0;
    fprintf('length of frequency list must divide into flicker period');
    break
end

% fix shrinkFactor if user is a smooth brain
if opt.shrinkFactor < 1; 
    opt.shrinkFactor = 1;  
end

%% Generate Stimulus and Flip Vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define escape fxn and open Screen
esckey = KbName('Escape');
[w, r] = Screen('OpenWindow', opt.whichScreen, 0, [], 32, 2, [], [], 1);
rect = r/opt.shrinkFactor;
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
HideCursor;

% calculate flicker frequency vector for the run
time.total = opt.FF_period * opt.FF_cycles;
vec.freq = 1:time.total;
numSecPerHz = opt.FF_period/length(opt.HzList);

% create frequency vector
counter.HzSample = 1;   % position in list
counter.HzSecond = 1;   % runs from 1 - time.total

% create frequency vector for a single cycle
vec.freq = 1:opt.FF_period;

for thePos = 1:length(vec.freq);
    vec.freq(thePos) = opt.HzList(counter.HzSample);
    counter.HzSecond = counter.HzSecond + 1;
    if counter.HzSecond > numSecPerHz;
        counter.HzSecond = 1;
        counter.HzSample = counter.HzSample + 1;
    end
end

% create a frame vector
flick = 1;                 % 1 = phase, 2 = counterphase
counter.freq = 1;         % freqs  / run
counter.secs = 1;
thePosTarget = opt.HzProjector/vec.freq(counter.freq); % frames / sec

% frame vector from actual projector refresh rate for a single cycle
ifi = Screen('GetFlipInterval', w);
vec.flip = 1:length(vec.freq) * 1/ifi;

for thePos = 1:length(vec.flip);
    vec.flip(thePos) = flick;
    counter.secs = counter.secs + 1;
    
    if thePos > thePosTarget;
        thePosTarget = thePosTarget + (opt.HzProjector/vec.freq(counter.freq));
        
        flick = flick + 1;
        if flick > 2; 
            flick = 1; 
        end
        
        if counter.secs > opt.HzProjector;
            counter.freq = counter.freq + 1;
            counter.secs = 1;
        end
    end
end

% now REPMAT frame vector to length of full run
vec.flip = repmat(vec.flip, [1 opt.FF_cycles]);

% finaly add in rewind frames
vec.rewind = 1:opt.rewind * 1/ifi;
vec.rewind = vec.flip(length(vec.flip) - length(vec.rewind) : length(vec.flip));
vec.flick(1:length(vec.rewind)) = vec.rewind;
vec.flick(length(vec.rewind)+1:length(vec.rewind)+length(vec.flip)) = vec.flip;

%% Begin Stimulation of Human Eye %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Screen-dependent circle dimensions
xc = r(3)/2;
yc = r(4)/2;

% b/w circle pairs (8) / angular segments (24) must be integral
if opt.shrinkFactor > 1
    checkerCycles.r = 8/2;
    checkerCycles.t = 24/2;
else
    checkerCycles.r = 8;
    checkerCycles.t = 24;
end	

% make flickulus
index.white=WhiteIndex(w);
index.black=BlackIndex(w);
index.hi = index.white;
index.lo = index.black;
index.bg = (index.white + index.black) / 2;

xysize = rect(4);
s = xysize/sqrt(2);
xylim = 2*pi*checkerCycles.r;

[x,y] = meshgrid(-xylim:2*xylim/(xysize-1):xylim, ...
                 -xylim:2*xylim/(xysize-1):xylim);
at = atan2(y,x);

checks = ((1+sign(sin(at*checkerCycles.t)+eps) .* sign(sin(sqrt(x.^2+y.^2))))/2) ...
       * (index.hi-index.lo) + index.lo;

circle = x.^2 + y.^2 <= xylim^2;
checks = circle .* checks + index.bg * ~circle;

% Create checkerboard / reversed contrast
t(1) = Screen('MakeTexture', w, checks);               % phase
t(2) = Screen('MakeTexture', w, index.white - checks); % counterphase

% Create Fixation Point
Screen('FillRect', w, index.black);
Screen('FillRect', w, index.bg, [xc-rect(3)/2 yc-rect(4)/2 ...
                                 xc+rect(3)/2 yc+rect(4)/2]);
%if opt.shrinkFactor <= 2;
    Screen('FillRect', w, index.black, [xc-3 yc-3 xc+3 yc+3]);
    Screen('FillRect', w, index.white, [xc-2 yc-2 xc+2 yc+2]);
    Screen('FillRect', w, index.black, [xc-1 yc-1 xc+1 yc+1]);
%else
%    Screen('FillRect', w, [255,0,0], [xc-1 yc-1 xc+1 yc+1]);
%end

% Draw Text & Wait for Keypress
txt = 'Fixate';
if opt.shrinkFactor <= 2;
    Screen('TextSize', w, 24);
else
    Screen('TextSize', w, 14);
end
normBoundsRect = Screen('TextBounds', w, txt);
txtloc = [xc - normBoundsRect(3)/2, yc + normBoundsRect(4)/2];
Screen('DrawText', w, txt, txtloc(1), txtloc(2), index.white);
Screen('Flip', w);

KbWait(-1);

%% Animation Loop: Now With Eye Tracking %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check connection with the iView, exit program if this fails
%if iViewX('initialize', ivx)~=1;
%	return;
% end;

% Calibrate and drift correct (what does this do??)
%iViewX('calibration', ivx);
%iViewX('driftcorrection', ivx);

time.now = -opt.rewind;

iViewX('startrecording', ivx);
iViewX('message', ivx, 'Start of display');
time.start = GetSecs;
theFrame = 0;

% animation loop
while time.now <= opt.RF_period * opt.RF_cycles;  
    
    % draw the next frame's full field checkerboard
    if theFrame ~= length(vec.flick);
        theFrame = theFrame + 1;
    end
    
    flick = vec.flick(theFrame);
    Screen('FillRect', w, index.bg, [xc-rect(3)/2 yc-rect(4)/2 ...
                                     xc+rect(3)/2 yc+rect(4)/2]); 
    Screen('DrawTexture', w, t(flick));
     
    % draw rectangular mask for rotation
    theta = mod(time.now, opt.RF_period)/opt.RF_period*2*pi;
    st = sin(theta);
    ct = cos(theta);
    xy = s * [-st,-ct; -st-ct, -ct+st; st-ct, ct+st; st, ct];
    xy = xy + ones(4,1) * [xc yc];
    Screen('FillPoly', w, index.bg, xy);
    Screen('FillRect', w, index.black, [xc+rect(3)/2 r(2) r(3) r(4)]);
    Screen('FillRect', w, index.black, [r(1) yc+rect(4)/2 r(3) r(4)]);
    Screen('FillRect', w, index.black, [r(1) r(2) xc-rect(3)/2 r(4)]);
    Screen('FillRect', w, index.black, [r(1) r(2) r(3) yc-rect(4)/2]);
    
    % draw fixation point
    %if opt.shrinkFactor <= 2;
        Screen('FillRect', w, index.black, [xc-3 yc-3 xc+3 yc+3]);
        Screen('FillRect', w, index.white, [xc-2 yc-2 xc+2 yc+2]);
        Screen('FillRect', w, index.black, [xc-1 yc-1 xc+1 yc+1]);
    %else
    %    Screen('FillRect', w, [255,0,0], [xc-1 yc-1 xc+1 yc+1]);
    %end
    
    % abort if keypress
    [kdown,secs,keyCode] = KbCheck(-1);
    if keyCode(esckey); break; end
    
    % execute a flip
    [time.VBL, time.flick, time.stamp, missCheck] = Screen('Flip', w);
    time.now = time.stamp - time.start - opt.rewind;
    
    % if we miss a deadline, skip to the next frame so we don't fall behind
    if missCheck > 0; 
        theFrame = theFrame + 1; 
    end
    
    % write out an image if in print mode
    if opt.print == 1;
        img = Screen('GetImage', w, [xc-rect(3)/2 yc-rect(4)/2 ...
                                     xc+rect(3)/2 yc+rect(4)/2]);
        imwrite(img, 'C:\Users\Paperbag\Desktop\stim.png')
    end
end

% rest period loop
while time.now <= opt.RF_period * opt.RF_cycles + opt.rest;
    
    % draw background
    Screen('FillRect', w, index.black);
    Screen('FillRect', w, index.bg, [xc-rect(3)/2 yc-rect(4)/2 ...
                                     xc+rect(3)/2 yc+rect(4)/2]);
    % draw fixation point
    %if opt.shrinkFactor <= 2;
        Screen('FillRect', w, index.black, [xc-3 yc-3 xc+3 yc+3]);
        Screen('FillRect', w, index.white, [xc-2 yc-2 xc+2 yc+2]);
        Screen('FillRect', w, index.black, [xc-1 yc-1 xc+1 yc+1]);
    %else
    %    Screen('FillRect', w, [255,0,0], [xc-1 yc-1 xc+1 yc+1]);
    %end
    
    % abort if keypress
    [kdown,secs,keyCode] = KbCheck(-1);
    if keyCode(esckey); break; end
    
    % execute a flip
    [time.VBL, time.flick, time.stamp, missCheck] = Screen('Flip', w);
    time.now = time.stamp - time.start - opt.rewind;
    
end

iViewX('stoprecording', ivx);
iViewX('message', ivx, 'End of display');
time.end = GetSecs;
time.total = time.end - time.start;

%% Go To Sleep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw background
Screen('FillRect', w, index.bg);
Screen('FillRect', w, index.black, [xc+rect(3)/2 r(2) r(3) r(4)]);
Screen('FillRect', w, index.black, [r(1) yc+rect(4)/2 r(3) r(4)]);
Screen('FillRect', w, index.black, [r(1) r(2) xc-rect(3)/2 r(4)]);
Screen('FillRect', w, index.black, [r(1) r(2) r(3) yc-rect(4)/2]);

% draw Text
txt = 'Thank you';
if opt.shrinkFactor <= 2;
    Screen('TextSize', w, 24);
else
    Screen('TextSize', w, 14 );
end
normBoundsRect = Screen('TextBounds', w, txt);
txtloc = [xc - normBoundsRect(3)/2, yc + normBoundsRect(4)/2];
Screen('DrawText', w, txt, txtloc(1), txtloc(2), index.white);
Screen('Flip', w);
   
KbWait(-1);

ShowCursor;
Screen('CloseAll');

catch
ShowCursor;
Screen('CloseAll');
rethrow(lasterror);
end
