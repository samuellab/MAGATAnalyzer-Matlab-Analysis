function writeFrameOrPause(vidObj, avirect, fps, nframes)
%function writeFrameOrPause(vidObj, avirect, fps, nframes)

existsAndDefault('fps', 25);
existsAndDefault('nframes', 1);
existsAndDefault('avirect', []);
for j = 1:nframes
    if (~isempty(vidObj))
       if (isempty(avirect))
           currFrame = getframe(gcf);
       else
           currFrame = getframe(gcf,avirect);
       end
      % disp('writing a frame');
       writeVideo(vidObj,currFrame);
    else
        pause(1/fps);
    end
end