function makeEthogram(b, color, fs, yrange, opacity)
if nargin < 5, opacity = 0.75; end
time = xt(b,fs);
[bOn, bOff] = helper.getBehaviourEvents(b, 0);
x = [time(bOn); time(bOn); time(bOff); time(bOff)];
y = repmat([yrange(1); yrange(2); yrange(2); yrange(1);],[1, size(x,2)]);
patch(x, y, color, 'FaceAlpha', opacity, 'EdgeColor', 'None'), hold on

end

