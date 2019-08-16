%moving window averaging
clear; clf
xy = 1.1123*rand(100,2);
x = xy(:,1);
y = xy(:,2);
window = 0.2;
step   = 0.025;
newx = [min(x): step : max(x)];
left = newx(1);
for i = 1:length(newx)
    right = left + window;
    cut = y(find((x>=left) & (x<right)));
    newy(i) = mean(cut);
    left = left + step;
end
figure(1)
plot(x, y, 'o');
hold on
plot(newx, newy, 's');