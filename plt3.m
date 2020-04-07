function plt3(data)
plot3(1:length(data),real(data),imag(data),'*');
xlabel('samples');
ylabel('real');
zlabel('imag');
end

