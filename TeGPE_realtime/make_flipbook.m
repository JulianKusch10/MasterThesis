ffmpegPath = 'C:\ffmpeg-8.0-essentials_build\bin\ffmpeg.exe';
directory = 'C:\Realtime_simulations\TeGPE_realtime_comparison_5_110_110';
input_files = 'TeGPE_realtime_comparison_5_110_110_%03d.png';
output_file = 'flipbook_110_110.mp4';

fps = 3;

output_file_path = fullfile(directory, output_file);
%output_file_path = sprintf('%s/%s', directory, output_file);

if isfile(output_file_path)
    disp('Output file already exists. Overwrite? [y/n]: ');
    answer = input('', 's'); 
    if answer == 'y'
        cd(directory)
        cmd = sprintf('%s -y -framerate %d -i %s -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p %s', ffmpegPath, fps, input_files, output_file);
        status = system(cmd);
    elseif answer == 'n'
        disp("Stopping, please change the output filename.")
        return 
    end
else
    cd(directory)
    cmd = sprintf('%s -y -framerate %d -i %s -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p %s', ffmpegPath, fps, input_files, output_file);
    status = system(cmd);
end