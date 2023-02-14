ffmpeg -framerate 4 -i image%04d.png -c:v libx264 -crf 18 -r 24 output.mp4
