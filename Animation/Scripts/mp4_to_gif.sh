#!/bin/bash



ffmpeg \
  -i $1.mp4 \
  -r 15 \
  -vf scale=512:-1 \
  -ss 00:00:$3 -to 00:00:$4 \
  $2.gif
