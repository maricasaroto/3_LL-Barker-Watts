#!/bin/bash
 
# created     : 2018/06/20
# last update : 2018/06/20
# author      : Mariana Casaroto
# notes       : cria um vídeo com as imagens


ffmpeg -r 6 -i %d.png -c:v libx264 -profile:v baseline -level 3.0 -pix_fmt yuv420p -y video.mp4
 
