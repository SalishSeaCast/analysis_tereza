#!/bin/sh
# first ARG ($1) is working directory where files are and movie will be made
# ex: "~/pyCode/notebooks/figs/nowcastGreen"
WDIR=$1
# second ARG ($2) is base of filename, eg "int" or "surf" or "thw", where files are named int*.png etc
# links  will be mPATxxxxx.png, output file will be movie_PAT.avi
# ex: "integrated"
PAT=$2

FRATE=$3
# change directories to WDIR:
cd $WDIR

# delete existing symbolic links matching pattern in working directory so that link creation won't fail due to past runs
find -type l -name m${PAT}"*.png" -delete

# delete movie if already exists
rm movie_${PAT}.mp4

# create links
i=0
for f in `ls ${PAT}*.png`;
do
        ln -s $f $(printf "m"${PAT}"%05d.png" $i);
        i=$((i+1));
done 

# use links to create movie
avconv -framerate $FRATE -f image2 -i m${PAT}%05d.png -qscale 2 movie_${PAT}.mp4
