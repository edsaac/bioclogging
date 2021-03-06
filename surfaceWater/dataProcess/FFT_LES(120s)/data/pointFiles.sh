## How to use:
## ./pointFiles.sh <pathToFolder> <nPointsProbed>

FOLDER=$1
LINES=$2

# Note to self: Remember to chech the lines, it aint always 98

#mkdir ./organizedData/

#for LINE in $(seq -s ' ' 1 $LINES)
for LINE in $LINES
do
  resultFile="probeAt_$(printf "%04d" $LINE).txt"
  
  for TIME in $(ls -v $FOLDER)
  do
    # currentFile="./$FOLDER/$TIME/line_U.xy"
	# echo $TIME $(sed "$LINE!d" $currentFile | \
	# sed -e "s/\s\+/ /g" | \
	# cut -d " " -f 2-4) >> \
	# ./organizedData/$FOLDER/$resultFile
	
	currentFile="./$FOLDER/$TIME/points_p.xy"
	echo $TIME $(sed "$LINE!d" $currentFile | \
	sed -e "s/\s\+/ /g" | \
	cut -d " " -f 4) >> \
	./organizedData/$resultFile

	#echo $TIME
  done
done