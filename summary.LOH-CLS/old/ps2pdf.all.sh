for fl in $( ls *.ps )
do
 echo "ps2pdf $fl"
 ps2pdf $fl
done
exit 0
