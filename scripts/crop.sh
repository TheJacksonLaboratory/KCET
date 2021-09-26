files=`ls *.pdf`
for f in $files
do
    pdfcrop $f
    rm $f
done
files=`ls *.pdf`
for f in $files
do
    mv "$f"  "${f/-crop.pdf/.pdf}"
done
