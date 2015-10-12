
#
echo Downloading data...
mkdir data
cd data
python ../scripts/data_download.py

#
echo Decompressing files...
gunzip -vzxf *.gz

#
echo Converting file formats...
python ../scripts/gff2bed.py *.gff3
python ../scripts/mid_finder.py *.bed

#
echo Creating figure...
cd ..
python beaf_figure.py data/*.bed.mid

