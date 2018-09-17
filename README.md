# motif_annotation
Annotate ms2lda motifs using MAGMa

# Dependencies
- python2
- sqlalchemy
- rdkit
- chempipy
- magma

# Install / Run annotate_motifs.py
```
git clone https://github.com/iomega/motif_annotation.git
git clone https://github.com/sdrogers/motifdb.git
mkdir gnps
cp [path]/MSMS-GNPS-Curated-Pos-MfKit.msp gnps
cp [path]/MS2LDA.csv gnps

python motif_annotation/annotate_motifs.py gnps/MSMS-GNPS-Curated-Pos-MfKit.msp gnps/gnps.db 'Links:' gnps/MS2LDA.csv motifdb/motifs/gnps_binned_005/gnps_ >annotation.json
```
# Run create_html_view.py
```
wget https://web.chemdoodle.com/downloads/ChemDoodleWeb-8.0.0.zip
unzip ChemDoodleWeb-8.0.0.zip
cp ChemDoodleWeb-8.0.0/install/ChemDoodleWeb.js .

python motif_annotation/create_html_view.py <annotation.json >annotattion.html
```
