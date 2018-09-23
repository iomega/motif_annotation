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
cp -r [path]/spectra_gnps gnps
cp [path]/MS2LDA.csv gnps

python motif_annotation/annotate_motifs.py gnps/spectra_gnps gnps/gnps.db 0.05 gnps/MS2LDA.csv motifdb/motifs/gnps_binned_005/gnps_ >gnps/annotation.json
```
# Run create_html_view.py
```
wget https://web.chemdoodle.com/downloads/ChemDoodleWeb-8.0.0.zip
unzip ChemDoodleWeb-8.0.0.zip
cp ChemDoodleWeb-8.0.0/install/ChemDoodleWeb.js .

python motif_annotation/create_html_view.py <gnps/annotation.json >gnps/annotattion.html
```
