FROM nlesc/magma

MAINTAINER Lars Ridder <l.ridder@esciencecenter.nl>

RUN conda install -c conda-forge chemspipy
ADD annotate_docs.py /motif_annotation/annotate_docs.py
ADD ms2lda_feature_extraction.py /ms2ldaviz/lda/code/ms2lda_feature_extraction.py

WORKDIR /
ENTRYPOINT ["/motif_annotation/annotate_docs.py"]
