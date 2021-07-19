from ccdlstaging/dr_affymetrix:latest

COPY scrape_brainarray_genes.R .

USER root
RUN mkdir /out && chown user /out 
USER user

ENTRYPOINT [ "Rscript", "scrape_brainarray_genes.R" ]