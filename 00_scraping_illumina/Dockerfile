from ccdlstaging/dr_illumina:latest

COPY scrape_illumina_duplicate_probe_ids.R .

USER root
RUN mkdir /out && chown user /out 
USER user

ENTRYPOINT [ "Rscript", "scrape_illumina_duplicate_probe_ids.R" ]