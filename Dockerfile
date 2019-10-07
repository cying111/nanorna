FROM nfcore/base
LABEL authors="Laura Wratten" \
      description="Docker image containing all requirements for nf-core/nanorna pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-nanorna-1.0dev/bin:$PATH
