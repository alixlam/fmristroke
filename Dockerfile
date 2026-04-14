ARG FMRIPREP_VERSION=23.0.2
FROM nipreps/fmriprep:${FMRIPREP_VERSION}

LABEL org.opencontainers.image.title="fMRIStroke" \
      org.opencontainers.image.description="Container image for the fMRIStroke BIDS application" \
      org.opencontainers.image.source="https://github.com/alixlam/fmristroke"

ENV NO_ET=1 \
    NIPYPE_NO_ET=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1

USER root

WORKDIR /opt/fmristroke

COPY . /opt/fmristroke

RUN python -m pip install --upgrade pip && \
    python -m pip install .

ENTRYPOINT ["fmristroke"]
CMD ["--help"]
