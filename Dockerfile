# RootCTrait container image.
#   docker build -t rootctrait:2.4.0 .
#   docker run --rm -v "$PWD":/work rootctrait:2.4.0              # runs run_pipeline.py in /work
#   docker run --rm -v "$PWD":/work rootctrait:2.4.0 python -m tools.qc_rank --results /work/results
# /work must hold params.txt; DATA_ROOT and RESULTS_ROOT are read relative to /work.
FROM python:3.11-slim

LABEL org.opencontainers.image.title="RootCTrait" \
      org.opencontainers.image.source="https://github.com/MA-Chiasson/RootCTrait" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg \
    PYTHONPATH=/opt/rootctrait

WORKDIR /opt/rootctrait
COPY pyproject.toml README.md LICENSE requirements.txt ./
COPY rootctrait ./rootctrait
RUN pip install --no-cache-dir ".[figures,formats]"
COPY run_pipeline.py pipeline_api.py ./
COPY tools ./tools
COPY validation ./validation
COPY tests ./tests
COPY example ./example

WORKDIR /work
CMD ["python", "/opt/rootctrait/run_pipeline.py"]
