FROM condaforge/miniforge3:latest

WORKDIR /app

COPY environment.yml /app/environment.yml
RUN conda env create -f /app/environment.yml \
  && conda clean -a -y

ENV PATH=/opt/conda/envs/chapi/bin:$PATH \
    CHAPI_PYTHON=/opt/conda/envs/chapi/bin/python \
    COOT_PREFIX=/opt/conda/envs/chapi \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

COPY api /app/api
COPY start.sh /app/start.sh

EXPOSE 8000

CMD ["/app/start.sh"]
