FROM condaforge/miniforge3:latest AS coot-runtime

WORKDIR /app

COPY environment.yml /app/environment.yml
RUN conda env create -f /app/environment.yml \
  && conda clean -a -y

ENV PATH=/opt/conda/envs/chapi/bin:$PATH \
    CHAPI_PYTHON=/opt/conda/envs/chapi/bin/python \
    COOT_PREFIX=/opt/conda/envs/chapi \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

COPY scripts/patch_arpeggio.py /app/scripts/patch_arpeggio.py
RUN python /app/scripts/patch_arpeggio.py --installed
ENV ROAMI_ARPEGGIO_PATCHED=1

# Compile against the runtime's existing libraries. Build tools and headers are
# installed separately so adding a compiler cannot re-solve the chapi environment.
FROM coot-runtime AS coot-ribbon-build
RUN apt-get update \
  && apt-get install -y --no-install-recommends g++ cmake ninja-build \
  && rm -rf /var/lib/apt/lists/* \
  && conda create -y -p /opt/coot-build-headers --override-channels -c conda-forge libboost-headers=1.86.0 \
  && conda clean -a -y

COPY coot-ribbon-source.json /build/coot-ribbon-source.json
COPY scripts/build_coot_ribbon.py /build/build_coot_ribbon.py
COPY data/9gnq.cif /build/9gnq.cif
RUN python /build/build_coot_ribbon.py build \
    --lock /build/coot-ribbon-source.json \
    --prefix /opt/conda/envs/chapi \
    --boost-include /opt/coot-build-headers/include \
    --structure /build/9gnq.cif \
    --output /build/output

# The compiler, source checkout and test structure stay in the build stage.
FROM coot-runtime AS runtime
COPY --from=coot-ribbon-build /build/output/libcootmoleculestotriangles.so.1.1 /opt/conda/envs/chapi/lib/libcootmoleculestotriangles.so.1.1
COPY --from=coot-ribbon-build /build/output/coot-ribbon-build.json /app/coot-ribbon-build.json
COPY --from=coot-ribbon-build /build/output/mesh-regression.json /app/coot-ribbon-mesh-regression.json
COPY --from=coot-ribbon-build /build/output/licenses /app/licenses/coot
COPY scripts/build_coot_ribbon.py /app/scripts/build_coot_ribbon.py
RUN python /app/scripts/build_coot_ribbon.py verify \
    --prefix /opt/conda/envs/chapi \
    --manifest /app/coot-ribbon-build.json

COPY api /app/api
COPY start.sh /app/start.sh

EXPOSE 8000

CMD ["/app/start.sh"]
