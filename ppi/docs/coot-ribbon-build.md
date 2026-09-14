# Coot ribbon build and activation

Roami's Docker image builds the ribbon library from [the Coot fork at commit `49170049247335d59a28ae3c8144f943ce5a8143`](https://github.com/sami-chaaban/coot/tree/49170049247335d59a28ae3c8144f943ce5a8143). This corrects ribbon coverage at chain starts and after physical chain breaks while retaining Coot's mesh API.

[Native Linux CI passed for this exact commit](https://github.com/sami-chaaban/coot/actions/runs/34884156906), including the original-bug reproduction, source build, spline tests and nine patched CHAPI meshes. The local macOS build also passed 18 meshes through Roami's actual worker paths. The complete backend Docker image has not yet been built or deployed; its integration workflow runs when these backend changes are pushed.

`../coot-ribbon-source.json` pins the immutable source commit and archive SHA256. The matching runtime is Bioconda `coot-headless=1.1.20=*_1`, with NumPy 1.26.4 and Biopython 1.86. Creating the conda environment alone installs the original Coot package; the Docker build then compiles and installs the corrected `libcootmoleculestotriangles.so.1.1`.

Run from the `ppi` directory:

```sh
docker build --progress=plain --tag roami-backend:ribbon .
docker run --rm -p 8000:8000 roami-backend:ribbon
```

To activate this change on the hosted backend, rebuild and redeploy its Docker image from the updated `ppi` directory. The site workflow `.github/workflows/roami-backend.yml` validates the image and uploads evidence; it does not publish or deploy it.

The build verifies the source archive, runs the native spline regression and checks nine 9GNQ meshes across all three secondary-structure modes. It checks both chain starts and starts after breaks, including A:47 and K:199. The final image also verifies the checksum and actual loaded path of the ribbon library through CHAPI. Compilers, source checkout, build headers and the test fixture stay in the builder stage.

The image records provenance in `/app/coot-ribbon-build.json` and mesh results in `/app/coot-ribbon-mesh-regression.json`. Repeat the final load check with:

```sh
docker run --rm --entrypoint python roami-backend:ribbon \
  /app/scripts/build_coot_ribbon.py verify \
  --prefix /opt/conda/envs/chapi --manifest /app/coot-ribbon-build.json
```

For a future update, validate the replacement fork commit first, then update both the commit and archive SHA256 in the source lock and rebuild the image. The shared library must remain paired with Coot 1.1.20, build 1; changing that dependency requires rebuilding and validating the integration against its new ABI. Isolated single-residue segments still require a separate visible primitive and are outside this fix.
