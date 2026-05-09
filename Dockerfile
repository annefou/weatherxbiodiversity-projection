FROM mambaorg/micromamba:1.5-jammy

LABEL org.opencontainers.image.source="https://github.com/{{REPO_ORG}}/{{REPO_NAME}}"
LABEL org.opencontainers.image.description="Replication study container for {{REPO_NAME}}"
LABEL org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# Mount any required credentials at runtime, e.g.:
#   docker run -v ~/.cdsapirc:/home/mambauser/.cdsapirc {{REPO_NAME}}
# See data/README.md for per-dataset credential setup.

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["snakemake", "--cores", "1"]
