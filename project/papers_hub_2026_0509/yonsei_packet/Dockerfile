FROM python:3.12-slim

LABEL maintainer="Lumenix · Seungho Cook · kukshomr@gmail.com"
LABEL description="Yonsei PDAC drop-in cohort integration · KRAS allele × Moffitt × meta"

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

# Bundled scripts + the existing meta-analysis JSON the new cohort gets pooled into
COPY 21_yonsei_dropin.py /app/
COPY META_ANALYSIS.json /data/pdac_poc/results/method_panel/META_ANALYSIS.json

# Output dir mounted by the user at runtime
RUN mkdir -p /data/pdac_poc/results/yonsei

ENTRYPOINT ["python", "/app/21_yonsei_dropin.py"]
CMD ["--help"]
