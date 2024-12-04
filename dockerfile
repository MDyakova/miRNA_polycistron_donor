# Load base image
FROM python:3.8 AS builder

# Install libraries
COPY settings /settings
RUN chmod -R 777 /settings
RUN bash /settings/install_viennarna.sh && rm /settings/install_viennarna.sh

RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -r /settings/requirements.txt && \
    python /settings/load_databases.py

FROM builder AS main

# Copy files
COPY main /main
RUN chmod -R 777 /main
WORKDIR /main

# Run unit-tests
RUN python3 -m pytest unit_tests.py

# Launch code
ENTRYPOINT ["python"]
CMD ["mirna_analysis.py"]
