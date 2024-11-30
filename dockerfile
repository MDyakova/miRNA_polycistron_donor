# Load base image
FROM python:3.8 AS builder

# Install Vienna RNA tool
COPY install_viennarna.sh /tmp/
RUN bash /tmp/install_viennarna.sh && rm /tmp/install_viennarna.sh

# Copy files
COPY settings /settings
RUN chmod -R 777 /settings


# Install libraries and prepare runtime
RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -r /settings/requirements.txt && \
    python /settings/load_databases.py

FROM builder AS main

# Copy files
COPY main /main
RUN chmod -R 777 /main
WORKDIR /main

# Launch code
#CMD python miRNA_analysis.py
ENTRYPOINT ["python"]
CMD ["miRNA_analysis.py"]
