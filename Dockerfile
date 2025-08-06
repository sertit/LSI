FROM sertit/eo-containers:geo_latest

LABEL org.opencontainers.image.authors="Bastien CORIAT"

RUN apt-get update && apt-get install -y wget libgdal-dev build-essential
RUN pip install gdal==3.6.2
RUN wget https://f-tep.com/sites/default/files/pythonapi/ftep_api-0.0.16-py3-none-any.whl \
    && wget https://f-tep.com/sites/default/files/ftep_util-0.0.8-py3-none-any.whl \
    && PYTHONPATH=/usr/lib/python3/dist-packages/ pip install ftep_api-0.0.16-py3-none-any.whl ftep_util-0.0.8-py3-none-any.whl \
    && rm -f ftep_api-0.0.16-py3-none-any.whl ftep_util-0.0.8-py3-none-any.whl

RUN pip install --index-url https://git.unistra.fr/api/v4/projects/40355/packages/pypi/simple lsi
RUN pip install sertit==1.46.4
RUN pip install numpy==1.26.4
RUN pip install xarray-spatial
