FROM continuumio/miniconda3

RUN pip install pandas
RUN pip install gprofiler-official

RUN pip install dist/cdgprofilergenestoterm*whl

ENTRYPOINT ["/opt/miniconda/bin/cdgprofilergenestoterm.py"]
CMD ["--help"]
