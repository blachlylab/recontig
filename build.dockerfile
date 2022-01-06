FROM charesgregory/dlang-htslib-static

ADD . /home/recontig

WORKDIR /home/recontig

RUN dub build --compiler ldc2 -c static-alpine -b release

RUN cp recontig /usr/local/bin

ENTRYPOINT ["/usr/local/bin/recontig"]