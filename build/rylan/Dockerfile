FROM ubuntu:bionic

RUN mkdir /rylan

RUN apt-get update
RUN apt-get install -y software-properties-common vim

RUN apt-get install -y build-essential python3 python3-pip python3-dev python3-venv
RUN apt-get install -y git curl nano
RUN apt-get update

COPY ./build/rylan/bins/usr/local/lib/* /usr/local/lib/
COPY ./build/rylan/bins/usr/local/bin/* /usr/local/bin/
RUN ln -s /usr/local/lib/libhts.so /usr/local/lib/libhts.so.2

RUN mkdir -p /rylan/dtmtool
RUN mkdir -p /rylan/querytools
run mkdir -p /rylan/databases
RUN mkdir -p /rylan/tmp

COPY ./dtmtool/*.py /rylan/dtmtool/
COPY ./dtmtool/*.sh /rylan/dtmtool/
COPY ./dtmtool/requirements.txt /rylan/dtmtool/
COPY ./databases/* /rylan/databases/
COPY ./queryTools/* /rylan/querytools/

COPY ./build/rylan/deploy/entrypointv2.sh /rylan/
RUN chmod +x /rylan/entrypointv2.sh
RUN pip3 install --upgrade pip
RUN pip3 install -r /rylan/dtmtool/requirements.txt

CMD ["/rylan/entrypointv2.sh"]
