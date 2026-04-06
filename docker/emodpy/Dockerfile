FROM ubuntu:24.04
#
RUN apt-get update
RUN apt-get -y upgrade
#
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN apt-get -y install python3-venv
#
RUN apt-get clean
#
ENV PY_VENV=/py_env
RUN python3 -m venv $PY_VENV
ENV PATH="$PY_VENV/bin:$PATH"
#
RUN pip3 install pip --upgrade
RUN pip3 install emodpy~=3.1
