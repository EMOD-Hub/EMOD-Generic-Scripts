FROM ubuntu:24.04
#
RUN apt-get update
RUN apt-get upgrade -y
#
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN apt-get install -y python3-venv
RUN apt-get install -y mpich
#
RUN apt-get clean
#
RUN python3 -m venv py_env
#
RUN py_env/bin/pip3 install pip --upgrade
RUN py_env/bin/pip3 install emod-api~=3.1
