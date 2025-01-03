FROM ubuntu:24.04
#
RUN apt-get update
RUN apt-get -y install software-properties-common
RUN apt-get update
RUN add-apt-repository universe
RUN apt-get -y upgrade
#
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN apt-get -y install python3-venv
RUN apt-get -y install mpich
#
RUN apt-get clean
#
RUN python3 -m venv py_env
#
RUN py_env/bin/pip3 install pip             --upgrade
RUN py_env/bin/pip3 install emod-api==2.0.1 --extra-index-url https://packages.idmod.org/api/pypi/pypi-production/simple
