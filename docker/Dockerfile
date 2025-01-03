FROM almalinux:9.5
#
RUN dnf upgrade -y
#
RUN dnf -y install python
RUN dnf -y install python-pip
RUN dnf clean all
#
RUN pip install pip          --upgrade
RUN pip install emodpy==2.1  --extra-index-url https://packages.idmod.org/api/pypi/pypi-production/simple