FROM rockylinux/rockylinux:9.5
#
RUN dnf clean all
RUN dnf upgrade -y
RUN dnf install python-devel -y
RUN dnf install mpich-devel -y
RUN dnf install sqlite-devel -y
RUN dnf install git -y
RUN dnf groupinstall 'Development Tools' -y
RUN dnf clean all
#
RUN python3 -m pip install pip --upgrade
RUN python3 -m pip install scons
#
#
ARG USER=DockerComposeUser
ARG PATH=$PATH:/usr/lib64/mpich/bin
#
RUN git clone https://github.com/EMOD-Hub/EMOD-Generic
#
RUN git --git-dir=/EMOD-Generic/.git fetch
RUN git --git-dir=/EMOD-Generic/.git checkout Generic-Ongoing
#
RUN scons --Release --jobs=4 --Disease=Generic --directory=/EMOD-Generic/
#
RUN mkdir outputs
RUN cp EMOD-Generic/build/x64/Release/Eradication/Eradication outputs/Eradication
#
RUN mkdir outputs/reporter_plugins
RUN cp EMOD-Generic/build/x64/Release/reporter_plugins/*.so outputs/reporter_plugins
#
RUN outputs/Eradication --get-schema --schema-path outputs/schema.json -P EMOD-Generic/Scripts
#
RUN rm -r EMOD-Generic
