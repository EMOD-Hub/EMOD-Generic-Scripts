FROM ubuntu:24.04
#
RUN apt-get update
RUN apt-get -y install software-properties-common
RUN apt-get update
RUN add-apt-repository universe
RUN apt-get -y upgrade
#
RUN apt-get -y install python3-dev
RUN apt-get -y install libmpich-dev
RUN apt-get -y install libsqlite3-dev
RUN apt-get -y install git
RUN apt-get -y install libc-dev
RUN apt-get -y install g++
#
RUN apt-get -y install scons
#
RUN apt-get clean
#
RUN git clone https://github.com/InstituteforDiseaseModeling/EMOD-Generic
#
RUN git --git-dir=/EMOD-Generic/.git fetch
RUN git --git-dir=/EMOD-Generic/.git checkout Generic-Ongoing
#
ARG USER=DockerComposeUser
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
