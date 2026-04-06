FROM ubuntu:24.04
#
RUN apt-get update
RUN apt-get upgrade -y
#
RUN apt-get install -y python3-dev
RUN apt-get install -y libmpich-dev
RUN apt-get install -y libsqlite3-dev
RUN apt-get install -y git
RUN apt-get install -y libc-dev
RUN apt-get install -y g++
#
RUN apt-get install -y scons
#
RUN apt-get clean
#
RUN git clone https://github.com/EMOD-Hub/EMOD-Generic
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
RUN outputs/Eradication --get-schema --schema-path outputs/schema.json -P EMOD-Generic/Scripts
#
RUN rm -r EMOD-Generic
