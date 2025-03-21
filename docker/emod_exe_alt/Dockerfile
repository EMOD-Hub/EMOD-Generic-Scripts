FROM dtk-centos-buildenv:3.0
#
ARG USER=DockerComposeUser
#
RUN git clone https://github.com/EMOD-Hub/EMOD-Generic
#
RUN git --git-dir=/EMOD-Generic/.git fetch
RUN git --git-dir=/EMOD-Generic/.git checkout Generic-Ongoing
#
RUN scons --Release --jobs=4 --Disease=Generic --directory=/EMOD-Generic/
#
RUN strip EMOD-Generic/build/x64/Release/Eradication/Eradication
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
