#author russell jarvis rjjarvis@asu.edu
#NEURON Dockerfile
#Set the base image to Ubuntu
FROM scidash/scipy-notebook-plus

#Get a whole lot of GNU core development tools
#version control java development, maven
#Libraries required for building MPI from source
#Libraries required for building NEURON from source
#Also DO this part as root.
USER root

RUN apt-get update && apt-get install -y wget bzip2 ca-certificates automake libtool  \
                       python libpython-dev libncurses5-dev libreadline-dev libgsl0-dev cython \
                       cmake ssh
RUN conda install -n python2 pandas numpy scipy matplotlib -y

#Do the rest of the build  as user:
#This will create a more familiar environment to continue developing in.
#with less of a need to chown and chmod everything done as root at dockerbuild completion

USER jovyan
RUN sudo chown -R jovyan /home/jovyan
ENV HOME /home/jovyan
ENV PATH /opt/conda/bin:/opt/conda/bin/conda:/opt/conda/bin/python:$PATH

WORKDIR $HOME
RUN \
  wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.5/nrn-7.5.tar.gz && \
  tar -xzf nrn-7.5.tar.gz && \
  rm nrn-7.5.tar.gz

WORKDIR $HOME/nrn-7.5
RUN ./configure --prefix=`pwd` --with-paranrn --without-iv --with-nrnpython=/opt/conda/envs/python2/bin/python
RUN sudo make all && \
   make install
RUN make all && \
   make install

WORKDIR src/nrnpython
ENV PATH /opt/conda/envs/python2/bin/:$PATH
RUN python2.7 setup.py install
ENV NEURON_HOME $HOME/nrn-7.5/x86_64
ENV PATH $NEURON_HOME/bin:$PATH


# Get shunt level
WORKDIR $HOME
RUN git clone https://github.com/Chris-Currin/shunt-level
WORKDIR shunt-level
