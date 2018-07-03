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

RUN apt-get update && apt-get install -y wget bzip2 ca-certificates default-jre default-jdk maven automake libtool  \
                       wget python3 libpython3-dev libncurses5-dev libreadline-dev libgsl0-dev cython3 \
                       python3-pip python3-numpy python3-scipy python3-matplotlib python3-mock \
                       ipython3 python3-docutils \
                       python3-mpi4py cmake ssh


#Do the rest of the build  as user:
#This will create a more familiar environment to continue developing in.
#with less of a need to chown and chmod everything done as root at dockerbuild completion

USER jovyan
# Use numpy 1.12.1 until quantities is compatible with 1.13.
RUN conda install -y scipy numpy==1.12.1 matplotlib
RUN sudo chown -R jovyan /home/jovyan
ENV HOME /home/jovyan
ENV PATH /opt/conda/bin:/opt/conda/bin/conda:/opt/conda/bin/python:$PATH

#Test matplotlib
#RUN python -c "import matplotlib"
#Install General MPI, such that mpi4py can later bind with it.

WORKDIR $HOME
RUN \
  wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.5/nrn-7.5.tar.gz && \
  tar -xzf nrn-7.5.tar.gz && \
  rm nrn-7.5.tar.gz

WORKDIR $HOME/nrn-7.5
ENV PATH /usr/bin/python3/python:/opt/conda/bin:/opt/conda/bin/conda:/opt/conda/bin/python:$PATH
RUN ./configure --prefix=`pwd` --with-paranrn --without-iv --with-nrnpython=/opt/conda/bin/python3
RUN sudo make all && \
   make install

RUN make all
RUN make install

WORKDIR src/nrnpython
RUN python setup.py install
ENV NEURON_HOME $HOME/nrn-7.5/x86_64
ENV PATH $NEURON_HOME/bin:$PATH


# Get JNeuroML
#Change to this when PR is accepted
WORKDIR $HOME
RUN git clone https://github.com/NeuroML/jNeuroML
WORKDIR jNeuroML
#RUN python getNeuroML.py
WORKDIR $WORK_HOME
