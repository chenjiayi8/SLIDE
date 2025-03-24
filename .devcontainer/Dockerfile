FROM ubuntu:24.04 AS slide
ENV TZ=UTC
ENV PYTHON_VERSION=3.12.0
ENV DOCKER_GID=1000
ENV DOCKER_UID=1000


# install common tools
RUN apt update
RUN apt install -y bash

# replace shell with bash so we can source files
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update && apt install -y git nano htop wget tree
RUN apt install -y build-essential libssl-dev cmake

RUN apt install -y python3 python3-pip python3-venv

USER ubuntu
WORKDIR /home/ubuntu

# Create a new directory for ubuntu config
RUN mkdir -p /home/ubuntu/.config

# Create a new virtual environment
RUN python3 -m venv /home/ubuntu/env
RUN source /home/ubuntu/env/bin/activate

# Enable venv
ENV PATH="/home/ubuntu/env/bin:$PATH"



USER root

# install common tools
RUN apt update
RUN apt install -y sudo
RUN apt install -y nmap net-tools iputils-ping

# add user ubuntu to sudoers
RUN echo "ubuntu ALL=(ALL) NOPASSWD: ALL" | tee /etc/sudoers.d/ubuntu
RUN chmod 0440 /etc/sudoers.d/ubuntu


# add helpful aliases
RUN echo "alias ll='ls -alFh'" >> /home/ubuntu/.bashrc
RUN echo "alias show-tree='tree -I \".git|__pycache__|*.pyc|*.pyo|*.so|*.egg-info|*.egg|*.dist-info\"'" >> /home/ubuntu/.bashrc
RUN echo "alias isort-notebook='nbqa isort --float-to-top '" >> /home/ubuntu/.bashrc
RUN echo "alias isort-notebook='nbqa isort --float-to-top '" >> /home/ubuntu/.bashrc
RUN echo "alias nbrun='jupyter nbconvert --to notebook --execute '" >> /home/ubuntu/.bashrc
RUN echo "source /home/ubuntu/env/bin/activate"  >> /home/ubuntu/.bashrc

USER ubuntu

# COPY --chown=ubuntu:ubuntu ./requirements.tx[t] /home/ubuntu/requirements.txt
# COPY --chown=ubuntu:ubuntu .devcontaine[r]/requirements.tx[t] /home/ubuntu/requirements.txt

# RUN source /home/ubuntu/env/bin/activate && pip install -r /home/ubuntu/requirements.txt

ENTRYPOINT ["sleep", "infinity"]