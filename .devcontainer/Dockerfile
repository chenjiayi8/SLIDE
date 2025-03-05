FROM ubuntu:22.04 AS compiler
ENV TZ=Europe/Amsterdam
ENV PYTHON_VERSION=3.12.0

# timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install common tools
RUN apt update
RUN apt install -y make build-essential libssl-dev zlib1g-dev \
libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm \
libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev \
liblzma-dev git bash

# remove sh and replace with bash
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Install pyenv
ENV PYENV_ROOT="/pyenv"
ENV PATH "$PYENV_ROOT/bin:$PATH"
ENV eval "$(pyenv init -)"
RUN curl https://pyenv.run | bash

# Install python version
RUN pyenv install $PYTHON_VERSION
ENV PATH="/pyenv/versions/${PYTHON_VERSION}/bin:$PATH"

# Install virtualenv
RUN pip install virtualenv

FROM ubuntu:24.04 AS default
COPY --from=compiler /pyenv /pyenv
ENV TZ=UTC
ENV PYTHON_VERSION=3.12.0
ENV NODE_VERSION=20.11.0
ENV DOCKER_GID=1000
ENV DOCKER_UID=1000

ENV PATH="/pyenv/versions/${PYTHON_VERSION}/bin:$PATH"

# install common tools
RUN apt update
RUN apt install -y bash

# replace shell with bash so we can source files
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update && apt install -y git nano htop wget tree
RUN apt install -y build-essential libssl-dev cmake

# add user ubuntu with same uid and gid
# RUN useradd ubuntu -u $DOCKER_UID -g $DOCKER_GID --create-home -m -s /bin/bash

USER ubuntu
WORKDIR /home/ubuntu

# Create a new directory for ubuntu config
RUN mkdir -p /home/ubuntu/.config

# Create a new virtual environment
RUN virtualenv -p=/pyenv/versions/${PYTHON_VERSION}/bin/python /home/ubuntu/env
RUN source /home/ubuntu/env/bin/activate

# Enable venv
ENV PATH="/home/ubuntu/env/bin:$PATH"


# Install nodejs
ENV NVM_DIR=/home/ubuntu/nvm
RUN mkdir -p $NVM_DIR
RUN git clone https://github.com/nvm-sh/nvm.git "$NVM_DIR" \
    && cd "$NVM_DIR" \
    && git checkout `git describe --abbrev=0 --tags --match "v[0-9]*" $(git rev-list --tags --max-count=1)` \
    && \. "$NVM_DIR/nvm.sh"
# install node and npm
RUN source $NVM_DIR/nvm.sh \
    && nvm install $NODE_VERSION \
    && nvm alias default $NODE_VERSION \
    && nvm use default

# add node and npm to path so the commands are available
ENV NODE_PATH $NVM_DIR/v$NODE_VERSION/lib/node_modules
ENV PATH $NVM_DIR/versions/node/v$NODE_VERSION/bin:$PATH



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


USER ubuntu

# COPY --chown=ubuntu:ubuntu ./requirements.tx[t] /home/ubuntu/requirements.txt
# COPY --chown=ubuntu:ubuntu .devcontaine[r]/requirements.tx[t] /home/ubuntu/requirements.txt

# RUN source /home/ubuntu/env/bin/activate && pip install -r /home/ubuntu/requirements.txt

ENTRYPOINT ["sleep", "infinity"]