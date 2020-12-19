# set the environment for compilation
FROM ubuntu:20.04

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
	build-essential \
	tzdata \
	cmake	\
	git \
	wget \
	pkg-config \
	libgomp1

ARG version=1.0.1

# Set the working directory
RUN mkdir -p ./traitrelax/exec \
    && wget https://github.com/halabikeren/TraitRELAX/releases/download/$version/traitrelax_v$version.tar.gz \
    && tar -xf traitrelax_v$version.tar.gz -C ./traitrelax

COPY traitrelax /traitrelax

WORKDIR /traitrelax/

# compile bpp-core
RUN mkdir -p /traitrelax/bpp-core/build/ \
    && cd /traitrelax/bpp-core/build/ \
    && cmake -DCMAKE_INSTALL_PREFIX=/traitrelax/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .. \
    && make \
    && make install

# compile bpp-seq
RUN mkdir -p /traitrelax/bpp-seq/build/ \
    && cd /traitrelax/bpp-seq/build/ \
    && cmake -DCMAKE_INSTALL_PREFIX=/traitrelax/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .. \
    && make \
    && make install

# compile bpp-phyl
RUN mkdir -p /traitrelax/bpp-phyl/build/ \
    && cd /traitrelax/bpp-phyl/build/ \
    && cmake -DCMAKE_INSTALL_PREFIX=/traitrelax/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE OMP_NUM_THREADS=4 .. \
    && make \
    && make install

# compile TraitRELAX
RUN cd /traitrelax/TraitRELAX/ \
    && cmake -DCMAKE_INSTALL_PREFIX=/traitrelax/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE . \
    && make \
    && make install


# execute TraitRELAX with argument given when the docker is run
WORKDIR /traitrelax/exec/
ENTRYPOINT ["/traitrelax/TraitRELAX/TraitRELAX/TraitRELAX"]
CMD ["param=./parameters.bpp"]


