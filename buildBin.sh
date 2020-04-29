#! /bin/sh
arch=`uname -m`
version=1.0.0

strip TraitRELAX/traitrelax
tar cvzf traitrelax-${arch}-bin-static-${version}.tar.gz TraitRELAX/traitrelax
