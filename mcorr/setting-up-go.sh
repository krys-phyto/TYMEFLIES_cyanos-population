### Krys Kibler
### Setting up go since mcorr is written in go

# https://github.com/kussell-lab/mcorr # tested with go1.16
# https://go.dev/doc/install/source # go1.16 requires go1.4


# setting up go1.4
mkdir GOROOT_BOOTSTRAP
cd GOROOT_BOOTSTRAP/
git clone https://go.googlesource.com/go goroot
cd goroot/
git checkout release-branch.go1.4
CGO_ENABLED=0
cd src/
./make.bash

# install go1.16 NOT IN THE GOROOT_BOOTSTRAP directory
cd
git clone https://go.googlesource.com/go goroot
cd goroot/
git checkout go1.16
cd src/
export GOROOT_BOOTSTRAP=/home/glbrc.org/kjkibler/GOROOT_BOOTSTRAP/goroot
./all.bash


# install mcorr
export PATH=$PATH:$HOME/goroot/bin:$HOME/.local/bin

go get -u github.com/kussell-lab/mcorr/cmd/mcorr-bam
cd go/pkg/mod/github.com/kussell-lab/mcorr@v0.0.0-20220107174400-4adea90557f1/cmd/mcorr-fit/

python3 setup.py install --user

# error: could not create 'mcorr.egg-info': Permission denied
# GLBRC installed mcorr in an image
