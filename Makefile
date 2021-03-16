PROJECT = structural_communities

.PHONY: all test todo install develop doc version clean

all: version

install: version
	python setup.py install

user: version
	python setup.py install --user

clean:
	rm -rf ${PROJECT}/*pyc ${PROJECT}/*/*pyc ${PROJECT}/*/*so build
