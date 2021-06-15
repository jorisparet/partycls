PROJECT = partycls

.PHONY: all test install coverage clean

all: version

install: version
	python setup.py install

user: version
	python setup.py install --user

test:
	mv partycls partycls_tmp
	python -m unittest discover -s tests; mv partycls_tmp partycls

coverage:
	mv partycls partycls_tmp
	coverage run --source partycls -m unittest discover -s tests; mv partycls_tmp partycls 
	coverage report -m

clean:
	rm -rf ${PROJECT}/*pyc ${PROJECT}/*/*pyc ${PROJECT}/*/*so build
