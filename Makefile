PROJECT = partycls

.PHONY: all fortran test install coverage version clean

all: version

install: version
	python setup.py install

user: version
	python setup.py install --user

fortran:
	cd partycls/descriptor/; f2py -c -m realspace_wrap realspace.f90; cd ../../

test: fortran
	python -m unittest discover -s tests

coverage: fortran
	coverage run --source partycls -m unittest discover -s tests 
	coverage report -m

docs:
	pdoc -o docs --force --html partycls

book:
	jupyter-book build ../partycls

clean:
	rm -rf ${PROJECT}/*pyc ${PROJECT}/*/*pyc ${PROJECT}/*/*so build
